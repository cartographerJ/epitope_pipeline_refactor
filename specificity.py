"""
Step 7: Human Proteome Specificity — screen epitope patches for cross-
reactivity with other human proteins using whole-patch evaluation.

Uses NCBI BLASTp (remote API via BioPython) to search the full protein
sequence against the human swissprot database. Each patch is evaluated
as an atomic unit: patches with ≤15% non-specific residues (covered by
≥70% identity off-target HSPs) pass filtering.

After filtering, adjacent patches are merged into larger contiguous
epitope regions.

Results are cached to cache/blast/ for efficient re-runs.
"""

import hashlib
import io
import json
import logging
import os
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from epitope_pipeline import config
from epitope_pipeline.config import (
    BLAST_DATABASE,
    BLAST_DELAY_SECONDS,
    BLAST_EVALUE_CUTOFF,
    BLAST_HITLIST_SIZE,
    BLAST_ORGANISM_FILTER,
    BLAST_WORD_SIZE,
    CACHE_DIR,
    LOCAL_BLAST_DB_PATH,
    MERGE_DISTANCE_THRESHOLD_A,
    RETRY_BLAST,
    USE_LOCAL_BLAST,
)

logger = logging.getLogger(__name__)

# Resolved BLAST DB path — if the configured path contains spaces (which
# blastp can't handle), we create a symlink under /tmp and use that instead.
_resolved_blast_db_path: Optional[str] = None


def _get_blast_db_path() -> str:
    """Return a blast-safe database path, creating a temp symlink if needed."""
    global _resolved_blast_db_path
    if _resolved_blast_db_path is not None:
        return _resolved_blast_db_path

    db_str = str(LOCAL_BLAST_DB_PATH)
    if " " not in db_str:
        _resolved_blast_db_path = db_str
        return _resolved_blast_db_path

    # blastp can't handle spaces in -db path — symlink the parent directory
    link_dir = Path(tempfile.gettempdir()) / "blast_db_link"
    target_dir = LOCAL_BLAST_DB_PATH.parent
    # Recreate symlink if it doesn't point to the right place
    if link_dir.is_symlink():
        if link_dir.resolve() != target_dir.resolve():
            link_dir.unlink()
            link_dir.symlink_to(target_dir)
    elif not link_dir.exists():
        link_dir.symlink_to(target_dir)
    else:
        # Something non-symlink exists at this path — remove and recreate
        import shutil
        shutil.rmtree(str(link_dir), ignore_errors=True)
        link_dir.symlink_to(target_dir)

    _resolved_blast_db_path = str(link_dir / LOCAL_BLAST_DB_PATH.name)
    logger.info("  BLAST DB path contains spaces — using symlink: %s",
                _resolved_blast_db_path)
    return _resolved_blast_db_path


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class SpecificityResult:
    """Specificity analysis for all candidate patches of one target."""
    uniprot_id: str
    gene_name: str
    specific_patches: list              # SurfacePatch objects passing specificity
    rejected_patches: list              # [(patch, worst_paralog_accession, worst_match_fraction), ...]
    unscreened_patches: list            # Patches where BLAST failed/timed out
    blast_results: dict                 # {patch_id: [hits]} for audit trail
    residue_specificity: dict = field(default_factory=dict)
    # {resnum (1-based): bool or None — True=unique (specific), False=at
    #  least one paralog matches here (non-specific), None=not in any
    #  assessed ectodomain patch}. Derived from paralog_matches.
    full_blast_hits: list = field(default_factory=list)
    # All HSPs from _blast_full_sequence() — one dict per HSP
    paralog_matches: dict = field(default_factory=dict)
    # {accession: set[int]} — for every paralog that passed the HSP quality
    # floor (>=40% id, >=30 aa), the set of 1-based target residue positions
    # where that paralog has the same amino acid as the target. This is the
    # raw input to the patch-level per-paralog max-match-fraction rule.
    patch_worst_paralog: dict = field(default_factory=dict)
    # {patch_id: (accession, match_fraction)} — for each patch evaluated,
    # the single worst paralog and its match fraction. match_fraction is in
    # [0, 1].


class SpecificityError(Exception):
    """Specificity analysis failed."""
    pass


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def filter_specificity(target, conservation_result, ectodomain_patches=None,
                       ca_coords=None):
    """
    Screen the full protein against the human proteome and evaluate patches.

    Steps:
      1. BLAST the full protein sequence against human refseq_protein
      2. Build per-residue specificity map from alignment details
      3. For each conserved patch: check overlap with off-target hits -> accept/reject

    Args:
        target: TargetInfo with sequence and uniprot_id.
        conservation_result: ConservationResult with conserved_patches.

    Returns:
        SpecificityResult with per-residue map, filtered patches, and BLAST data.
    """
    # Pre-flight: validate local BLAST database exists if enabled
    if USE_LOCAL_BLAST:
        db_path = Path(LOCAL_BLAST_DB_PATH)
        # Check for .phr file (one of the BLAST database index files)
        if not (db_path.parent / f"{db_path.name}.phr").exists():
            error_msg = (
                f"Local BLAST database not found at {LOCAL_BLAST_DB_PATH}. "
                f"Run 'bash epitope_pipeline/scripts/setup_blast_db.sh' to set up, "
                f"or set USE_LOCAL_BLAST=False in config.py to use remote NCBI BLAST."
            )
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        # Resolve path (creates temp symlink if path contains spaces)
        resolved = _get_blast_db_path()
        logger.debug("  Using local BLAST database: %s", resolved)
    else:
        logger.debug("  Using remote NCBI BLAST API")

    # Step 1: Full-sequence BLAST (single call, cached)
    logger.info("  BLASTing full %s sequence (%d aa) against human proteome...",
                target.gene_name, len(target.sequence))

    full_hits = []
    paralog_matches = {}
    try:
        full_hits, paralog_matches = _blast_full_sequence(
            target.sequence, target.uniprot_id,
            exclude_gene_name=target.gene_name,
        )
        logger.info(
            "  Full-sequence BLAST: %d HSPs, %d paralogs above quality floor",
            len(full_hits), len(paralog_matches),
        )
    except Exception as e:
        logger.warning("  Full-sequence BLAST failed: %s — falling back to per-patch", e)

    # Derive the per-residue binary map (for figures/exports only — the
    # pass/fail decision uses per-paralog fractions, not this map).
    residue_specificity = _derive_residue_specificity(
        paralog_matches, ectodomain_patches or [], len(target.sequence),
    )
    if ectodomain_patches:
        n_nonspec = sum(1 for v in residue_specificity.values() if v is False)
        n_assessed = sum(1 for v in residue_specificity.values() if v is not None)
        logger.info(
            "  Ectodomain surface: %d/%d residues shared with at least one paralog",
            n_nonspec, n_assessed,
        )

    # Step 2: Per-patch evaluation
    if not conservation_result.conserved_patches:
        logger.info("  %s: No conserved patches to screen for specificity", target.gene_name)
        return SpecificityResult(
            uniprot_id=target.uniprot_id,
            gene_name=target.gene_name,
            specific_patches=[],
            rejected_patches=[],
            unscreened_patches=[],
            blast_results={},
            residue_specificity=residue_specificity,
            full_blast_hits=full_hits,
            paralog_matches=paralog_matches,
        )

    specific_patches = []
    rejected_patches = []
    unscreened_patches = []
    all_blast_results = {}
    patch_worst_paralog = {}

    total = len(conservation_result.conserved_patches)
    for idx, patch in enumerate(conservation_result.conserved_patches):
        logger.info(
            "  Screening patch %d/%d (patch_id=%d, %d residues)...",
            idx + 1, total, patch.patch_id, len(patch.residue_numbers),
        )

        # Use full-sequence BLAST hits filtered to patch region
        patch_hits = _filter_hits_to_patch(full_hits, patch.residue_numbers)
        all_blast_results[patch.patch_id] = patch_hits

        # Fallback to per-patch BLAST if full BLAST failed
        if not full_hits:
            segments = _extract_patch_sequences(target.sequence, patch.residue_numbers)
            blast_failed = False
            for seg_idx, (seq_segment, seg_start, seg_end) in enumerate(segments):
                if len(seq_segment) < 10:
                    continue
                logger.info("    Fallback BLAST: segment %d, residues %d-%d (%d aa)",
                            seg_idx, seg_start, seg_end, len(seq_segment))
                try:
                    hits = _blast_sequence(seq_segment, exclude_accession=target.uniprot_id)
                    patch_hits.extend(hits)
                    if not USE_LOCAL_BLAST and (seg_idx < len(segments) - 1 or idx < total - 1):
                        time.sleep(BLAST_DELAY_SECONDS)
                except Exception as e:
                    logger.warning("    BLAST failed for segment %d: %s", seg_idx, e)
                    blast_failed = True
            all_blast_results[patch.patch_id] = patch_hits
            if blast_failed and not patch_hits:
                unscreened_patches.append(patch)
                continue

        # Per-paralog max-match-fraction rule
        passes, worst_pct, threshold, worst_acc = evaluate_patch_specificity(
            patch, paralog_matches)
        patch_worst_paralog[patch.patch_id] = (worst_acc, worst_pct / 100.0)

        if passes:
            specific_patches.append(patch)
            logger.info(
                "    Patch %d: PASSED (worst paralog %s @ %.1f%% <= %.1f%%, %d residues)",
                patch.patch_id, worst_acc or "—", worst_pct, threshold,
                len(patch.residue_numbers),
            )
        else:
            rejected_patches.append((patch, worst_acc, worst_pct / 100.0))
            logger.info(
                "    Patch %d: REJECTED (worst paralog %s @ %.1f%% > %.1f%%, %d residues)",
                patch.patch_id, worst_acc or "—", worst_pct, threshold,
                len(patch.residue_numbers),
            )

    # Step 7.5: Merge adjacent patches
    if specific_patches:
        logger.info(
            "  Merging adjacent patches (threshold=%.1f A)...",
            MERGE_DISTANCE_THRESHOLD_A,
        )

        pre_merge_count = len(specific_patches)
        specific_patches = merge_adjacent_patches(
            specific_patches,
            distance_threshold_a=MERGE_DISTANCE_THRESHOLD_A,
        )
        post_merge_count = len(specific_patches)

        logger.info(
            "    %d patches merged into %d regions",
            pre_merge_count, post_merge_count,
        )

        # Log merged patch details
        for patch in specific_patches:
            logger.info(
                "    Merged patch %d: %d residues, %.0f A²",
                patch.patch_id, len(patch.residue_numbers), patch.total_sasa_a2,
            )

    logger.info(
        "  %s specificity: %d pass, %d rejected, %d unscreened",
        target.gene_name,
        len(specific_patches), len(rejected_patches), len(unscreened_patches),
    )

    return SpecificityResult(
        uniprot_id=target.uniprot_id,
        gene_name=target.gene_name,
        specific_patches=specific_patches,
        rejected_patches=rejected_patches,
        unscreened_patches=unscreened_patches,
        blast_results=all_blast_results,
        residue_specificity=residue_specificity,
        full_blast_hits=full_hits,
        paralog_matches=paralog_matches,
        patch_worst_paralog=patch_worst_paralog,
    )


# ---------------------------------------------------------------------------
# Full-sequence BLAST with per-residue parsing
# ---------------------------------------------------------------------------

def _run_local_blast(sequence: str) -> str:
    """
    Run local BLAST using blastp command-line tool.

    Args:
        sequence: Amino acid sequence to query.

    Returns:
        XML string output from blastp (same format as NCBIWWW.qblast()).

    Raises:
        subprocess.CalledProcessError: If blastp command fails.
        FileNotFoundError: If blastp executable not found.
    """
    # Write query sequence to temporary FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(f">query\n{sequence}\n")
        query_file = f.name

    try:
        # Build blastp command (use resolved path to avoid spaces)
        cmd = [
            "blastp",
            "-db", _get_blast_db_path(),
            "-query", query_file,
            "-evalue", str(BLAST_EVALUE_CUTOFF),
            "-word_size", str(BLAST_WORD_SIZE),
            "-max_target_seqs", str(BLAST_HITLIST_SIZE),
            "-outfmt", "5",  # XML output (same format as NCBIWWW)
        ]

        # Execute BLAST
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=300,  # 5-minute timeout
        )

        return result.stdout  # XML string

    except subprocess.CalledProcessError as e:
        logger.error("Local BLAST failed: %s", e.stderr)
        raise
    except FileNotFoundError:
        logger.error(
            "blastp executable not found. Ensure NCBI BLAST+ is installed "
            "(brew install blast / apt-get install ncbi-blast+)"
        )
        raise
    finally:
        # Clean up temp file
        try:
            os.unlink(query_file)
        except Exception:
            pass


def _blast_full_sequence(full_sequence, exclude_accession, exclude_gene_name=None):
    """
    BLAST the full protein sequence and return per-paralog match sets.

    Walks each qualifying HSP alignment character-by-character and, for
    every paralog, records the set of target residue positions where that
    paralog has the SAME amino acid as the target. The patch-level
    specificity filter then asks, per patch: does any single paralog match
    more than MAX_NONSPECIFIC_PERCENT of the patch residues?

    Args:
        full_sequence: Full protein amino acid sequence.
        exclude_accession: UniProt accession to exclude (self).
        exclude_gene_name: Gene name to exclude (e.g. "CEACAM5") for
            RefSeq databases where UniProt IDs aren't in the title.

    Returns:
        Tuple of:
          - hits: List of hit summary dicts (same format as _blast_sequence)
          - paralog_matches: Dict {accession: set[int]} — for each paralog
            above the HSP quality floor (>=40% id, >=30 aa), the set of
            1-based target residue positions where that paralog matches the
            target amino acid exactly.
    """
    seq_hash = hashlib.md5(full_sequence.encode()).hexdigest()
    cache_path = CACHE_DIR / "blast" / "fullseq_{}.json".format(seq_hash)

    if cache_path.exists():
        logger.debug("  Full-sequence BLAST cache hit")
        with open(cache_path) as f:
            cached = json.load(f)
        # New cache format is keyed on "paralog_matches"; older caches keyed
        # on "residue_specificity" are schema-incompatible and ignored.
        if "paralog_matches" in cached:
            hits = cached.get("hits", [])
            paralog_matches = {
                acc: set(positions)
                for acc, positions in cached["paralog_matches"].items()
            }
            return hits, paralog_matches
        logger.debug("  Cache schema mismatch — recomputing")

    from Bio.Blast import NCBIWWW, NCBIXML

    max_attempts = RETRY_BLAST["attempts"]
    backoff = RETRY_BLAST["backoff_seconds"]

    blast_kwargs = {
        "expect": BLAST_EVALUE_CUTOFF,
        "word_size": BLAST_WORD_SIZE,
        "hitlist_size": BLAST_HITLIST_SIZE,
    }
    if BLAST_ORGANISM_FILTER:
        blast_kwargs["entrez_query"] = BLAST_ORGANISM_FILTER

    for attempt in range(1, max_attempts + 1):
        try:
            if USE_LOCAL_BLAST:
                # Local BLAST: run blastp subprocess
                xml_string = _run_local_blast(full_sequence)
                result_handle = io.StringIO(xml_string)
            else:
                # Remote BLAST: use NCBI API
                result_handle = NCBIWWW.qblast(
                    "blastp", BLAST_DATABASE, full_sequence, **blast_kwargs,
                )

            record = next(NCBIXML.parse(result_handle))
            result_handle.close()
            break

        except Exception as e:
            if attempt < max_attempts:
                wait = backoff * attempt
                logger.warning("  Full BLAST attempt %d/%d failed: %s. Retrying in %ds...",
                               attempt, max_attempts, e, wait)
                time.sleep(wait)
            else:
                raise

    # Per-paralog match sets: {accession: set of 1-based positions where
    # paralog AA matches target AA}. We walk qseq/sseq for every HSP that
    # passes the quality floor (>=40% identity, >=30 aa alignment) to
    # filter out fold-level noise (Ig-fold matches etc.) while keeping
    # every credible paralog on the table for the patch-level decision.
    paralog_matches = {}  # {accession: set[int]}

    import re
    self_gene_pattern = None
    if exclude_gene_name:
        self_gene_pattern = re.compile(
            r'\b{}\b'.format(re.escape(exclude_gene_name.upper())))

    hits = []
    for alignment in record.alignments:
        title_upper = alignment.title.upper()

        # Post-filter: only human hits (swissprot titles include organism)
        if "HOMO SAPIENS" not in title_upper and "HUMAN" not in title_upper:
            continue

        # Skip self-hits by accession or gene name
        if exclude_accession and exclude_accession.upper() in title_upper:
            continue
        if self_gene_pattern and self_gene_pattern.search(title_upper):
            continue

        accession = alignment.accession
        for hsp in alignment.hsps:
            identity = hsp.identities / hsp.align_length if hsp.align_length > 0 else 0.0

            hits.append({
                "accession": accession,
                "title": alignment.title[:200],
                "identity": identity,
                "identities": hsp.identities,
                "align_length": hsp.align_length,
                "query_start": hsp.query_start,
                "query_end": hsp.query_end,
                "evalue": hsp.expect,
            })

            # Quality floor: skip HSPs that are too divergent or too short
            # to meaningfully inform a cross-reactivity decision.
            if identity < 0.40 or hsp.align_length < 30:
                continue

            # Walk aligned pair position-by-position. query_pos is 1-based
            # in the full target sequence. A gap in the query doesn't
            # consume a target position; a gap in the subject is just a
            # non-match at that target position.
            matched = paralog_matches.setdefault(accession, set())
            query_pos = hsp.query_start
            for q_char, s_char in zip(hsp.query, hsp.sbjct):
                if q_char == "-":
                    continue
                if s_char != "-" and q_char.upper() == s_char.upper():
                    matched.add(query_pos)
                query_pos += 1

    # Cache — but skip caching if 0 hits for a large protein (likely
    # transient NCBI failure). Sets serialize as sorted lists.
    if hits or len(full_sequence) < 100:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump({
                "hits": hits,
                "paralog_matches": {
                    acc: sorted(positions)
                    for acc, positions in paralog_matches.items()
                },
            }, f, indent=2)
    else:
        logger.warning("  BLAST returned 0 hits for %d aa protein — not caching (likely transient)",
                        len(full_sequence))

    return hits, paralog_matches


def _filter_hits_to_patch(full_hits, patch_residue_numbers):
    """
    Filter full-sequence BLAST hits to those overlapping a given patch.

    Args:
        full_hits: List of hit dicts from full-sequence BLAST.
        patch_residue_numbers: List of 1-based residue numbers.

    Returns:
        List of hit dicts whose query alignment overlaps the patch.
    """
    if not full_hits or not patch_residue_numbers:
        return []

    patch_min = min(patch_residue_numbers)
    patch_max = max(patch_residue_numbers)

    overlapping = []
    for hit in full_hits:
        q_start = hit.get("query_start", 0)
        q_end = hit.get("query_end", 0)
        # Check if the hit alignment overlaps the patch region
        if q_start <= patch_max and q_end >= patch_min:
            overlapping.append(hit)

    return overlapping


# ---------------------------------------------------------------------------
# Sequence extraction
# ---------------------------------------------------------------------------

def _extract_patch_sequences(full_sequence, residue_numbers):
    """
    Extract contiguous sequence segments that span the patch residues.

    Surface patches may not be contiguous in primary sequence. This function
    groups patch residues into contiguous runs (allowing small gaps <= 20
    residues) and extracts the sequence for each run.

    Args:
        full_sequence: Full protein sequence (0-based indexing internally,
                       but residue_numbers are 1-based).
        residue_numbers: List of 1-based residue numbers in the patch.

    Returns:
        List of (sequence_segment, start_resnum, end_resnum) tuples.
    """
    if not residue_numbers:
        return []

    sorted_res = sorted(residue_numbers)

    # Group into segments (split at gaps > 20 residues)
    segments = []
    seg_start = sorted_res[0]
    seg_end = sorted_res[0]

    for resnum in sorted_res[1:]:
        if resnum - seg_end > 20:
            # New segment
            segments.append((seg_start, seg_end))
            seg_start = resnum
        seg_end = resnum

    segments.append((seg_start, seg_end))

    # Extract sequences (convert 1-based to 0-based for slicing)
    result = []
    for start, end in segments:
        seq = full_sequence[start - 1:end]  # 1-based to 0-based
        result.append((seq, start, end))

    return result


# ---------------------------------------------------------------------------
# BLAST
# ---------------------------------------------------------------------------

def _blast_sequence(query_sequence, exclude_accession=None):
    """
    BLASTp a query sequence against human refseq_protein using NCBI remote API.

    Uses BioPython's NCBIWWW.qblast() which handles the NCBI polling
    internally. Results are cached by sequence hash.

    Args:
        query_sequence: Amino acid sequence string to search.
        exclude_accession: UniProt/RefSeq accession to exclude (self).

    Returns:
        List of hit dicts: {accession, title, identity, evalue,
                            query_start, query_end, align_length}
    """
    # Check cache
    seq_hash = hashlib.md5(query_sequence.encode()).hexdigest()
    cache_path = CACHE_DIR / "blast" / "{}.json".format(seq_hash)

    if cache_path.exists():
        logger.debug("    BLAST cache hit: %s", cache_path.name)
        with open(cache_path) as f:
            return json.load(f)

    from Bio.Blast import NCBIWWW, NCBIXML

    max_attempts = RETRY_BLAST["attempts"]
    backoff = RETRY_BLAST["backoff_seconds"]

    blast_kwargs = {
        "expect": BLAST_EVALUE_CUTOFF,
        "word_size": BLAST_WORD_SIZE,
        "hitlist_size": BLAST_HITLIST_SIZE,
    }
    if BLAST_ORGANISM_FILTER:
        blast_kwargs["entrez_query"] = BLAST_ORGANISM_FILTER

    for attempt in range(1, max_attempts + 1):
        try:
            if USE_LOCAL_BLAST:
                # Local BLAST: run blastp subprocess
                xml_string = _run_local_blast(query_sequence)
                result_handle = io.StringIO(xml_string)
            else:
                # Remote BLAST: use NCBI API
                result_handle = NCBIWWW.qblast(
                    "blastp", BLAST_DATABASE, query_sequence, **blast_kwargs,
                )

            blast_records = NCBIXML.parse(result_handle)
            record = next(blast_records)

            hits = []
            for alignment in record.alignments:
                title_upper = alignment.title.upper()
                # Post-filter: only human hits
                if "HOMO SAPIENS" not in title_upper and "HUMAN" not in title_upper:
                    continue
                # Skip self-hits
                accession = alignment.accession
                if exclude_accession and exclude_accession.upper() in title_upper:
                    continue

                for hsp in alignment.hsps:
                    identity = hsp.identities / hsp.align_length if hsp.align_length > 0 else 0.0
                    hits.append({
                        "accession": accession,
                        "title": alignment.title[:200],
                        "identity": identity,
                        "identities": hsp.identities,
                        "align_length": hsp.align_length,
                        "query_start": hsp.query_start,
                        "query_end": hsp.query_end,
                        "evalue": hsp.expect,
                    })

            result_handle.close()

            # Cache results
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_path, "w") as f:
                json.dump(hits, f, indent=2)

            logger.debug("    BLAST returned %d hits", len(hits))
            return hits

        except Exception as e:
            if attempt < max_attempts:
                wait = backoff * attempt
                logger.warning(
                    "    BLAST attempt %d/%d failed: %s. Retrying in %ds...",
                    attempt, max_attempts, e, wait,
                )
                time.sleep(wait)
            else:
                raise


# ---------------------------------------------------------------------------
# Per-paralog patch specificity
# ---------------------------------------------------------------------------

def evaluate_patch_specificity(patch, paralog_matches):
    """
    Evaluate a patch by per-paralog match fraction.

    For each paralog that passed the HSP quality floor in BLAST, compute
    what fraction of the patch residues it matches exactly. The patch
    fails if any single paralog matches more than MAX_NONSPECIFIC_PERCENT
    of the patch — because an antibody raised to that patch surface would
    cross-react with the matching paralog.

    Rationale: cross-reactivity is a per-paralog problem. An antibody
    cross-reacts with a specific off-target protein, not a chimera of all
    paralogs. Taking the max over paralogs gives each paralog an
    independent veto while averaging across paralogs would let one
    near-identical off-target hide behind several distant ones.

    Args:
        patch: SurfacePatch object.
        paralog_matches: Dict {accession: set[int]} — for each paralog,
            the set of target residues it matches exactly. From
            _blast_full_sequence().

    Returns:
        Tuple (passes, worst_match_percent, threshold, worst_paralog):
            - passes: True if max-over-paralogs match fraction <= threshold.
            - worst_match_percent: Float 0-100 — worst paralog's match %.
            - threshold: MAX_NONSPECIFIC_PERCENT (0-100).
            - worst_paralog: Accession of the worst paralog (or None).
    """
    patch_set = set(patch.residue_numbers)
    n_total = len(patch_set)
    threshold = config.MAX_NONSPECIFIC_PERCENT

    if n_total == 0 or not paralog_matches:
        return True, 0.0, threshold, None

    worst_acc = None
    worst_frac = 0.0
    for accession, matched_positions in paralog_matches.items():
        frac = len(patch_set & matched_positions) / n_total
        if frac > worst_frac:
            worst_frac = frac
            worst_acc = accession

    worst_percent = worst_frac * 100.0
    passes = worst_percent <= threshold
    return passes, worst_percent, threshold, worst_acc


def _derive_residue_specificity(paralog_matches, ectodomain_patches, seq_len):
    """
    Build the per-residue binary specificity map used by figures/exports.

    A residue is marked:
      - True  (specific)     — in an assessed patch, no paralog matches here
      - False (non-specific) — in an assessed patch, at least one paralog matches
      - None  (not assessed) — outside any ectodomain surface patch

    This is purely for visualization/reporting — the pass/fail decision
    comes from evaluate_patch_specificity's per-paralog max fraction,
    not from aggregating these flags.
    """
    assessed = set()
    for patch in ectodomain_patches or []:
        assessed.update(patch.residue_numbers)

    any_paralog_matches = set()
    for positions in paralog_matches.values():
        any_paralog_matches.update(positions)

    result = {}
    for i in range(1, seq_len + 1):
        if i not in assessed:
            result[i] = None
        else:
            result[i] = i not in any_paralog_matches
    return result


def _find_best_hit(blast_hits, patch_residue_numbers):
    """Find the highest-identity BLAST hit overlapping the patch.

    Used only for the rejection-message annotation — tells the user which
    off-target hit most clearly explains why the patch failed.
    """
    patch_set = set(patch_residue_numbers)
    best_hit = None
    best_ident = 0.0
    for hit in blast_hits:
        if hit["identity"] > best_ident:
            q_s = hit.get("query_start", 0)
            q_e = hit.get("query_end", 0)
            if any(q_s <= r <= q_e for r in patch_set):
                best_ident = hit["identity"]
                best_hit = hit
    return best_hit or {
        "accession": "unknown",
        "identity": 0.0,
        "title": "off-target",
    }


def merge_adjacent_patches(patches, distance_threshold_a=15.0):
    """
    Merge spatially adjacent patches into larger contiguous regions.

    Uses centroid-based distance to identify adjacent patches, then merges
    them into larger epitope regions. Patches within distance_threshold_a
    of each other are merged.

    Args:
        patches: List of SurfacePatch objects.
        distance_threshold_a: Distance threshold in Angstroms (default 15.0).

    Returns:
        List of merged SurfacePatch objects.
    """
    import numpy as np
    from scipy.spatial.distance import pdist, squareform
    from .surface import SurfacePatch

    if len(patches) == 0:
        return []

    if len(patches) == 1:
        return patches

    # Calculate pairwise centroid distances
    centroids = np.array([p.centroid for p in patches])
    dist_matrix = squareform(pdist(centroids))

    # Single-linkage clustering: merge patches within threshold
    merged_groups = []
    assigned = set()

    for i, patch_i in enumerate(patches):
        if i in assigned:
            continue

        # Start new group
        group = [i]
        assigned.add(i)

        # Find all patches within threshold (transitively)
        to_check = [i]
        while to_check:
            current = to_check.pop(0)
            for j, patch_j in enumerate(patches):
                if j in assigned:
                    continue
                if dist_matrix[current, j] <= distance_threshold_a:
                    group.append(j)
                    assigned.add(j)
                    to_check.append(j)

        merged_groups.append(group)

    # Create merged patches
    merged_patches = []
    for group_idx, group in enumerate(merged_groups):
        if len(group) == 1:
            # Single patch, keep as-is
            merged_patches.append(patches[group[0]])
        else:
            # Merge multiple patches
            group_patches = [patches[i] for i in group]

            # Combine residue numbers
            all_residues = []
            for p in group_patches:
                all_residues.extend(p.residue_numbers)
            merged_residues = sorted(set(all_residues))

            # Combine residue AAs
            merged_aas = {}
            for p in group_patches:
                merged_aas.update(p.residue_aas)

            # Sum surface areas
            total_sasa = sum(p.total_sasa_a2 for p in group_patches)

            # Recalculate centroid (weighted by SASA)
            weighted_centroid = sum(
                p.centroid * p.total_sasa_a2 for p in group_patches
            ) / total_sasa

            # Recalculate max dimension
            coords = np.array([p.centroid for p in group_patches])
            max_dist = pdist(coords).max() if len(coords) > 1 else 0.0

            # Average distance from membrane
            avg_membrane_dist = np.mean([
                p.avg_distance_from_membrane for p in group_patches
            ])

            # Create merged patch (use first patch's ID as base)
            merged_patch = SurfacePatch(
                patch_id=patches[group[0]].patch_id,
                residue_numbers=merged_residues,
                residue_aas=merged_aas,
                total_sasa_a2=total_sasa,
                centroid=weighted_centroid,
                max_dimension_a=max_dist,
                avg_distance_from_membrane=avg_membrane_dist,
            )
            merged_patches.append(merged_patch)

    return merged_patches


# ---------------------------------------------------------------------------
