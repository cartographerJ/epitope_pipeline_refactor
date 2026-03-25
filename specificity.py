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
    MAX_NONSPECIFIC_PER_600A2,  # Deprecated, kept for backward compatibility
    MERGE_DISTANCE_THRESHOLD_A,
    RETRY_BLAST,
    SPECIFICITY_IDENTITY_THRESHOLD,
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
    rejected_patches: list              # [(patch, hit_protein_accession, identity), ...]
    unscreened_patches: list            # Patches where BLAST failed/timed out
    blast_results: dict                 # {patch_id: [hits]} for audit trail
    residue_specificity: dict = field(default_factory=dict)
    # {resnum (1-based): bool or None — True=specific, False=non-specific
    #  (in a >=600A² ectodomain patch with >=95% off-target), None=not assessed}
    full_blast_hits: list = field(default_factory=list)
    # All HSPs from _blast_full_sequence() — one dict per HSP
    residue_identity_scores: dict = field(default_factory=dict)
    # {resnum (1-based): float} — max off-target HSP identity covering each position


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
    residue_identity_scores = {}
    try:
        full_hits, residue_identity_scores = _blast_full_sequence(
            target.sequence, target.uniprot_id,
            exclude_gene_name=target.gene_name,
        )
        logger.info("  Full-sequence BLAST: %d off-target hits", len(full_hits))
    except Exception as e:
        logger.warning("  Full-sequence BLAST failed: %s — falling back to per-patch", e)

    # Binary specificity from ectodomain 3D patches (FIXED: uses per-residue identity scores)
    if ectodomain_patches and residue_identity_scores:
        residue_specificity = _screen_patches_binary(
            ectodomain_patches, residue_identity_scores, len(target.sequence),
        )
        n_nonspec = sum(1 for v in residue_specificity.values() if v is False)
        n_assessed = sum(1 for v in residue_specificity.values() if v is not None)
        logger.info(
            "  Ectodomain patch screen: %d/%d assessed residues non-specific (>=%.0f%% off-target)",
            n_nonspec, n_assessed, SPECIFICITY_IDENTITY_THRESHOLD * 100,
        )
    else:
        residue_specificity = {i: None for i in range(1, len(target.sequence) + 1)}

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
            residue_identity_scores=residue_identity_scores,
        )

    specific_patches = []
    rejected_patches = []
    unscreened_patches = []
    all_blast_results = {}

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

        # Whole-patch evaluation (threshold scales with patch size)
        passes, nonspecific_pct, effective_thresh = evaluate_patch_specificity(patch, residue_specificity)

        if passes:
            specific_patches.append(patch)
            logger.info(
                "    Patch %d: PASSED (%.1f%% non-specific <= %.1f%% threshold, %d residues)",
                patch.patch_id, nonspecific_pct, effective_thresh, len(patch.residue_numbers),
            )
        else:
            offending = _find_best_hit(patch_hits, patch.residue_numbers)
            rejected_patches.append((
                patch,
                offending["accession"],
                offending["identity"],
            ))
            logger.info(
                "    Patch %d: REJECTED (%.1f%% non-specific > %.1f%% threshold, %d residues) — "
                "%.1f%% identity with %s",
                patch.patch_id, nonspecific_pct, effective_thresh, len(patch.residue_numbers),
                offending["identity"] * 100, offending["accession"],
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
        residue_identity_scores=residue_identity_scores,
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
    BLAST the full protein sequence and compute per-residue specificity.

    Walks through each HSP alignment character-by-character to determine
    per-position identity with off-target proteins.

    Args:
        full_sequence: Full protein amino acid sequence.
        exclude_accession: UniProt accession to exclude (self).
        exclude_gene_name: Gene name to exclude (e.g. "CEACAM5") for
            RefSeq databases where UniProt IDs aren't in the title.

    Returns:
        Tuple of:
          - hits: List of hit summary dicts (same format as _blast_sequence)
          - residue_specificity: Dict {resnum: float} — max off-target HSP identity
            (0.0 = no coverage, higher = more cross-reactivity risk)
    """
    seq_hash = hashlib.md5(full_sequence.encode()).hexdigest()
    cache_path = CACHE_DIR / "blast" / "fullseq_{}.json".format(seq_hash)

    if cache_path.exists():
        logger.debug("  Full-sequence BLAST cache hit")
        with open(cache_path) as f:
            cached = json.load(f)
        hits = cached.get("hits", [])
        residue_spec = {int(k): v for k, v in cached.get("residue_specificity", {}).items()}
        return hits, residue_spec

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

    # Per-residue: track max off-target HSP identity covering each position.
    # For each position, store the highest identity of any off-target HSP
    # that spans it. Only consider HSPs with >=40% identity and >=30 aa
    # alignment — below that is fold-level similarity, not cross-reactivity.
    max_identity_at = {}  # {resnum: max_hsp_identity}

    hits = []
    for alignment in record.alignments:
        title_upper = alignment.title.upper()

        # Post-filter: only human hits (swissprot titles include organism)
        if "HOMO SAPIENS" not in title_upper and "HUMAN" not in title_upper:
            continue

        # Skip self-hits by accession or gene name
        is_self = False
        if exclude_accession and exclude_accession.upper() in title_upper:
            is_self = True
        if not is_self and exclude_gene_name:
            import re
            pattern = r'\b{}\b'.format(re.escape(exclude_gene_name.upper()))
            if re.search(pattern, title_upper):
                is_self = True
        if is_self:
            continue

        for hsp in alignment.hsps:
            identity = hsp.identities / hsp.align_length if hsp.align_length > 0 else 0.0

            hits.append({
                "accession": alignment.accession,
                "title": alignment.title[:200],
                "identity": identity,
                "identities": hsp.identities,
                "align_length": hsp.align_length,
                "query_start": hsp.query_start,
                "query_end": hsp.query_end,
                "evalue": hsp.expect,
            })

            # Only use meaningful alignments for per-residue mapping
            # >=40% identity filters out structural fold noise (Ig-fold etc.)
            # >=30 aa avoids spurious short matches
            if identity < 0.40 or hsp.align_length < 30:
                continue

            # Mark every query position spanned by this HSP with its identity
            query_pos = hsp.query_start  # 1-based
            for q_char in hsp.query:
                if q_char != "-":
                    if identity > max_identity_at.get(query_pos, 0.0):
                        max_identity_at[query_pos] = identity
                    query_pos += 1

    # Build per-residue map: 0.0 = no off-target coverage, float = max identity
    residue_specificity = {}
    for resnum in range(1, len(full_sequence) + 1):
        residue_specificity[resnum] = max_identity_at.get(resnum, 0.0)

    # Cache — but skip caching if 0 hits for a large protein (likely transient NCBI failure)
    if hits or len(full_sequence) < 100:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump({
                "hits": hits,
                "residue_specificity": {str(k): v for k, v in residue_specificity.items()},
            }, f, indent=2)
    else:
        logger.warning("  BLAST returned 0 hits for %d aa protein — not caching (likely transient)",
                        len(full_sequence))

    return hits, residue_specificity


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
# Binary specificity screening (ectodomain 3D patches)
# ---------------------------------------------------------------------------

def _screen_patches_binary(ectodomain_patches, residue_identity_scores, seq_len):
    """
    Screen ectodomain patches using per-residue BLAST identity scores.

    FIXED: Uses actual per-residue identity scores from BLAST, not HSP ranges.
    A residue is "non-specific" if its identity to any human protein is
    >= SPECIFICITY_IDENTITY_THRESHOLD.

    Args:
        ectodomain_patches: List of SurfacePatch objects (all exposed ectodomain).
        residue_identity_scores: Dict {resnum: float} — max BLAST identity per residue.
        seq_len: Total sequence length.

    Returns:
        Dict {resnum: bool or None} — True=specific, False=non-specific,
        None=not in any assessed patch.
    """
    # Collect all assessed residues (those in any ectodomain patch)
    assessed = set()
    for patch in ectodomain_patches:
        assessed.update(patch.residue_numbers)

    # Build per-residue binary map using identity scores
    result = {}
    for i in range(1, seq_len + 1):
        if i in assessed:
            identity = residue_identity_scores.get(i, 0.0)
            # True = specific (identity < threshold)
            # False = non-specific (identity >= threshold)
            result[i] = (identity < SPECIFICITY_IDENTITY_THRESHOLD)
        else:
            result[i] = None  # Not in any patch

    return result


# ---------------------------------------------------------------------------
# Specificity checking (qualifying patches)
# ---------------------------------------------------------------------------

def _check_patch_specificity(blast_hits, patch, residue_specificity,
                             ca_coords=None):
    """
    Check whether a qualifying patch should be rejected based on specificity.

    Uses a sliding window approach (same as cyno conservation). When the
    window threshold is exceeded, trims failing residues and re-clusters
    survivors into sub-patches.

    Args:
        blast_hits: List of hit dicts from BLAST.
        patch: SurfacePatch object to evaluate.
        residue_specificity: Dict {resnum: bool or None} from binary screening.
        ca_coords: Optional dict {resnum: np.array} for sliding window.

    Returns:
        None if patch passes cleanly.
        List of SurfacePatch sub-patches if trimmed (may be empty).
        Dict (offending hit) for direct rejection (no coords / fallback).
    """
    patch_residue_numbers = patch.residue_numbers

    if not blast_hits:
        return None

    patch_set = set(patch_residue_numbers)

    # Check if we have binary screening data for this patch
    assessed_count = sum(1 for r in patch_set
                         if residue_specificity.get(r) is not None)

    if assessed_count > 0:
        n_nonspec = sum(1 for r in patch_set
                        if residue_specificity.get(r) is False)

        if n_nonspec == 0:
            return None

        # Sliding window: check local non-specific density
        if ca_coords:
            # Build a boolean map: True = specific (OK), False = non-specific
            spec_as_conservation = {
                r: (residue_specificity.get(r) is not False)
                for r in patch_residue_numbers
            }
            worst_pos, worst_count = check_local_mismatch_density(
                patch_residue_numbers, spec_as_conservation, ca_coords,
            )
            if worst_count > MAX_NONSPECIFIC_PER_600A2:
                logger.info(
                    "      Specificity sliding window: residue %d has %d "
                    "non-specific in ~600A² (max %d), trimming...",
                    worst_pos, worst_count, MAX_NONSPECIFIC_PER_600A2,
                )
                # Trim failing residues and re-cluster survivors
                failing = identify_failing_residues(
                    patch_residue_numbers, spec_as_conservation, ca_coords,
                    MAX_NONSPECIFIC_PER_600A2,
                )
                survivors = [r for r in patch_residue_numbers
                             if r not in failing]
                logger.info(
                    "      Trimmed %d non-specific residues, %d survivors",
                    len(failing), len(survivors),
                )
                return _recluster_specificity_survivors(
                    survivors, ca_coords, patch,
                )
            return None
        else:
            # Fallback without coords: flat count
            import math
            from epitope_pipeline.config import VHH_FOOTPRINT_MIN_A2
            max_allowed = max(
                MAX_NONSPECIFIC_PER_600A2,
                int(math.ceil(len(patch_set) / 15.0))
                * MAX_NONSPECIFIC_PER_600A2,
            )
            if n_nonspec > max_allowed:
                return _find_best_hit(blast_hits, patch_residue_numbers)
            return None

    # Fallback: HSP-level identity check (when binary data unavailable)
    for hit in blast_hits:
        if hit["identity"] >= SPECIFICITY_IDENTITY_THRESHOLD:
            if hit["align_length"] >= len(patch_set) * 0.3:
                return hit

    return None


def _find_best_hit(blast_hits, patch_residue_numbers):
    """Find the highest-identity BLAST hit overlapping the patch."""
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
        "identity": SPECIFICITY_IDENTITY_THRESHOLD,
        "title": "off-target",
    }


# ---------------------------------------------------------------------------
# Whole-patch evaluation (NEW)
# ---------------------------------------------------------------------------

def evaluate_patch_specificity(patch, residue_specificity):
    """
    Evaluate a patch based on percentage of non-specific residues.

    The threshold scales with patch size (same logic as conservation):

        effective_threshold = min(base_threshold * sqrt(n_residues / 20), 30%)

    Args:
        patch: SurfacePatch object.
        residue_specificity: Dict {resnum: bool or None} —
            True=specific, False=non-specific, None=not assessed.

    Returns:
        Tuple (passes, nonspecific_percent, effective_threshold):
            - passes: True if patch meets criteria.
            - nonspecific_percent: Float percentage (0-100).
            - effective_threshold: Size-adjusted threshold used (0-100).
    """
    import math

    patch_set = set(patch.residue_numbers)
    n_total = len(patch_set)

    # Count assessed residues
    n_assessed = sum(1 for r in patch_set
                     if residue_specificity.get(r) is not None)

    if n_assessed == 0:
        # No BLAST data — pass by default
        return True, 0.0, 0.0

    # Count non-specific residues
    n_nonspecific = sum(1 for r in patch_set
                        if residue_specificity.get(r) is False)

    nonspecific_percent = (n_nonspecific / n_total) * 100.0

    # Scale threshold by patch size: sqrt(n/20) with 30% ceiling
    effective_threshold = min(
        config.MAX_NONSPECIFIC_PERCENT * math.sqrt(n_total / 20.0),
        30.0,
    )
    passes = (nonspecific_percent <= effective_threshold)

    return passes, nonspecific_percent, effective_threshold


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
# DEPRECATED: Sliding window functions (kept for reference)
# ---------------------------------------------------------------------------

def _recluster_specificity_survivors(survivors, ca_coords, original_patch):
    """
    Re-cluster survivors after trimming non-specific residues.

    Args:
        survivors: List of residue numbers that survived trimming.
        ca_coords: Dict {resnum: np.array} Cα coordinates.
        original_patch: The original SurfacePatch.

    Returns:
        List of SurfacePatch objects (may be empty if no sub-clusters >= 600A²).
    """
    from epitope_pipeline.surface import cluster_surface_patches, SurfacePatch
    from epitope_pipeline.config import VHH_FOOTPRINT_MIN_A2
    import numpy as np

    survivors_with_ca = [r for r in survivors if r in ca_coords]
    if len(survivors_with_ca) < 2:
        return []

    clusters = cluster_surface_patches(survivors_with_ca, ca_coords)

    result = []
    for cluster_residues in clusters:
        # Estimate area from the original patch's per-residue SASA
        # We don't have residue_sasa here directly, but we can use the
        # original patch's total_sasa_a2 proportionally
        n_original = len(original_patch.residue_numbers)
        n_cluster = len(cluster_residues)
        est_area = original_patch.total_sasa_a2 * (n_cluster / n_original) if n_original > 0 else 0

        if est_area < VHH_FOOTPRINT_MIN_A2:
            continue

        coords = np.array([ca_coords[r] for r in cluster_residues])
        centroid = np.mean(coords, axis=0)

        if len(coords) > 1:
            from scipy.spatial.distance import pdist
            max_dim = float(np.max(pdist(coords)))
        else:
            max_dim = 0.0

        sp = SurfacePatch(
            patch_id=0,  # Will be assigned by caller
            residue_numbers=sorted(cluster_residues),
            residue_aas={r: original_patch.residue_aas.get(r, "X")
                         for r in cluster_residues},
            total_sasa_a2=est_area,
            centroid=centroid,
            max_dimension_a=max_dim,
            avg_distance_from_membrane=original_patch.avg_distance_from_membrane,
        )
        result.append(sp)

    return result
