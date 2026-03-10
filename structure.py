"""
Step 2: Structure Acquisition — obtain 3D structures for target proteins.

Prefers validated experimental structures from the PDB (X-ray or cryo-EM,
<=3.5A resolution). Falls back to AlphaFold DB pre-computed structures for
proteins with low experimental coverage, and to Tamarind AlphaFold prediction
as a last resort.

Structures are ranked by coverage of the full-length protein weighted by
resolution. The best structure is downloaded to the run's structures/ directory.
"""

import io
import json
import logging
import os
import sys
import time
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import requests

from epitope_pipeline.config import (
    ALLOWED_METHODS,
    ALPHAFOLD_DB_URL,
    CACHE_DIR,
    FOLDING_NUM_MODELS,
    FOLDING_NUM_RECYCLES,
    FOLDING_TOOL,
    FOLDING_TIMEOUT,
    MIN_COVERAGE_FRACTION,
    MIN_ECTODOMAIN_COVERAGE,
    RAW_JOBS_DIR,
    RCSB_API,
    RCSB_FILES,
    RCSB_SEARCH_API,
    RESOLUTION_THRESHOLD_A,
    RETRY_RCSB,
)

# Import TamarindClient from sibling project
sys.path.insert(0, str(Path(__file__).parent.parent))
from tamarind.client import TamarindClient, TamarindError

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class StructureResult:
    """Result of structure acquisition for one target."""
    uniprot_id: str
    gene_name: str
    pdb_path: str               # Local path to downloaded/predicted PDB file
    source: str                 # "pdb_experimental" or "tamarind_alphafold"
    pdb_id: str                 # PDB ID or Tamarind job name
    method: str                 # "X-RAY DIFFRACTION", "ELECTRON MICROSCOPY", or "AlphaFold"
    resolution: Optional[float] # Resolution in Angstroms (None for predicted)
    chain_id: str               # Chain ID of the target protein in the PDB file
    coverage_start: int         # First residue covered (1-based)
    coverage_end: int           # Last residue covered (1-based)


class StructureAcquisitionError(Exception):
    """No suitable structure found or folding failed."""
    pass


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def acquire_structure(target, structures_dir, force_experimental=False):
    """
    Obtain a 3D structure for a single target protein.

    Default: AlphaFold DB (full-length, gap-free coverage).
    Opt-in: experimental PDB (X-ray/cryo-EM, >=80% ectodomain coverage).
    Last resort: Tamarind AlphaFold prediction (de novo).

    Args:
        target: TargetInfo object from target_input module.
        structures_dir: Directory to save downloaded PDB files.
        force_experimental: If True, prefer experimental PDB over AlphaFold
            when available with sufficient ectodomain coverage.

    Returns:
        StructureResult with path to the PDB file and metadata.

    Raises:
        StructureAcquisitionError: If no structure can be obtained.
    """
    structures_dir = Path(structures_dir)
    structures_dir.mkdir(parents=True, exist_ok=True)

    # Experimental PDB path (opt-in only)
    if force_experimental and target.pdb_ids:
        logger.info("  Evaluating %d PDB entries for %s...", len(target.pdb_ids), target.gene_name)
        ranked = _rank_pdb_structures(target.pdb_ids, target)
        if ranked:
            best = ranked[0]
            coverage_frac = best.get("coverage_fraction", 0.0)
            logger.info(
                "  Best PDB: %s (%s, %.1fA, chain %s, residues %d-%d, %.0f%% ectodomain coverage)",
                best["pdb_id"], best["method"], best["resolution"],
                best["chain_id"], best["coverage_start"], best["coverage_end"],
                coverage_frac * 100,
            )

            # If ectodomain coverage is too low, prefer AlphaFold
            if coverage_frac < MIN_ECTODOMAIN_COVERAGE:
                logger.info(
                    "  Ectodomain coverage %.0f%% < %.0f%% threshold — deferring to AlphaFold",
                    coverage_frac * 100, MIN_ECTODOMAIN_COVERAGE * 100,
                )
            else:
                pdb_path = _download_pdb(best["pdb_id"], structures_dir)
                return StructureResult(
                    uniprot_id=target.uniprot_id,
                    gene_name=target.gene_name,
                    pdb_path=str(pdb_path),
                    source="pdb_experimental",
                    pdb_id=best["pdb_id"],
                    method=best["method"],
                    resolution=best["resolution"],
                    chain_id=best["chain_id"],
                    coverage_start=best["coverage_start"],
                    coverage_end=best["coverage_end"],
                )
        else:
            logger.info("  No PDB entries pass resolution/method filters for %s", target.gene_name)

    # Default: AlphaFold DB pre-computed structure (most human proteins)
    afdb_result = _try_alphafold_db(target, structures_dir)
    if afdb_result is not None:
        return afdb_result

    # Fallback: Tamarind AlphaFold prediction (de novo)
    logger.info("  Folding %s via Tamarind %s...", target.gene_name, FOLDING_TOOL)
    return _fold_with_tamarind(target, structures_dir)


def acquire_structures(targets, structures_dir):
    """
    Acquire structures for multiple targets.

    Args:
        targets: List of TargetInfo objects.
        structures_dir: Directory to save PDB files.

    Returns:
        Dict mapping uniprot_id -> StructureResult.
        Targets that fail are logged and skipped.
    """
    results = {}
    for target in targets:
        try:
            result = acquire_structure(target, structures_dir)
            results[target.uniprot_id] = result
        except StructureAcquisitionError as e:
            logger.error("  Structure acquisition failed for %s: %s", target.gene_name, e)
    return results


def download_alphafold_supplement(target, structures_dir):
    """
    Download AlphaFold DB structure as a supplement for experimental PDB.

    When using an experimental PDB structure (ectodomain-only), the
    AlphaFold model provides full-length coverage for SASA backfill
    of unresolved regions in the visualization.

    Args:
        target: TargetInfo with uniprot_id.
        structures_dir: Directory to save the PDB file.

    Returns:
        Path string to the AlphaFold PDB file, or None if unavailable.
    """
    structures_dir = Path(structures_dir)
    structures_dir.mkdir(parents=True, exist_ok=True)

    dest = structures_dir / "{}_alphafold_supplement.pdb".format(
        target.gene_name.lower()
    )

    # Check if already downloaded
    if dest.exists():
        logger.info("  Using cached AlphaFold supplement: %s", dest.name)
        return str(dest)

    # Query AlphaFold DB API to get latest version
    api_url = "https://alphafold.ebi.ac.uk/api/prediction/{}".format(
        target.uniprot_id
    )
    logger.info("  Downloading AlphaFold supplement for SASA backfill...")
    try:
        resp = requests.get(api_url, timeout=15)
        if resp.status_code != 200:
            logger.info("  No AlphaFold DB entry for %s", target.uniprot_id)
            return None
        entries = resp.json()
        if not entries:
            return None
        latest_version = entries[0].get("latestVersion", 4)
    except (requests.RequestException, ValueError) as e:
        logger.debug("  AlphaFold DB API failed: %s", e)
        return None

    pdb_url = "{}/AF-{}-F1-model_v{}.pdb".format(
        ALPHAFOLD_DB_URL, target.uniprot_id, latest_version
    )
    try:
        resp = requests.get(pdb_url, timeout=30)
        if resp.status_code == 200:
            with open(dest, "w") as f:
                f.write(resp.text)
            logger.info("  Downloaded AlphaFold supplement (v%d): %s",
                        latest_version, dest.name)
            return str(dest)
        else:
            logger.info("  AlphaFold supplement download failed (HTTP %d)",
                        resp.status_code)
    except requests.RequestException as e:
        logger.debug("  AlphaFold supplement download failed: %s", e)

    return None


# ---------------------------------------------------------------------------
# PDB ranking
# ---------------------------------------------------------------------------

def _rank_pdb_structures(pdb_ids, target):
    """
    Rank PDB structures by suitability for surface epitope analysis.

    For each PDB ID:
      1. Query RCSB API for resolution, method, and entity information
      2. Filter: method must be X-ray or cryo-EM, resolution <= 3.5A
      3. Score = coverage_fraction * (1 / resolution)
      4. Sort by score descending

    Args:
        pdb_ids: List of PDB ID strings from UniProt cross-references.
        target: TargetInfo for sequence matching.

    Returns:
        List of dicts sorted by score, each containing:
        pdb_id, method, resolution, chain_id, coverage_start, coverage_end, score
    """
    candidates = []
    # Limit to first 20 PDB entries to avoid excessive API calls
    # (entries are typically ordered by relevance in UniProt)
    check_ids = pdb_ids[:20]
    for pdb_id in check_ids:
        info = _fetch_pdb_info(pdb_id)
        if info is None:
            continue

        method = info.get("method", "")
        resolution = info.get("resolution")

        # Filter by method
        if method not in ALLOWED_METHODS:
            continue

        # Filter by resolution
        if resolution is None or resolution > RESOLUTION_THRESHOLD_A:
            continue

        # Find the chain matching our target protein
        chain_id, cov_start, cov_end = _find_matching_chain(info, target)
        if chain_id is None:
            continue

        # Calculate coverage fraction against ectodomain (not full protein)
        if target.ectodomain_length > 0:
            pdb_range = set(range(cov_start, cov_end + 1))
            ectodomain_set = set()
            for ec_start, ec_end in target.ectodomain_ranges:
                ectodomain_set.update(range(ec_start, ec_end + 1))
            overlap = len(pdb_range & ectodomain_set)
            coverage_fraction = overlap / target.ectodomain_length
        else:
            coverage_fraction = (cov_end - cov_start + 1) / target.sequence_length
        score = coverage_fraction * (1.0 / resolution)

        candidates.append({
            "pdb_id": pdb_id,
            "method": method,
            "resolution": resolution,
            "chain_id": chain_id,
            "coverage_start": cov_start,
            "coverage_end": cov_end,
            "coverage_fraction": coverage_fraction,
            "score": score,
        })

    # Sort by score descending (best first)
    candidates.sort(key=lambda c: c["score"], reverse=True)
    return candidates


def _fetch_pdb_info(pdb_id):
    """
    Fetch structure metadata from RCSB REST API.

    The RCSB API has separate endpoints for entry-level data (resolution,
    method) and entity-level data (chains, UniProt mappings). This function
    queries both and combines them.

    Returns dict with method, resolution, and entity/chain information,
    or None on failure.
    """
    cache_path = CACHE_DIR / "pdb" / "{}.json".format(pdb_id.lower())
    if cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)

    max_attempts = RETRY_RCSB["attempts"]
    backoff = RETRY_RCSB["backoff_seconds"]

    # Fetch entry-level data (resolution, method)
    entry_data = None
    url = "{}/core/entry/{}".format(RCSB_API, pdb_id)
    for attempt in range(1, max_attempts + 1):
        try:
            resp = requests.get(url, timeout=15)
            if resp.status_code == 200:
                entry_data = resp.json()
                break
            elif resp.status_code == 404:
                logger.debug("  PDB %s not found in RCSB", pdb_id)
                return None
            elif resp.status_code in (429, 503):
                time.sleep(backoff * attempt)
            else:
                return None
        except requests.RequestException:
            if attempt < max_attempts:
                time.sleep(backoff * attempt)
            else:
                return None

    if entry_data is None:
        return None

    # Extract method and resolution
    methods = entry_data.get("exptl", [])
    method = methods[0].get("method", "") if methods else ""
    entry_info = entry_data.get("rcsb_entry_info", {})
    resolution = entry_info.get("resolution_combined", [None])[0]

    # Get entity IDs from the entry
    identifiers = entry_data.get("rcsb_entry_container_identifiers", {})
    entity_ids = identifiers.get("polymer_entity_ids", [])

    # Fetch entity-level data for each polymer entity
    entities = []
    for entity_id in entity_ids:
        entity_url = "{}/core/polymer_entity/{}/{}".format(RCSB_API, pdb_id, entity_id)
        try:
            resp = requests.get(entity_url, timeout=15)
            if resp.status_code != 200:
                continue
            edata = resp.json()

            container = edata.get("rcsb_polymer_entity_container_identifiers", {})
            chains = container.get("auth_asym_ids", [])
            uniprot_ids = container.get("uniprot_ids", [])

            ep = edata.get("entity_poly", {})
            sequence = ep.get("pdbx_seq_one_letter_code_can", "")

            entities.append({
                "entity_id": entity_id,
                "chains": chains,
                "uniprot_ids": uniprot_ids if uniprot_ids else [],
                "sequence": sequence,
            })
        except requests.RequestException:
            continue

    info = {
        "pdb_id": pdb_id,
        "method": method,
        "resolution": resolution,
        "entities": entities,
    }

    # Cache
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(info, f, indent=2)

    return info


def _find_matching_chain(pdb_info, target):
    """
    Find the chain in the PDB structure that corresponds to our target protein.

    Matches by UniProt accession in the entity cross-references, or by
    sequence similarity if UniProt mapping is unavailable.

    Uses local sequence alignment to determine the exact coverage range
    of the PDB entity on the target UniProt sequence.

    Returns:
        (chain_id, coverage_start, coverage_end) or (None, None, None).
        coverage_start and coverage_end are 1-based UniProt residue numbers.
    """
    for entity in pdb_info.get("entities", []):
        # Match by UniProt accession
        if target.uniprot_id in entity.get("uniprot_ids", []):
            chain_id = entity["chains"][0] if entity["chains"] else None
            if chain_id and entity.get("sequence"):
                entity_seq = entity["sequence"].replace("(", "").replace(")", "")
                cov_start, cov_end = _align_entity_to_target(
                    entity_seq, target.sequence
                )
                if cov_start is not None:
                    return chain_id, cov_start, cov_end

    # No UniProt match found — try sequence-based matching
    for entity in pdb_info.get("entities", []):
        if not entity.get("sequence"):
            continue
        entity_seq = entity["sequence"].replace("(", "").replace(")", "")
        if len(entity_seq) < 20:
            continue
        chain_id = entity["chains"][0] if entity["chains"] else None
        if chain_id:
            cov_start, cov_end = _align_entity_to_target(
                entity_seq, target.sequence
            )
            if cov_start is not None:
                return chain_id, cov_start, cov_end

    return None, None, None


def _align_entity_to_target(entity_seq, target_seq):
    """
    Align a PDB entity sequence to the target UniProt sequence using
    local alignment to find the coverage range.

    Returns:
        (cov_start, cov_end) as 1-based UniProt positions, or (None, None)
        if alignment quality is too low.
    """
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(target_seq, entity_seq)
    if not alignments:
        return None, None

    best = alignments[0]

    # Extract the aligned intervals in the target (UniProt) sequence
    target_intervals = best.aligned[0]  # List of (start, end) 0-based half-open
    if not len(target_intervals):
        return None, None

    # Coverage range: first to last aligned target position (1-based inclusive)
    cov_start = int(target_intervals[0][0]) + 1
    cov_end = int(target_intervals[-1][1])  # half-open → 1-based inclusive

    # Quality check: aligned length should be >= 50% of entity sequence
    aligned_length = sum(e - s for s, e in target_intervals)
    if aligned_length < len(entity_seq) * 0.5:
        return None, None

    return cov_start, cov_end


# ---------------------------------------------------------------------------
# PDB download
# ---------------------------------------------------------------------------

def _download_pdb(pdb_id, dest_dir):
    """
    Download a PDB file from RCSB.

    Tries .pdb format first, falls back to .cif if unavailable.

    Args:
        pdb_id: 4-character PDB ID.
        dest_dir: Path to save directory.

    Returns:
        Path to the downloaded file.
    """
    dest_dir = Path(dest_dir)
    pdb_path = dest_dir / "{}.pdb".format(pdb_id.lower())

    if pdb_path.exists():
        logger.debug("  PDB %s already downloaded: %s", pdb_id, pdb_path)
        return pdb_path

    # Try PDB format
    url = "{}/{}.pdb".format(RCSB_FILES, pdb_id.upper())
    resp = requests.get(url, timeout=30)
    if resp.status_code == 200:
        with open(pdb_path, "w") as f:
            f.write(resp.text)
        logger.info("  Downloaded PDB: %s", pdb_path)
        return pdb_path

    # Try mmCIF format as fallback
    cif_path = dest_dir / "{}.cif".format(pdb_id.lower())
    url = "{}/{}.cif".format(RCSB_FILES, pdb_id.upper())
    resp = requests.get(url, timeout=30)
    if resp.status_code == 200:
        with open(cif_path, "w") as f:
            f.write(resp.text)
        logger.info("  Downloaded CIF: %s", cif_path)
        return cif_path

    raise StructureAcquisitionError(
        "Failed to download PDB {}: HTTP {}".format(pdb_id, resp.status_code)
    )


# ---------------------------------------------------------------------------
# AlphaFold DB pre-computed structure
# ---------------------------------------------------------------------------

def _try_alphafold_db(target, structures_dir):
    """
    Try to download a pre-computed AlphaFold structure from the EBI database.

    Queries the AlphaFold DB API to discover the latest model version,
    then downloads the PDB file.

    Args:
        target: TargetInfo with uniprot_id and sequence_length.
        structures_dir: Directory to save the PDB file.

    Returns:
        StructureResult if successful, None if not available.
    """
    dest = Path(structures_dir) / "{}_alphafold.pdb".format(target.gene_name.lower())

    # Check if already downloaded
    if dest.exists():
        logger.info("  Using cached AlphaFold DB structure: %s", dest)
        return StructureResult(
            uniprot_id=target.uniprot_id,
            gene_name=target.gene_name,
            pdb_path=str(dest),
            source="alphafold_db",
            pdb_id="AF-{}".format(target.uniprot_id),
            method="AlphaFold",
            resolution=None,
            chain_id="A",
            coverage_start=1,
            coverage_end=target.sequence_length,
        )

    # Query AlphaFold DB API to get latest version
    api_url = "https://alphafold.ebi.ac.uk/api/prediction/{}".format(target.uniprot_id)
    logger.info("  Checking AlphaFold DB for %s...", target.uniprot_id)
    try:
        resp = requests.get(api_url, timeout=15)
        if resp.status_code != 200:
            logger.info("  No AlphaFold DB entry for %s (HTTP %d)", target.uniprot_id, resp.status_code)
            return None
        entries = resp.json()
        if not entries:
            logger.info("  No AlphaFold DB entry for %s", target.uniprot_id)
            return None
        latest_version = entries[0].get("latestVersion", 4)
    except (requests.RequestException, ValueError) as e:
        logger.debug("  AlphaFold DB API failed: %s", e)
        return None

    # Download PDB using discovered version
    pdb_url = "{}/AF-{}-F1-model_v{}.pdb".format(
        ALPHAFOLD_DB_URL, target.uniprot_id, latest_version
    )
    logger.info("  Downloading AlphaFold DB v%d: %s", latest_version, pdb_url)
    try:
        resp = requests.get(pdb_url, timeout=30)
        if resp.status_code == 200:
            with open(dest, "w") as f:
                f.write(resp.text)
            logger.info(
                "  Downloaded AlphaFold DB structure (v%d) for %s: %s",
                latest_version, target.gene_name, dest,
            )
            return StructureResult(
                uniprot_id=target.uniprot_id,
                gene_name=target.gene_name,
                pdb_path=str(dest),
                source="alphafold_db",
                pdb_id="AF-{}".format(target.uniprot_id),
                method="AlphaFold",
                resolution=None,
                chain_id="A",
                coverage_start=1,
                coverage_end=target.sequence_length,
            )
        else:
            logger.info("  AlphaFold DB PDB download failed (HTTP %d)", resp.status_code)
    except requests.RequestException as e:
        logger.debug("  AlphaFold DB download failed: %s", e)

    return None


# ---------------------------------------------------------------------------
# Tamarind AlphaFold fallback
# ---------------------------------------------------------------------------

def _fold_with_tamarind(target, structures_dir):
    """
    Submit a structure prediction job to Tamarind AlphaFold.

    Uses the TamarindClient from the tamarind project. Job results
    (PDB files) are saved to raw_jobs/ and copied to structures_dir.

    Args:
        target: TargetInfo with sequence.
        structures_dir: Directory for the output PDB.

    Returns:
        StructureResult with predicted structure.

    Raises:
        StructureAcquisitionError: If folding fails or times out.
    """
    job_name = "epitope_{}_fold".format(target.gene_name.lower())
    raw_dir = RAW_JOBS_DIR / job_name
    raw_dir.mkdir(parents=True, exist_ok=True)

    # Check if we already have results from a previous run
    existing_pdb = list(raw_dir.glob("*.pdb"))
    if existing_pdb:
        pdb_path = existing_pdb[0]
        dest = Path(structures_dir) / "{}_alphafold.pdb".format(target.gene_name.lower())
        if not dest.exists():
            import shutil
            shutil.copy2(str(pdb_path), str(dest))
        logger.info("  Using cached AlphaFold prediction: %s", dest)
        return StructureResult(
            uniprot_id=target.uniprot_id,
            gene_name=target.gene_name,
            pdb_path=str(dest),
            source="tamarind_alphafold",
            pdb_id=job_name,
            method="AlphaFold",
            resolution=None,
            chain_id="A",  # AlphaFold predictions use chain A
            coverage_start=1,
            coverage_end=target.sequence_length,
        )

    # Submit to Tamarind
    try:
        client = TamarindClient()

        settings = {
            "sequence": target.sequence,
            "numModels": str(FOLDING_NUM_MODELS),  # API expects string options
            "numRecycles": FOLDING_NUM_RECYCLES,
        }

        try:
            client.submit_job(job_name, FOLDING_TOOL, settings)
            logger.info("  Submitted AlphaFold job: %s", job_name)
        except TamarindError as e:
            if "already exists" in str(e).lower():
                logger.info("  Job %s already exists, waiting for completion...", job_name)
            else:
                raise

        # Wait for completion
        client.wait_for_job(job_name, timeout=FOLDING_TIMEOUT, verbose=True)
        logger.info("  AlphaFold job completed: %s", job_name)

        # Download results
        result = client.get_result(job_name)
        pdb_path = _download_fold_result(result, raw_dir, target, structures_dir)

        return StructureResult(
            uniprot_id=target.uniprot_id,
            gene_name=target.gene_name,
            pdb_path=str(pdb_path),
            source="tamarind_alphafold",
            pdb_id=job_name,
            method="AlphaFold",
            resolution=None,
            chain_id="A",
            coverage_start=1,
            coverage_end=target.sequence_length,
        )

    except (TamarindError, Exception) as e:
        raise StructureAcquisitionError(
            "Tamarind AlphaFold folding failed for {}: {}".format(target.gene_name, e)
        )


def _download_fold_result(result, raw_dir, target, structures_dir):
    """
    Download and extract AlphaFold prediction result.

    The Tamarind API returns either a download URL (for zip archives)
    or a direct data dict. Handles both cases.

    Returns:
        Path to the extracted PDB file in structures_dir.
    """
    if isinstance(result, str):
        # Result is a download URL
        resp = requests.get(result, timeout=60)
        resp.raise_for_status()

        # Try as zip archive
        try:
            with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
                zf.extractall(str(raw_dir))
                # Find the PDB file in extracted contents
                pdb_files = [f for f in zf.namelist() if f.endswith(".pdb")]
                if pdb_files:
                    src = raw_dir / pdb_files[0]
                    dest = Path(structures_dir) / "{}_alphafold.pdb".format(
                        target.gene_name.lower()
                    )
                    import shutil
                    shutil.copy2(str(src), str(dest))
                    return dest
        except zipfile.BadZipFile:
            # Not a zip — probably raw PDB content
            dest = Path(structures_dir) / "{}_alphafold.pdb".format(
                target.gene_name.lower()
            )
            with open(dest, "wb") as f:
                f.write(resp.content)
            return dest

    elif isinstance(result, dict):
        # Direct data dict — look for PDB content or URL
        if "url" in result:
            return _download_fold_result(result["url"], raw_dir, target, structures_dir)
        elif "pdb" in result:
            dest = Path(structures_dir) / "{}_alphafold.pdb".format(
                target.gene_name.lower()
            )
            with open(dest, "w") as f:
                f.write(result["pdb"])
            return dest

    raise StructureAcquisitionError(
        "Could not extract PDB from AlphaFold result for {}".format(target.gene_name)
    )
