"""
Step 6: Cynomolgus Monkey Conservation — map sequence conservation onto
surface patches and filter by percentage of mismatched residues.

Performs pairwise Needleman-Wunsch alignment between human and cyno
ortholog sequences, maps conservation per-residue onto the 3D structure,
and evaluates each surface patch as an atomic unit. Patches with ≤15%
mismatched residues pass filtering (default threshold).
"""

import logging
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from epitope_pipeline import config
from epitope_pipeline.config import (
    MAX_CYNO_MISMATCHES_PER_600A2,  # Deprecated, kept for backward compatibility
    PATCH_CLUSTERING_DISTANCE_A,
    VHH_FOOTPRINT_MIN_A2,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ConservationResult:
    """Per-residue and per-patch conservation analysis."""
    uniprot_id: str
    gene_name: str
    alignment_human: str                    # Aligned human sequence (with gaps)
    alignment_cyno: str                     # Aligned cyno sequence (with gaps)
    overall_identity: float                 # Overall sequence identity fraction
    residue_conservation: dict              # {human_resnum: True/False (identical)}
    conserved_patches: list                 # SurfacePatch objects passing threshold
    rejected_patches: list                  # [(patch, identity_fraction), ...]
    patch_conservation: dict                # {patch_id: identity_fraction}


class ConservationError(Exception):
    """Conservation analysis failed."""
    pass


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def analyze_conservation(target, surface_analysis, ca_coords=None):
    """
    Map cynomolgus monkey conservation onto surface patches and filter.

    Steps:
      1. Align human and cyno sequences (Needleman-Wunsch, BLOSUM62)
      2. Map alignment to per-residue identity
      3. For each surface patch: calculate percentage of mismatched residues
      4. Accept patches with ≤15% mismatches (configurable threshold)

    Args:
        target: TargetInfo with sequence and cyno_sequence.
        surface_analysis: SurfaceAnalysis with patches to evaluate.
        ca_coords: Optional dict {resnum: np.array} (deprecated, not used).

    Returns:
        ConservationResult with filtered patches and conservation data.

    Raises:
        ConservationError: If no cyno ortholog sequence is available.
    """
    if not target.cyno_sequence:
        raise ConservationError(
            "No cynomolgus monkey ortholog sequence available for {}. "
            "Cannot assess conservation.".format(target.gene_name)
        )

    # Step 1: Pairwise alignment (always run — per-residue data needed for figure)
    aligned_human, aligned_cyno = _align_sequences(target.sequence, target.cyno_sequence)

    # Step 2: Map alignment positions to human residue conservation
    residue_conservation = _map_alignment_to_residues(aligned_human, aligned_cyno)

    # Overall identity
    total_aligned = sum(1 for h, c in zip(aligned_human, aligned_cyno) if h != "-" and c != "-")
    total_identical = sum(1 for h, c in zip(aligned_human, aligned_cyno) if h == c and h != "-")
    overall_identity = total_identical / total_aligned if total_aligned > 0 else 0.0

    logger.info(
        "  %s: Overall cyno identity: %.1f%% (%d/%d aligned positions)",
        target.gene_name, overall_identity * 100, total_identical, total_aligned,
    )

    if not surface_analysis.patches:
        logger.info("  %s: No patches to assess for conservation", target.gene_name)
        return ConservationResult(
            uniprot_id=target.uniprot_id,
            gene_name=target.gene_name,
            alignment_human=aligned_human,
            alignment_cyno=aligned_cyno,
            overall_identity=overall_identity,
            residue_conservation=residue_conservation,
            conserved_patches=[],
            rejected_patches=[],
            patch_conservation={},
        )

    # Step 3: Evaluate each patch (whole-patch mode, no trimming)
    conserved_patches = []
    rejected_patches = []
    patch_conservation = {}

    for patch in surface_analysis.patches:
        # Calculate patch-level metrics
        n_residues = len(patch.residue_numbers)
        n_conserved = sum(
            1 for r in patch.residue_numbers
            if residue_conservation.get(r, False)
        )
        n_mismatches = n_residues - n_conserved
        identity_fraction = n_conserved / n_residues if n_residues > 0 else 0.0
        mismatch_percent = (n_mismatches / n_residues) * 100.0

        patch_conservation[patch.patch_id] = identity_fraction

        logger.info(
            "    Patch %d: %d residues (%.0f A²), cyno identity %.1f%% "
            "(%d conserved, %d mismatches = %.1f%%)",
            patch.patch_id, n_residues, patch.total_sasa_a2,
            identity_fraction * 100, n_conserved, n_mismatches, mismatch_percent,
        )

        # Whole-patch evaluation
        passes, _ = evaluate_patch_conservation(patch, residue_conservation)

        if passes:
            conserved_patches.append(patch)
            logger.info("      -> PASSED (%.1f%% mismatches < %.1f%% threshold)",
                       mismatch_percent, config.MAX_CYNO_MISMATCH_PERCENT)
        else:
            rejected_patches.append((patch, identity_fraction))
            logger.info("      -> REJECTED (%.1f%% mismatches > %.1f%% threshold)",
                       mismatch_percent, config.MAX_CYNO_MISMATCH_PERCENT)

    logger.info(
        "  %s conservation: %d patches pass, %d rejected",
        target.gene_name, len(conserved_patches), len(rejected_patches),
    )

    return ConservationResult(
        uniprot_id=target.uniprot_id,
        gene_name=target.gene_name,
        alignment_human=aligned_human,
        alignment_cyno=aligned_cyno,
        overall_identity=overall_identity,
        residue_conservation=residue_conservation,
        conserved_patches=conserved_patches,
        rejected_patches=rejected_patches,
        patch_conservation=patch_conservation,
    )


# ---------------------------------------------------------------------------
# Sequence alignment
# ---------------------------------------------------------------------------

def _align_sequences(human_seq, cyno_seq):
    """
    Pairwise global alignment using BioPython PairwiseAligner with BLOSUM62.

    Args:
        human_seq: Human protein sequence string.
        cyno_seq: Cyno ortholog sequence string.

    Returns:
        Tuple (aligned_human, aligned_cyno) with gap characters.
    """
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(human_seq, cyno_seq)
    best = alignments[0]

    # Extract aligned sequences as strings
    # BioPython PairwiseAligner alignment formatting
    formatted = str(best).split("\n")
    # The formatted output has lines for target, match, query
    # But it's safer to reconstruct from the alignment object
    aligned_human = ""
    aligned_cyno = ""

    # Use the alignment's aligned property to reconstruct gapped sequences
    target_aligned = best.aligned[0]  # intervals in target (human)
    query_aligned = best.aligned[1]   # intervals in query (cyno)

    # Reconstruct aligned sequences from interval pairs
    human_pos = 0
    cyno_pos = 0
    ah = []
    ac = []

    for (t_start, t_end), (q_start, q_end) in zip(target_aligned, query_aligned):
        # Add gaps for any positions skipped in human
        while human_pos < t_start:
            ah.append(human_seq[human_pos])
            ac.append("-")
            human_pos += 1
        # Add gaps for any positions skipped in cyno
        while cyno_pos < q_start:
            ah.append("-")
            ac.append(cyno_seq[cyno_pos])
            cyno_pos += 1
        # Add aligned region
        t_len = t_end - t_start
        q_len = q_end - q_start
        # Both should be the same length in a global alignment block
        for i in range(max(t_len, q_len)):
            if i < t_len:
                ah.append(human_seq[t_start + i])
            else:
                ah.append("-")
            if i < q_len:
                ac.append(cyno_seq[q_start + i])
            else:
                ac.append("-")
        human_pos = t_end
        cyno_pos = q_end

    # Add any remaining trailing residues
    while human_pos < len(human_seq):
        ah.append(human_seq[human_pos])
        ac.append("-")
        human_pos += 1
    while cyno_pos < len(cyno_seq):
        ah.append("-")
        ac.append(cyno_seq[cyno_pos])
        cyno_pos += 1

    aligned_human = "".join(ah)
    aligned_cyno = "".join(ac)

    return aligned_human, aligned_cyno


def _map_alignment_to_residues(aligned_human, aligned_cyno):
    """
    Walk through alignment columns and map to human residue numbers.

    For each non-gap position in the human sequence, determines whether
    the corresponding cyno position is identical.

    Args:
        aligned_human: Aligned human sequence string (with gaps).
        aligned_cyno: Aligned cyno sequence string (with gaps).

    Returns:
        Dict {residue_number: is_identical (bool)}.
        Residue numbers are 1-based, matching UniProt convention.
    """
    conservation = {}
    human_resnum = 0

    for h_aa, c_aa in zip(aligned_human, aligned_cyno):
        if h_aa == "-":
            # Gap in human — no residue to map
            continue
        human_resnum += 1

        if c_aa == "-":
            # Gap in cyno — not conserved
            conservation[human_resnum] = False
        else:
            conservation[human_resnum] = (h_aa == c_aa)

    return conservation


# ---------------------------------------------------------------------------
# Whole-patch evaluation (NEW)
# ---------------------------------------------------------------------------

def evaluate_patch_conservation(patch, residue_conservation):
    """
    Evaluate a patch based on percentage of mismatched residues.

    Args:
        patch: SurfacePatch object.
        residue_conservation: Dict {resnum: bool} — True=conserved, False=mismatch.

    Returns:
        Tuple (passes, mismatch_percent):
            - passes: True if patch meets criteria.
            - mismatch_percent: Float percentage (0-100).
    """
    from .config import MAX_CYNO_MISMATCH_PERCENT

    n_residues = len(patch.residue_numbers)
    if n_residues == 0:
        return False, 100.0

    n_conserved = sum(
        1 for r in patch.residue_numbers
        if residue_conservation.get(r, False)
    )
    n_mismatches = n_residues - n_conserved
    mismatch_percent = (n_mismatches / n_residues) * 100.0

    passes = (mismatch_percent <= MAX_CYNO_MISMATCH_PERCENT)

    return passes, mismatch_percent


# ---------------------------------------------------------------------------
# DEPRECATED: Sliding window mismatch density check (kept for reference)
# ---------------------------------------------------------------------------

def check_local_mismatch_density(patch_residues, residue_conservation, ca_coords):
    """
    Sliding window: for every residue in the patch, examine its ~600A²
    neighborhood (radius = sqrt(600/π) ≈ 13.8A) and count mismatches.

    Returns (worst_resnum, worst_mismatch_count) — the residue whose
    neighborhood has the most mismatches.

    Args:
        patch_residues: List of residue numbers in the patch.
        residue_conservation: Dict {resnum: bool} — True=OK, False=mismatch.
        ca_coords: Dict {resnum: np.array} Cα coordinates.

    Returns:
        Tuple (worst_resnum, worst_count). (0, 0) if no coords available.
    """
    from scipy.spatial.distance import cdist

    footprint_radius = math.sqrt(VHH_FOOTPRINT_MIN_A2 / math.pi)

    res_with_coords = [r for r in patch_residues if r in ca_coords]
    if not res_with_coords:
        return (0, 0)

    coords = np.array([ca_coords[r] for r in res_with_coords])
    is_bad = np.array([
        not residue_conservation.get(r, False) for r in res_with_coords
    ])

    dists = cdist(coords, coords)

    worst_pos = 0
    worst_count = 0
    for i, resnum in enumerate(res_with_coords):
        neighbors = dists[i] <= footprint_radius
        local_bad = int(np.sum(is_bad[neighbors]))
        if local_bad > worst_count:
            worst_count = local_bad
            worst_pos = resnum

    return (worst_pos, worst_count)


# Backward-compatible alias for the old private name
_check_local_mismatch_density = check_local_mismatch_density


def identify_failing_residues(patch_residues, boolean_map, ca_coords,
                              max_bad_per_window):
    """
    Identify all residues whose sliding window exceeds the mismatch threshold.

    A residue "fails" if its ~600A² neighborhood contains more than
    max_bad_per_window mismatches. These are the residues to trim.

    Args:
        patch_residues: List of residue numbers in the patch.
        boolean_map: Dict {resnum: bool} — True=OK, False=bad.
        ca_coords: Dict {resnum: np.array} Cα coordinates.
        max_bad_per_window: Maximum allowed bad residues per window.

    Returns:
        Set of residue numbers that should be trimmed.
    """
    from scipy.spatial.distance import cdist

    footprint_radius = math.sqrt(VHH_FOOTPRINT_MIN_A2 / math.pi)

    res_with_coords = [r for r in patch_residues if r in ca_coords]
    if not res_with_coords:
        return set()

    coords = np.array([ca_coords[r] for r in res_with_coords])
    is_bad = np.array([
        not boolean_map.get(r, False) for r in res_with_coords
    ])

    dists = cdist(coords, coords)

    # A residue fails if its window exceeds the threshold
    failing = set()
    for i, resnum in enumerate(res_with_coords):
        neighbors = dists[i] <= footprint_radius
        local_bad = int(np.sum(is_bad[neighbors]))
        if local_bad > max_bad_per_window:
            failing.add(resnum)

    return failing


# ---------------------------------------------------------------------------
# Trim-and-recluster
# ---------------------------------------------------------------------------

def _recluster_survivors(survivors, ca_coords, residue_sasa, original_patch,
                         residue_conservation, next_patch_id):
    """
    Re-cluster survivor residues after trimming and build new SurfacePatch
    objects for sub-clusters >= 600A².

    Also applies the conserved sub-patch check (Gate 2) to each sub-patch.

    Args:
        survivors: List of residue numbers that survived trimming.
        ca_coords: Dict {resnum: np.array} Cα coordinates.
        residue_sasa: Dict {resnum: sasa_A2}.
        original_patch: The original SurfacePatch (for inheriting properties).
        residue_conservation: Dict {resnum: bool}.
        next_patch_id: Starting patch ID for new sub-patches.

    Returns:
        List of SurfacePatch objects that pass area and sub-patch checks.
    """
    from epitope_pipeline.surface import cluster_surface_patches, SurfacePatch

    survivors_with_ca = [r for r in survivors if r in ca_coords]
    if len(survivors_with_ca) < 2:
        return []

    clusters = cluster_surface_patches(survivors_with_ca, ca_coords)

    result = []
    pid = next_patch_id
    for cluster_residues in clusters:
        total_area = sum(residue_sasa.get(r, 0.0) for r in cluster_residues)
        if total_area < VHH_FOOTPRINT_MIN_A2:
            continue

        # Gate 2: conserved sub-patch check on this sub-cluster
        conserved_in_cluster = [
            r for r in cluster_residues
            if residue_conservation.get(r, False) and r in ca_coords
        ]
        if not _check_conserved_subpatch(conserved_in_cluster, ca_coords, residue_sasa):
            logger.info(
                "      Sub-patch (%d residues, %.0f A²): conserved residues "
                "don't form contiguous >= %.0f A², skipping",
                len(cluster_residues), total_area, VHH_FOOTPRINT_MIN_A2,
            )
            continue

        # Build the sub-patch
        coords = np.array([ca_coords[r] for r in cluster_residues])
        centroid = np.mean(coords, axis=0)

        if len(coords) > 1:
            from scipy.spatial.distance import pdist
            max_dim = float(np.max(pdist(coords)))
        else:
            max_dim = 0.0

        sp = SurfacePatch(
            patch_id=pid,
            residue_numbers=sorted(cluster_residues),
            residue_aas={r: original_patch.residue_aas.get(r, "X")
                         for r in cluster_residues},
            total_sasa_a2=total_area,
            centroid=centroid,
            max_dimension_a=max_dim,
            avg_distance_from_membrane=original_patch.avg_distance_from_membrane,
        )
        result.append(sp)
        pid += 1

    return result


# ---------------------------------------------------------------------------
# 3D conserved sub-patch check
# ---------------------------------------------------------------------------

def _check_conserved_subpatch(conserved_residues, ca_coords, residue_sasa):
    """
    Verify that conserved residues within a patch form a contiguous 3D
    region large enough for a VHH CDR footprint.

    Re-clusters only the conserved residues and checks that the largest
    sub-cluster has sufficient surface area.

    Args:
        conserved_residues: List of conserved residue numbers in the patch.
        ca_coords: Dict {resnum: np.array} Cα coordinates.
        residue_sasa: Dict {resnum: sasa_A2} for area calculation.

    Returns:
        True if the largest conserved sub-patch >= VHH footprint minimum.
    """
    if len(conserved_residues) < 2:
        # A single conserved residue can't form a meaningful patch
        area = sum(residue_sasa.get(r, 0.0) for r in conserved_residues)
        return area >= VHH_FOOTPRINT_MIN_A2

    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import fcluster, linkage

    # Build coordinate matrix
    res_list = sorted(conserved_residues)
    coords = np.array([ca_coords[r] for r in res_list])

    # Cluster at the same distance threshold
    dist_matrix = pdist(coords)
    Z = linkage(dist_matrix, method="single")
    labels = fcluster(Z, t=PATCH_CLUSTERING_DISTANCE_A, criterion="distance")

    # Find the largest sub-cluster by area
    clusters = {}
    for res, label in zip(res_list, labels):
        clusters.setdefault(label, []).append(res)

    max_area = 0.0
    for cluster_residues in clusters.values():
        area = sum(residue_sasa.get(r, 0.0) for r in cluster_residues)
        if area > max_area:
            max_area = area

    return max_area >= VHH_FOOTPRINT_MIN_A2
