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

        # Whole-patch evaluation (threshold scales with patch size)
        passes, _, effective_thresh = evaluate_patch_conservation(patch, residue_conservation)

        if passes:
            conserved_patches.append(patch)
            logger.info("      -> PASSED (%.1f%% mismatches <= %.1f%% threshold, %d residues)",
                       mismatch_percent, effective_thresh, n_residues)
        else:
            rejected_patches.append((patch, identity_fraction))
            logger.info("      -> REJECTED (%.1f%% mismatches > %.1f%% threshold, %d residues)",
                       mismatch_percent, effective_thresh, n_residues)

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

    The threshold scales with patch size: larger patches tolerate a higher
    mismatch percentage because the VHH only contacts ~20 residues within
    the patch. Formula:

        effective_threshold = min(base_threshold * sqrt(n_residues / 20), 30%)

    A 20-residue patch uses the base threshold (default 15%). A 100-residue
    patch gets ~34% → capped at 30%. This reflects the reality that
    mismatches in a large patch can be avoided spatially.

    Args:
        patch: SurfacePatch object.
        residue_conservation: Dict {resnum: bool} — True=conserved, False=mismatch.

    Returns:
        Tuple (passes, mismatch_percent, effective_threshold):
            - passes: True if patch meets criteria.
            - mismatch_percent: Float percentage (0-100).
            - effective_threshold: Size-adjusted threshold used (0-100).
    """
    from epitope_pipeline.config import MAX_CYNO_MISMATCH_PERCENT

    n_residues = len(patch.residue_numbers)
    if n_residues == 0:
        return False, 100.0, 0.0

    n_conserved = sum(
        1 for r in patch.residue_numbers
        if residue_conservation.get(r, False)
    )
    n_mismatches = n_residues - n_conserved
    mismatch_percent = (n_mismatches / n_residues) * 100.0

    # Scale threshold by patch size: sqrt(n/20) with 30% ceiling
    effective_threshold = min(
        MAX_CYNO_MISMATCH_PERCENT * math.sqrt(n_residues / 20.0),
        30.0,
    )
    passes = (mismatch_percent <= effective_threshold)

    return passes, mismatch_percent, effective_threshold
