"""
Step 8: Epitope Scoring — rank surviving epitope patches by a composite
score reflecting their suitability for VHH antibody targeting.

Score components:
  - Area (25%): larger patches provide more binding options
  - Distance (20%): farther from membrane = better accessibility
  - Conservation (25%): higher cyno identity = better safety profile
  - Specificity (20%): more unique to target = lower cross-reactivity
  - Accessibility (10%): higher relative SASA = better solvent exposure
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional

from epitope_pipeline.config import (
    SCORE_WEIGHTS,
    VHH_FOOTPRINT_MAX_A2,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class EpitopeScore:
    """Scored and ranked epitope patch."""
    patch_id: int
    uniprot_id: str
    gene_name: str
    patch: object                       # SurfacePatch reference

    # Component scores (0-1, higher is better)
    area_score: float
    distance_score: float
    conservation_score: float
    specificity_score: float
    accessibility_score: float

    composite_score: float              # Weighted sum
    rank: int                           # 1-based rank within this target

    # Raw values for reporting
    patch_area_a2: float
    avg_distance_a: float
    cyno_identity: float
    max_off_target_identity: float
    avg_relative_sasa: float
    n_residues: int
    residue_range: str                  # e.g. "45-95, 120-145"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def score_epitopes(
    target,
    specific_patches,
    conservation_result,
    specificity_result,
    spatial_filter,
    surface_analysis,
    distance_mode="distal",
):
    """
    Score and rank all qualifying epitope patches for a target.

    Args:
        target: TargetInfo object.
        specific_patches: List of SurfacePatch objects passing all filters.
        conservation_result: ConservationResult with patch conservation data.
        specificity_result: SpecificityResult with BLAST hit data.
        spatial_filter: SpatialFilter with distance data.
        surface_analysis: SurfaceAnalysis with SASA data.
        distance_mode: "distal" (reward high distance) or "proximal" (flat 1.0).

    Returns:
        List of EpitopeScore objects sorted by composite_score (best first).
        Returns empty list if no patches qualify.
    """
    if not specific_patches:
        logger.info("  %s: No patches to score (score = 0)", target.gene_name)
        return []

    scores = []
    for patch in specific_patches:
        score = _score_patch(
            patch, target, conservation_result, specificity_result,
            spatial_filter, surface_analysis, distance_mode=distance_mode,
        )
        scores.append(score)

    # Sort by composite score descending
    scores.sort(key=lambda s: s.composite_score, reverse=True)

    # Assign ranks
    for rank, score in enumerate(scores, start=1):
        score.rank = rank

    # Log summary
    logger.info("  %s epitope scores:", target.gene_name)
    for s in scores:
        logger.info(
            "    Rank %d (patch %d): composite=%.3f, area=%.0fA², "
            "dist=%.0fA, cyno=%.1f%%, residues=%s",
            s.rank, s.patch_id, s.composite_score, s.patch_area_a2,
            s.avg_distance_a, s.cyno_identity * 100, s.residue_range,
        )

    return scores


def compute_target_epitope_metric(scores, surface_analysis):
    """
    Compute an overall epitope space metric for a target.

    The metric represents the fraction of total ectodomain surface that
    qualifies as druggable epitope space.

    Args:
        scores: List of EpitopeScore objects.
        surface_analysis: SurfaceAnalysis for the target.

    Returns:
        Dict with:
          total_epitope_area_a2: sum of qualifying patch areas
          total_ectodomain_area_a2: total ectodomain surface
          epitope_fraction: ratio
          n_patches: number of qualifying patches
          best_score: highest composite score
    """
    total_epitope = sum(s.patch_area_a2 for s in scores)

    # Total ectodomain surface = sum of SASA for all exposed residues
    total_ectodomain = sum(
        sasa for resnum, sasa in surface_analysis.residue_sasa.items()
        if resnum in surface_analysis.exposed_residues
    ) if surface_analysis.exposed_residues else 0.0

    fraction = total_epitope / total_ectodomain if total_ectodomain > 0 else 0.0

    return {
        "total_epitope_area_a2": total_epitope,
        "total_ectodomain_area_a2": total_ectodomain,
        "epitope_fraction": fraction,
        "n_patches": len(scores),
        "best_score": scores[0].composite_score if scores else 0.0,
    }


# ---------------------------------------------------------------------------
# Internal scoring
# ---------------------------------------------------------------------------

def _score_patch(patch, target, conservation_result, specificity_result,
                 spatial_filter, surface_analysis, distance_mode="distal"):
    """
    Compute component scores and composite score for a single patch.
    """
    # --- Area score ---
    # Normalize by the typical VHH footprint max (900 A²)
    area = patch.total_sasa_a2
    area_score = min(area / VHH_FOOTPRINT_MAX_A2, 1.0)

    # --- Distance score ---
    avg_dist = patch.avg_distance_from_membrane
    if distance_mode == "proximal":
        # Proximal mode: flat 1.0 — filter already ensures <= threshold
        distance_score = 1.0
    else:
        # Distal mode: reward HIGH distance from membrane
        distance_score = min(avg_dist / 300.0, 1.0)

    # --- Conservation score ---
    cyno_identity = conservation_result.patch_conservation.get(patch.patch_id, 0.0)
    conservation_score = cyno_identity

    # --- Specificity score ---
    # Check BLAST hits for this patch: 1 - max_off_target_identity
    max_off_target = 0.0
    patch_hits = specificity_result.blast_results.get(patch.patch_id, [])
    for hit in patch_hits:
        if hit["identity"] > max_off_target:
            max_off_target = hit["identity"]
    specificity_score = 1.0 - max_off_target

    # --- Accessibility score ---
    # Average relative SASA of patch residues
    rel_sasa_values = [
        surface_analysis.residue_relative_sasa.get(r, 0.0)
        for r in patch.residue_numbers
    ]
    avg_rel_sasa = sum(rel_sasa_values) / len(rel_sasa_values) if rel_sasa_values else 0.0
    accessibility_score = min(avg_rel_sasa / 0.6, 1.0)  # Normalize; >60% rel SASA is very exposed

    # --- Composite ---
    w = SCORE_WEIGHTS
    composite = (
        w["area"] * area_score +
        w["distance"] * distance_score +
        w["conservation"] * conservation_score +
        w["specificity"] * specificity_score +
        w["accessibility"] * accessibility_score
    )

    # Build residue range string
    residue_range = _format_residue_ranges(patch.residue_numbers)

    return EpitopeScore(
        patch_id=patch.patch_id,
        uniprot_id=target.uniprot_id,
        gene_name=target.gene_name,
        patch=patch,
        area_score=area_score,
        distance_score=distance_score,
        conservation_score=conservation_score,
        specificity_score=specificity_score,
        accessibility_score=accessibility_score,
        composite_score=composite,
        rank=0,  # Assigned after sorting
        patch_area_a2=area,
        avg_distance_a=avg_dist,
        cyno_identity=cyno_identity,
        max_off_target_identity=max_off_target,
        avg_relative_sasa=avg_rel_sasa,
        n_residues=len(patch.residue_numbers),
        residue_range=residue_range,
    )


def _format_residue_ranges(residue_numbers):
    """
    Format a list of residue numbers into a compact range string.

    Example: [45, 46, 47, 50, 51, 78, 79, 80] -> "45-47, 50-51, 78-80"
    """
    if not residue_numbers:
        return ""

    sorted_res = sorted(residue_numbers)
    ranges = []
    start = sorted_res[0]
    end = sorted_res[0]

    for r in sorted_res[1:]:
        if r == end + 1:
            end = r
        else:
            if start == end:
                ranges.append(str(start))
            else:
                ranges.append("{}-{}".format(start, end))
            start = r
            end = r

    # Last range
    if start == end:
        ranges.append(str(start))
    else:
        ranges.append("{}-{}".format(start, end))

    return ", ".join(ranges)
