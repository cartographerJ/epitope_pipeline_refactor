"""Tests for target-level epitope metrics and residue-range formatting."""

import math

from epitope_pipeline.compute.scoring import (
    compute_target_epitope_metric,
    _format_residue_ranges,
)
from tests._fixtures import make_score, make_surface_analysis


# ---------------------------------------------------------------------------
# compute_target_epitope_metric
# ---------------------------------------------------------------------------

def test_empty_scores_returns_zeros():
    surface = make_surface_analysis()
    metric = compute_target_epitope_metric([], surface)
    assert metric["n_patches"] == 0
    assert metric["total_epitope_area_a2"] == 0.0
    assert metric["best_score"] == 0.0


def test_single_score_basic():
    s = make_score(patch_area_a2=1000.0, cyno_identity=0.9,
                   n_residues=25, composite_score=0.7)
    surface = make_surface_analysis(
        residue_sasa={r: 100.0 for r in range(1, 11)},
        exposed_residues=list(range(1, 11)),
    )
    metric = compute_target_epitope_metric([s], surface)
    assert metric["n_patches"] == 1
    assert metric["total_epitope_area_a2"] == 1000.0
    assert metric["total_ectodomain_area_a2"] == 1000.0  # 10 * 100
    assert math.isclose(metric["epitope_fraction"], 1.0, rel_tol=1e-6)
    assert metric["best_score"] == 0.7


def test_quality_component_area_weighted():
    # Two patches of identical area → quality is plain mean of cyno_identity
    s1 = make_score(patch_area_a2=500.0, cyno_identity=0.8, n_residues=20)
    s2 = make_score(patch_area_a2=500.0, cyno_identity=1.0, n_residues=20)
    surface = make_surface_analysis()
    metric = compute_target_epitope_metric([s1, s2], surface)
    assert math.isclose(metric["quality_component"], 0.9, rel_tol=1e-6)


def test_quality_component_weighting_by_residues():
    # 30-residue patch with 100% identity dominates 10-residue patch with 0%
    s1 = make_score(patch_area_a2=400.0, cyno_identity=0.0, n_residues=10)
    s2 = make_score(patch_area_a2=1200.0, cyno_identity=1.0, n_residues=30)
    surface = make_surface_analysis()
    metric = compute_target_epitope_metric([s1, s2], surface)
    expected = (0.0 * 10 + 1.0 * 30) / 40.0
    assert math.isclose(metric["quality_component"], expected, rel_tol=1e-6)


def test_area_component_capped_at_5000():
    # area_component = min(total_epitope / 5000.0, 1.0)
    big = make_score(patch_area_a2=10000.0, cyno_identity=1.0, n_residues=50)
    surface = make_surface_analysis()
    metric = compute_target_epitope_metric([big], surface)
    assert metric["area_component"] == 1.0


def test_area_component_proportional_below_5000():
    s = make_score(patch_area_a2=2500.0, cyno_identity=1.0, n_residues=50)
    surface = make_surface_analysis()
    metric = compute_target_epitope_metric([s], surface)
    assert math.isclose(metric["area_component"], 0.5, rel_tol=1e-6)


def test_zero_ectodomain_safe():
    """No exposed residues → fraction = 0, no division-by-zero."""
    s = make_score(patch_area_a2=1000.0)
    surface = make_surface_analysis(residue_sasa={}, exposed_residues=[])
    metric = compute_target_epitope_metric([s], surface)
    assert metric["epitope_fraction"] == 0.0
    assert metric["total_ectodomain_area_a2"] == 0.0


# ---------------------------------------------------------------------------
# _format_residue_ranges
# ---------------------------------------------------------------------------

def test_empty_range():
    assert _format_residue_ranges([]) == ""


def test_single_residue():
    assert _format_residue_ranges([42]) == "42"


def test_contiguous_range():
    assert _format_residue_ranges([3, 4, 5, 6, 7]) == "3-7"


def test_disjoint_ranges():
    assert _format_residue_ranges([45, 46, 47, 50, 51, 78, 79, 80]) == \
        "45-47, 50-51, 78-80"


def test_unsorted_input():
    assert _format_residue_ranges([5, 1, 4, 2, 3]) == "1-5"


def test_singletons_mixed():
    assert _format_residue_ranges([1, 5, 10]) == "1, 5, 10"
