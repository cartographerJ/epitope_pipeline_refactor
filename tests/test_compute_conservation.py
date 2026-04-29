"""Tests for the size-scaled patch conservation rule.

Threshold formula: min(MAX_CYNO_MISMATCH_PERCENT * sqrt(n_residues / 20), 30%).
"""

import math

from epitope_pipeline import config
from epitope_pipeline.compute.conservation import evaluate_patch_conservation
from tests._fixtures import make_patch


def test_all_conserved_passes():
    patch = make_patch(residue_numbers=list(range(1, 21)))
    cons = {r: True for r in patch.residue_numbers}
    passes, mismatch_pct, threshold = evaluate_patch_conservation(patch, cons)
    assert passes is True
    assert mismatch_pct == 0.0
    # 20-residue patch → sqrt(20/20)=1, threshold == base
    assert math.isclose(threshold, config.MAX_CYNO_MISMATCH_PERCENT, rel_tol=1e-6)


def test_at_threshold_passes():
    # 20 residues, base 15% → threshold 15%; 3 mismatches = exactly 15%
    base = config.MAX_CYNO_MISMATCH_PERCENT
    n = 20
    n_mismatches = int(round(n * base / 100.0))
    cons = {r: (i >= n_mismatches) for i, r in enumerate(range(1, n + 1))}
    patch = make_patch(residue_numbers=list(range(1, n + 1)))
    passes, pct, threshold = evaluate_patch_conservation(patch, cons)
    assert passes is True
    assert math.isclose(pct, base, abs_tol=0.01)


def test_above_threshold_fails():
    # 20 residues, 4 mismatches = 20% > 15% threshold
    n = 20
    cons = {r: (i >= 4) for i, r in enumerate(range(1, n + 1))}
    patch = make_patch(residue_numbers=list(range(1, n + 1)))
    passes, pct, threshold = evaluate_patch_conservation(patch, cons)
    assert passes is False
    assert math.isclose(pct, 20.0, abs_tol=0.01)


def test_threshold_scales_with_size():
    # 80-residue patch should tolerate higher mismatch %: sqrt(80/20)=2 → 2x base
    patch = make_patch(residue_numbers=list(range(1, 81)))
    cons = {r: True for r in patch.residue_numbers}
    _, _, threshold = evaluate_patch_conservation(patch, cons)
    expected = config.MAX_CYNO_MISMATCH_PERCENT * 2.0
    assert math.isclose(threshold, expected, rel_tol=1e-6)


def test_threshold_capped_at_30_percent():
    # Very large patch: base * sqrt(500/20) ≈ base*5 ≫ 30 → must be capped at 30
    patch = make_patch(residue_numbers=list(range(1, 501)))
    cons = {r: True for r in patch.residue_numbers}
    _, _, threshold = evaluate_patch_conservation(patch, cons)
    assert threshold == 30.0


def test_empty_patch():
    patch = make_patch(residue_numbers=[])
    passes, pct, threshold = evaluate_patch_conservation(patch, {})
    assert passes is False
    assert pct == 100.0
    assert threshold == 0.0


def test_missing_residue_treated_as_mismatch():
    # Residues not present in cons dict (= no orthologue overlap) count as mismatches
    patch = make_patch(residue_numbers=[1, 2, 3, 4, 5])
    cons = {1: True, 2: True}  # 3, 4, 5 missing → mismatches
    passes, pct, _ = evaluate_patch_conservation(patch, cons)
    assert math.isclose(pct, 60.0, abs_tol=0.01)
    assert passes is False
