"""Tests for the per-paralog specificity rule, patch merging, and segment extraction."""

import math

from epitope_pipeline import config
from epitope_pipeline.compute.specificity import (
    evaluate_patch_specificity,
    merge_adjacent_patches,
    _extract_patch_sequences,
)
from tests._fixtures import make_patch


# ---------------------------------------------------------------------------
# evaluate_patch_specificity
# ---------------------------------------------------------------------------

def test_no_paralog_matches_passes():
    patch = make_patch(residue_numbers=[1, 2, 3, 4, 5])
    passes, pct, threshold, worst = evaluate_patch_specificity(patch, {})
    assert passes is True
    assert pct == 0.0
    assert worst is None
    assert math.isclose(threshold, config.MAX_NONSPECIFIC_PERCENT, rel_tol=1e-6)


def test_single_paralog_below_threshold_passes():
    # 5-residue patch, paralog matches 2 of 5 = 40% < 85% threshold
    patch = make_patch(residue_numbers=[1, 2, 3, 4, 5])
    paralogs = {"P00001": {1, 2}}
    passes, pct, _, worst = evaluate_patch_specificity(patch, paralogs)
    assert passes is True
    assert math.isclose(pct, 40.0, abs_tol=0.01)
    assert worst == "P00001"


def test_single_paralog_above_threshold_fails():
    # 5-residue patch, paralog matches all 5 = 100% > 85%
    patch = make_patch(residue_numbers=[1, 2, 3, 4, 5])
    paralogs = {"P00001": {1, 2, 3, 4, 5}}
    passes, pct, _, worst = evaluate_patch_specificity(patch, paralogs)
    assert passes is False
    assert pct == 100.0
    assert worst == "P00001"


def test_max_over_paralogs_wins():
    # Two paralogs: weak (20%) and strong (60%). Worst = strong.
    patch = make_patch(residue_numbers=[1, 2, 3, 4, 5])
    paralogs = {
        "WEAK": {1},                   # 20%
        "STRONG": {1, 2, 3},            # 60%
    }
    _, pct, _, worst = evaluate_patch_specificity(patch, paralogs)
    assert math.isclose(pct, 60.0, abs_tol=0.01)
    assert worst == "STRONG"


def test_paralog_outside_patch_ignored():
    # Paralog matches 100 positions, but only 1 inside the patch → 20% of patch
    patch = make_patch(residue_numbers=[1, 2, 3, 4, 5])
    paralogs = {"P00001": set(range(1, 101))}
    _, pct, _, _ = evaluate_patch_specificity(patch, paralogs)
    assert math.isclose(pct, 100.0, abs_tol=0.01)  # all 5 covered


def test_empty_patch():
    patch = make_patch(residue_numbers=[])
    passes, pct, _, worst = evaluate_patch_specificity(patch, {"X": {1, 2}})
    assert passes is True
    assert pct == 0.0
    assert worst is None


# ---------------------------------------------------------------------------
# merge_adjacent_patches
# ---------------------------------------------------------------------------

def test_merge_empty_returns_empty():
    assert merge_adjacent_patches([]) == []


def test_merge_single_returns_unchanged():
    p = make_patch(centroid=(0, 0, 0))
    result = merge_adjacent_patches([p])
    assert len(result) == 1
    assert result[0] is p


def test_merge_far_patches_unchanged():
    # Centroids 100 A apart, threshold 15 A
    p1 = make_patch(patch_id=0, centroid=(0, 0, 0), residue_numbers=[1, 2, 3])
    p2 = make_patch(patch_id=1, centroid=(100, 0, 0), residue_numbers=[50, 51, 52])
    result = merge_adjacent_patches([p1, p2], distance_threshold_a=15.0)
    assert len(result) == 2


def test_merge_close_patches_combined():
    # Centroids 10 A apart, threshold 15 → merge
    p1 = make_patch(patch_id=0, centroid=(0, 0, 0),
                    residue_numbers=[1, 2, 3], total_sasa_a2=500.0,
                    avg_distance_from_membrane=80.0)
    p2 = make_patch(patch_id=1, centroid=(10, 0, 0),
                    residue_numbers=[50, 51], total_sasa_a2=400.0,
                    avg_distance_from_membrane=90.0)
    result = merge_adjacent_patches([p1, p2], distance_threshold_a=15.0)
    assert len(result) == 1
    merged = result[0]
    assert merged.residue_numbers == [1, 2, 3, 50, 51]
    assert merged.total_sasa_a2 == 900.0
    # Avg distance averaged across input patches
    assert math.isclose(merged.avg_distance_from_membrane, 85.0, rel_tol=1e-6)


def test_merge_transitive_chain():
    # A-B close, B-C close, A-C far → all three merge via transitivity
    a = make_patch(patch_id=0, centroid=(0, 0, 0), residue_numbers=[1])
    b = make_patch(patch_id=1, centroid=(10, 0, 0), residue_numbers=[2])
    c = make_patch(patch_id=2, centroid=(20, 0, 0), residue_numbers=[3])
    result = merge_adjacent_patches([a, b, c], distance_threshold_a=15.0)
    assert len(result) == 1
    assert result[0].residue_numbers == [1, 2, 3]


# ---------------------------------------------------------------------------
# _extract_patch_sequences
# ---------------------------------------------------------------------------

def test_extract_contiguous_run():
    seq = "ABCDEFGHIJ"
    segments = _extract_patch_sequences(seq, [3, 4, 5, 6])
    assert segments == [("CDEF", 3, 6)]


def test_extract_split_at_large_gap():
    # Gap of 30 positions > 20 → split into two segments
    seq = "X" * 100
    segments = _extract_patch_sequences(seq, [5, 6, 7, 50, 51, 52])
    assert len(segments) == 2
    assert segments[0][1:] == (5, 7)
    assert segments[1][1:] == (50, 52)


def test_small_gap_not_split():
    # Gap of 10 positions ≤ 20 → kept as one segment spanning the gap
    seq = "X" * 100
    segments = _extract_patch_sequences(seq, [5, 6, 7, 18, 19, 20])
    assert len(segments) == 1
    assert segments[0][1:] == (5, 20)


def test_extract_empty():
    assert _extract_patch_sequences("ABC", []) == []
