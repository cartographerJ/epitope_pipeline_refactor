"""Tests for single-linkage surface-patch clustering."""

import numpy as np

from epitope_pipeline.compute.surface import cluster_surface_patches


def _coords(positions):
    """Build a {resnum: np.array([x, y, z])} dict from a list of (resnum, xyz)."""
    return {r: np.array(xyz, dtype=float) for r, xyz in positions}


def test_empty_input():
    assert cluster_surface_patches([], {}) == []


def test_single_residue():
    coords = _coords([(1, (0, 0, 0))])
    clusters = cluster_surface_patches([1], coords)
    assert clusters == [[1]]


def test_two_close_residues_one_cluster():
    coords = _coords([(1, (0, 0, 0)), (2, (5, 0, 0))])  # 5 A apart
    clusters = cluster_surface_patches([1, 2], coords, max_distance=8.0)
    assert len(clusters) == 1
    assert sorted(clusters[0]) == [1, 2]


def test_two_far_residues_two_clusters():
    coords = _coords([(1, (0, 0, 0)), (2, (50, 0, 0))])
    clusters = cluster_surface_patches([1, 2], coords, max_distance=8.0)
    assert len(clusters) == 2


def test_chain_single_linkage():
    # 1—2—3—4 each 5 A apart, max_distance 8 → all in one cluster
    coords = _coords([
        (1, (0, 0, 0)), (2, (5, 0, 0)),
        (3, (10, 0, 0)), (4, (15, 0, 0)),
    ])
    clusters = cluster_surface_patches([1, 2, 3, 4], coords, max_distance=8.0)
    assert len(clusters) == 1
    assert sorted(clusters[0]) == [1, 2, 3, 4]


def test_two_separate_clusters():
    coords = _coords([
        (1, (0, 0, 0)), (2, (5, 0, 0)),       # cluster A
        (10, (50, 0, 0)), (11, (55, 0, 0)),   # cluster B (50 A away)
    ])
    clusters = cluster_surface_patches([1, 2, 10, 11], coords, max_distance=8.0)
    assert len(clusters) == 2
    cluster_sets = sorted(sorted(c) for c in clusters)
    assert cluster_sets == [[1, 2], [10, 11]]
