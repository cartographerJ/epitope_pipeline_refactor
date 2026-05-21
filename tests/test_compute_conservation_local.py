# tests/test_compute_conservation_local.py
"""Unit tests for the restored local sliding-window cyno helpers."""

import numpy as np

from epitope_pipeline.compute.conservation import (
    check_local_mismatch_density,
    identify_failing_residues,
    _check_conserved_subpatch,
    _recluster_survivors,
)
from tests._fixtures import make_patch


def test_density_counts_clustered_mismatches():
    coords = {1: (0, 0, 0), 2: (3, 0, 0), 3: (6, 0, 0), 4: (9, 0, 0), 5: (12, 0, 0)}
    coords = {r: np.array(c, dtype=float) for r, c in coords.items()}
    cons = {1: False, 2: False, 3: False, 4: True, 5: True}
    worst_res, worst_count = check_local_mismatch_density([1, 2, 3, 4, 5], cons, coords)
    assert worst_count == 3
    assert worst_res == 1  # first residue to accumulate the worst neighborhood


def test_density_all_conserved_is_zero():
    coords = {r: np.array((r * 3.0, 0, 0)) for r in range(1, 6)}
    cons = {r: True for r in range(1, 6)}
    assert check_local_mismatch_density([1, 2, 3, 4, 5], cons, coords) == (0, 0)


def test_density_no_coords_returns_zero():
    # Residues absent from ca_coords -> cannot run the spatial test
    assert check_local_mismatch_density([1, 2, 3], {1: False}, {}) == (0, 0)


def test_identify_failing_no_coords_returns_empty():
    assert identify_failing_residues([1, 2, 3], {1: False}, {}, 2) == set()


def test_conserved_subpatch_empty_returns_false():
    assert _check_conserved_subpatch([], {}, {}) is False


def test_identify_failing_residues_local_cluster():
    coords = {
        1: (0, 0, 0), 2: (2, 0, 0), 3: (4, 0, 0),
        10: (100, 0, 0), 11: (102, 0, 0), 12: (104, 0, 0),
    }
    coords = {r: np.array(c, dtype=float) for r, c in coords.items()}
    cons = {1: False, 2: False, 3: False, 10: True, 11: True, 12: True}
    failing = identify_failing_residues([1, 2, 3, 10, 11, 12], cons, coords, 2)
    assert failing == {1, 2, 3}


def test_identify_failing_none_when_threshold_high():
    coords = {1: (0, 0, 0), 2: (2, 0, 0), 3: (4, 0, 0)}
    coords = {r: np.array(c, dtype=float) for r, c in coords.items()}
    cons = {1: False, 2: False, 3: False}
    assert identify_failing_residues([1, 2, 3], cons, coords, 5) == set()


def test_conserved_subpatch_contiguous_passes():
    coords = {r: np.array((i * 5.0, 0, 0)) for i, r in enumerate(range(1, 8))}
    sasa = {r: 100.0 for r in range(1, 8)}
    assert _check_conserved_subpatch(list(range(1, 8)), coords, sasa) is True


def test_conserved_subpatch_sparse_fails():
    coords = {1: np.array((0.0, 0, 0)), 2: np.array((100.0, 0, 0))}
    sasa = {1: 100.0, 2: 100.0}
    assert _check_conserved_subpatch([1, 2], coords, sasa) is False


def test_recluster_survivors_builds_subpatch():
    survivors = list(range(1, 8))
    coords = {r: np.array((i * 5.0, 0, 0)) for i, r in enumerate(survivors)}
    sasa = {r: 100.0 for r in survivors}
    cons = {r: True for r in survivors}
    original = make_patch(patch_id=0, residue_numbers=survivors,
                          residue_aas={r: "A" for r in survivors})
    subs = _recluster_survivors(survivors, coords, sasa, original, cons, next_patch_id=1)
    assert len(subs) == 1
    assert subs[0].patch_id == 1
    assert sorted(subs[0].residue_numbers) == survivors
    assert subs[0].total_sasa_a2 == 700.0
