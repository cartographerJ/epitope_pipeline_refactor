# tests/test_compute_conservation_dispatch.py
"""analyze_conservation dispatches between local and whole_patch modes."""

from types import SimpleNamespace

import numpy as np

from epitope_pipeline.compute.conservation import analyze_conservation
from tests._fixtures import make_patch


def _build_target_and_surface():
    """30-residue patch: residues 1-4 are clustered mismatches, 5-30 conserved.

    Sequences are equal length and differ only at positions 1-4 (1:1 alignment),
    so residue_conservation maps cleanly: 1-4 False, 5-30 True.
    """
    human = "A" * 30
    cyno = "C" * 4 + "A" * 26  # positions 1-4 differ
    target = SimpleNamespace(
        uniprot_id="TEST", gene_name="TEST",
        sequence=human, cyno_sequence=cyno,
    )

    residues = list(range(1, 31))
    patch = make_patch(patch_id=0, residue_numbers=residues,
                       residue_aas={r: "A" for r in residues})
    surface = SimpleNamespace(
        patches=[patch],
        residue_sasa={r: 30.0 for r in residues},  # 30 * 30 = 900 Å² total
    )

    # Coords: residues 1-4 tightly clustered at the origin; 5-30 a far, contiguous
    # chain spaced 5 Å (so they re-cluster into one sub-patch).
    ca = {1: (0, 0, 0), 2: (2, 0, 0), 3: (4, 0, 0), 4: (6, 0, 0)}
    for i, r in enumerate(range(5, 31)):
        ca[r] = (100.0 + i * 5.0, 0, 0)
    ca = {r: np.array(c, dtype=float) for r, c in ca.items()}
    return target, surface, ca


def test_whole_patch_keeps_full_patch():
    target, surface, ca = _build_target_and_surface()
    res = analyze_conservation(target, surface, ca, cyno_mode="whole_patch")
    # 30 residues, 4 mismatches = 13.3% <= scaled 18.4% threshold -> full patch kept
    assert len(res.conserved_patches) == 1
    assert len(res.conserved_patches[0].residue_numbers) == 30


def test_local_trims_mismatch_dense_residues():
    target, surface, ca = _build_target_and_surface()
    res = analyze_conservation(target, surface, ca, cyno_mode="local",
                               max_mismatches_per_window=2)
    # Local gate trims clustered mismatches 1-4, salvages residues 5-30
    assert len(res.conserved_patches) == 1
    assert set(res.conserved_patches[0].residue_numbers) == set(range(5, 31))
