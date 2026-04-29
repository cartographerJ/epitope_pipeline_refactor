"""Lightweight fixture builders for compute-layer unit tests.

These avoid running the full pipeline by constructing minimal SurfacePatch /
SurfaceAnalysis / EpitopeScore objects with just the attributes each
function under test actually reads.
"""

from types import SimpleNamespace

import numpy as np

from epitope_pipeline.compute.surface import SurfacePatch


def make_patch(
    patch_id=0,
    residue_numbers=None,
    total_sasa_a2=800.0,
    centroid=(0.0, 0.0, 0.0),
    avg_distance_from_membrane=80.0,
    max_dimension_a=20.0,
    residue_aas=None,
):
    """Build a minimal SurfacePatch. Defaults keep tests terse."""
    if residue_numbers is None:
        residue_numbers = [10, 11, 12, 13, 14]
    if residue_aas is None:
        residue_aas = {r: "A" for r in residue_numbers}
    return SurfacePatch(
        patch_id=patch_id,
        residue_numbers=list(residue_numbers),
        residue_aas=dict(residue_aas),
        total_sasa_a2=float(total_sasa_a2),
        centroid=np.array(centroid, dtype=float),
        max_dimension_a=float(max_dimension_a),
        avg_distance_from_membrane=float(avg_distance_from_membrane),
    )


def make_score(
    patch_area_a2=800.0,
    cyno_identity=0.95,
    n_residues=20,
    composite_score=0.6,
):
    """Build a duck-typed EpitopeScore stand-in for compute_target_epitope_metric."""
    return SimpleNamespace(
        patch_area_a2=patch_area_a2,
        cyno_identity=cyno_identity,
        n_residues=n_residues,
        composite_score=composite_score,
    )


def make_surface_analysis(residue_sasa=None, exposed_residues=None):
    """Build a duck-typed SurfaceAnalysis stand-in."""
    if residue_sasa is None:
        residue_sasa = {r: 100.0 for r in range(10, 20)}
    if exposed_residues is None:
        exposed_residues = list(residue_sasa.keys())
    return SimpleNamespace(
        residue_sasa=dict(residue_sasa),
        exposed_residues=list(exposed_residues),
    )
