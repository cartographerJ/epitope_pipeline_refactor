"""
Shared utility functions for the epitope pipeline.

Consolidates helpers that were duplicated across run.py, bispecific.py,
spatial.py, and surface.py.
"""

import logging

import numpy as np
from Bio.PDB import PDBParser


logger = logging.getLogger("epitope_pipeline")


# ---------------------------------------------------------------------------
# Structure helpers
# ---------------------------------------------------------------------------

def get_chain(model, chain_id):
    """Get chain by ID with fallback to first chain."""
    if chain_id in model:
        return model[chain_id]
    chains = list(model.get_chains())
    if chains:
        return chains[0]
    raise ValueError("No chains found in structure")


def extract_ca_coords(pdb_path, chain_id):
    """
    Extract Cα coordinates from a PDB file.

    Args:
        pdb_path: Path to PDB file.
        chain_id: Target chain ID.

    Returns:
        Dict {residue_number: np.array([x, y, z])}.
    """
    parser = PDBParser(QUIET=True)
    bio_struct = parser.get_structure("target", pdb_path)
    model = bio_struct[0]
    chain = get_chain(model, chain_id)

    ca_coords = {}
    for res in chain:
        res_id = res.get_id()
        if res_id[0] != " ":
            continue
        if "CA" in res:
            ca_coords[res_id[1]] = res["CA"].get_vector().get_array()

    return ca_coords


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(verbose):
    """Configure root logger for the pipeline."""
    logger.setLevel(logging.DEBUG)
    if verbose:
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(logging.Formatter("%(message)s"))
        # Avoid duplicate handlers on re-runs
        if not logger.handlers:
            logger.addHandler(ch)


# ---------------------------------------------------------------------------
# Empty metric
# ---------------------------------------------------------------------------

def empty_metric():
    """Return an empty metric dict for targets with no epitope space."""
    return {
        "total_epitope_area_a2": 0.0,
        "total_ectodomain_area_a2": 0.0,
        "epitope_fraction": 0.0,
        "n_patches": 0,
        "best_score": 0.0,
    }
