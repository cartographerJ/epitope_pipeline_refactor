"""
Shared utility functions for the epitope pipeline.

Logging and the empty-metric default. Structure helpers (extract_ca_coords,
get_chain) live in epitope_pipeline.io.pdb.
"""

import logging


logger = logging.getLogger("epitope_pipeline")


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
