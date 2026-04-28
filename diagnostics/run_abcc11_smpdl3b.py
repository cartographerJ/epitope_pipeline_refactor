#!/usr/bin/env python3
"""
One-off bispecific run: ABCC11 x SMPDL3B

SMPDL3B is GPI-anchored at Ala431 but UniProt lacks the Lipidation/GPI
annotation, so the pipeline's automatic GPI detection fails.

This script patches annotate_membrane to inject the GPI anchor for
SMPDL3B before the normal detection logic runs.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from epitope_pipeline.io import membrane as _membrane
from epitope_pipeline import bispecific as _bispecific

# ── Monkey-patch: inject GPI feature for SMPDL3B ───────────────────
GPI_OVERRIDES = {
    "SMPDL3B": 431,  # Ala431, GPI-anchored (literature)
}

_original_annotate = _membrane.annotate_membrane

def _patched_annotate(target, structure):
    """Inject synthetic GPI Lipidation feature for known missing targets."""
    if target.gene_name in GPI_OVERRIDES:
        target.features.append({
            "type": "Lipidation",
            "description": "GPI-anchor amidated alanine (manual override)",
            "start": GPI_OVERRIDES[target.gene_name],
            "end": GPI_OVERRIDES[target.gene_name],
        })
    return _original_annotate(target, structure)

# Patch in both the module and the bispecific module's local reference
_membrane.annotate_membrane = _patched_annotate
_bispecific.annotate_membrane = _patched_annotate
# ────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    results = _bispecific.run_bispecific(
        [("ABCC11", "SMPDL3B")],
        run_name="bispecific_abcc11_smpdl3b",
    )
