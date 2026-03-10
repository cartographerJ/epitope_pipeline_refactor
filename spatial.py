"""
Step 4: Ectodomain Spatial Filtering — measure distance from the membrane
plane for all residues and isolate extracellular residues >= threshold.

Distances are computed for ALL residues in the structure (not just ECD)
so the visualization shows a continuous distance trace. Qualifying residues
are restricted to extracellular positions meeting the distance threshold.
"""

import logging
from dataclasses import dataclass
from typing import Dict, List

import numpy as np
from Bio.PDB import PDBParser

from epitope_pipeline import config
from epitope_pipeline.utils import get_chain

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class SpatialFilter:
    """Residues passing the distance-from-membrane filter."""
    uniprot_id: str
    gene_name: str
    qualifying_residues: list           # Residue numbers passing the distance threshold
    residue_distances: dict             # {residue_num: distance_from_membrane_A} (ALL residues)
    total_extracellular: int            # Count of all extracellular residues
    total_qualifying: int               # Count of residues >= min_distance
    max_distance: float                 # Maximum distance of any extracellular residue
    min_distance_threshold: float       # Threshold used (Angstroms)
    max_distance_threshold: float = None  # Upper bound (Angstroms), None if distal-only


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def filter_ectodomain(target, structure, membrane, min_distance=None,
                      max_distance=None, ca_coords=None):
    """
    Measure distance from the membrane plane for all residues in the
    structure, then isolate extracellular residues meeting the distance
    threshold.

    Distances are computed for ALL residues (extracellular, TM, and
    cytoplasmic) so the visualization shows a continuous trace. Only
    extracellular residues can qualify for downstream patch analysis.

    Two modes:
      - Distal (default): extracellular residues >= min_distance
      - Proximal: extracellular residues <= max_distance (when max_distance set)

    Args:
        target: TargetInfo object.
        structure: StructureResult with pdb_path and chain_id.
        membrane: MembraneAnnotation with plane definition and topology.
        min_distance: Override min distance threshold (default: config value).
        max_distance: If set, use proximal mode (extracellular residues <= this).
        ca_coords: Optional pre-computed {resnum: np.array} Cα coordinates.

    Returns:
        SpatialFilter with qualifying residue list and distance data.
    """
    if max_distance is None and min_distance is None:
        min_distance = config.ECTODOMAIN_MIN_DISTANCE_A

    # Use pre-computed CA coords or extract from PDB
    if ca_coords is None:
        from epitope_pipeline.utils import extract_ca_coords
        ca_coords = extract_ca_coords(structure.pdb_path, structure.chain_id)

    # Also need to iterate the chain for topology classification
    parser = PDBParser(QUIET=True)
    bio_structure = parser.get_structure(target.gene_name, structure.pdb_path)
    model = bio_structure[0]
    chain = get_chain(model, structure.chain_id)

    residue_distances = {}
    extracellular_residues = []

    for res in chain:
        res_id = res.get_id()
        if res_id[0] != " ":  # Skip hetero atoms / water
            continue
        resnum = res_id[1]

        # Track extracellular residues (for qualifying filter)
        topo = membrane.residue_topology.get(resnum)
        if topo == "extracellular":
            extracellular_residues.append(resnum)

        # Compute distance using pre-computed CA coords
        if resnum not in ca_coords:
            continue

        ca_pos = ca_coords[resnum]
        # Absolute distance along membrane normal from membrane center.
        d = abs(np.dot(ca_pos - membrane.membrane_center, membrane.membrane_normal))
        residue_distances[resnum] = d

    # Re-zero distances to the membrane surface (top of bilayer).
    # raw distance is from membrane center; subtract half_thickness to get
    # distance from the bilayer surface. This is more robust than using the
    # most proximal ECD residue, which can be artifactually buried in
    # AlphaFold models that lack membrane context (e.g., ERBB2, EGFR).
    extracellular_set = set(extracellular_residues)
    surface_offset = membrane.membrane_half_thickness
    for r in residue_distances:
        residue_distances[r] = max(0.0, residue_distances[r] - surface_offset)

    # Filter to EXTRACELLULAR residues meeting distance criteria
    if max_distance is not None:
        # Proximal mode: extracellular residues <= max_distance
        qualifying = [r for r, d in residue_distances.items()
                      if d <= max_distance and r in extracellular_set]
    else:
        # Distal mode (default): extracellular residues >= min_distance
        qualifying = [r for r, d in residue_distances.items()
                      if d >= min_distance and r in extracellular_set]
    qualifying.sort()

    # Max distance of extracellular residues only (for logging)
    ecd_distances = [d for r, d in residue_distances.items() if r in extracellular_set]
    max_d = max(ecd_distances) if ecd_distances else 0.0

    if max_distance is not None:
        logger.info(
            "  %s: %d extracellular residues, %d qualify at <= %.0fA (proximal mode)",
            target.gene_name,
            len(extracellular_residues),
            len(qualifying),
            max_distance,
        )
    else:
        logger.info(
            "  %s: %d extracellular residues, %d qualify at >= %.0fA (max distance: %.1fA)",
            target.gene_name,
            len(extracellular_residues),
            len(qualifying),
            min_distance,
            max_d,
        )

    if not qualifying:
        if max_distance is not None:
            logger.info(
                "  %s: No extracellular residues <= %.0fA from membrane.",
                target.gene_name, max_distance,
            )
        else:
            logger.info(
                "  %s: No residues >= %.0fA from membrane. "
                "Ectodomain extends to %.1fA.",
                target.gene_name, min_distance, max_d,
            )

    return SpatialFilter(
        uniprot_id=target.uniprot_id,
        gene_name=target.gene_name,
        qualifying_residues=qualifying,
        residue_distances=residue_distances,
        total_extracellular=len(extracellular_residues),
        total_qualifying=len(qualifying),
        max_distance=max_d,
        min_distance_threshold=min_distance if min_distance is not None else 0.0,
        max_distance_threshold=max_distance,
    )


