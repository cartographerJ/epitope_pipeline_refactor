"""
Step 5: Surface Patch Identification — SASA calculation and spatial clustering.

Calculates per-residue solvent-accessible surface area (SASA) using FreeSASA,
identifies surface-exposed residues among the spatially qualified set, and
clusters them into contiguous surface patches using single-linkage
agglomerative clustering on Cα-Cα distances.

Patches below the VHH CDR footprint minimum (600 A²) are discarded.
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
from Bio.PDB import PDBParser

from epitope_pipeline.config import (
    AA_1TO3,
    MAX_SASA,
    PATCH_CLUSTERING_DISTANCE_A,
    SASA_PROBE_RADIUS,
    SURFACE_EXPOSURE_THRESHOLD,
    VHH_FOOTPRINT_MIN_A2,
)
from epitope_pipeline.io.pdb import get_chain, extract_ca_coords

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class SurfacePatch:
    """A contiguous patch of surface-exposed residues."""
    patch_id: int
    residue_numbers: list               # List of residue numbers in this patch
    residue_aas: dict                   # {resnum: one_letter_aa}
    total_sasa_a2: float                # Sum of SASA for patch residues (A²)
    centroid: np.ndarray                # 3D centroid of patch Cα atoms
    max_dimension_a: float              # Largest pairwise Cα distance in patch
    avg_distance_from_membrane: float   # Mean distance along membrane normal


@dataclass
class SurfaceAnalysis:
    """Complete surface analysis for one target."""
    uniprot_id: str
    gene_name: str
    residue_sasa: dict                  # {resnum: sasa_A2}
    residue_relative_sasa: dict         # {resnum: sasa / max_sasa_for_aa}
    exposed_residues: list              # Residue numbers with SASA > threshold
    patches: list                       # List of SurfacePatch objects
    total_epitope_surface_a2: float     # Sum of qualifying patch areas


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def analyze_surface(target, structure, spatial_filter, membrane=None,
                     ca_coords=None):
    """
    Identify surface-exposed patches in the spatially qualified ectodomain.

    Steps:
      1. Calculate per-residue SASA using FreeSASA
      2. Filter to qualifying residues from SpatialFilter
      3. Identify surface-exposed residues (relative SASA > threshold)
      4. Cluster into patches using single-linkage agglomerative clustering
      5. Filter patches by minimum VHH footprint area (600 A²)

    Args:
        target: TargetInfo object.
        structure: StructureResult with pdb_path and chain_id.
        spatial_filter: SpatialFilter from the spatial module.
        membrane: Optional MembraneAnnotation for distance calculations.
        ca_coords: Optional pre-computed {resnum: np.array} Cα coordinates.

    Returns:
        SurfaceAnalysis with patches and SASA data.
    """
    # Step 1: Always calculate SASA (needed for visualization even when
    # no residues pass the spatial distance filter)
    residue_sasa = _calculate_sasa(structure.pdb_path, structure.chain_id)
    logger.info("  Calculated SASA for %d residues", len(residue_sasa))

    if not spatial_filter.qualifying_residues:
        logger.info("  %s: No qualifying residues — skipping patch identification", target.gene_name)
        return SurfaceAnalysis(
            uniprot_id=target.uniprot_id,
            gene_name=target.gene_name,
            residue_sasa=residue_sasa,
            residue_relative_sasa={},
            exposed_residues=[],
            patches=[],
            total_epitope_surface_a2=0.0,
        )

    # Parse structure for residue names
    parser = PDBParser(QUIET=True)
    bio_structure = parser.get_structure(target.gene_name, structure.pdb_path)
    model = bio_structure[0]
    chain = get_chain(model, structure.chain_id)

    # Build amino acid mapping for the qualifying residues
    residue_aa = {}
    for res in chain:
        res_id = res.get_id()
        if res_id[0] != " ":
            continue
        resnum = res_id[1]
        resname = res.get_resname()
        residue_aa[resnum] = resname

    # Use passed-in CA coords or extract from structure
    if ca_coords is None:
        ca_coords = extract_ca_coords(structure.pdb_path, structure.chain_id)

    # Step 2: Compute relative SASA for qualifying residues
    qualifying_set = set(spatial_filter.qualifying_residues)
    residue_relative_sasa = {}
    for resnum in qualifying_set:
        sasa_val = residue_sasa.get(resnum, 0.0)
        resname = residue_aa.get(resnum, "ALA")
        max_sasa = MAX_SASA.get(resname, 150.0)
        residue_relative_sasa[resnum] = sasa_val / max_sasa if max_sasa > 0 else 0.0

    # Step 3: Identify exposed residues
    exposed = [
        r for r in qualifying_set
        if residue_relative_sasa.get(r, 0.0) > SURFACE_EXPOSURE_THRESHOLD
    ]
    exposed.sort()
    logger.info(
        "  %d of %d qualifying residues are surface-exposed (rel SASA > %.0f%%)",
        len(exposed), len(qualifying_set), SURFACE_EXPOSURE_THRESHOLD * 100,
    )

    if not exposed:
        return SurfaceAnalysis(
            uniprot_id=target.uniprot_id,
            gene_name=target.gene_name,
            residue_sasa=residue_sasa,
            residue_relative_sasa=residue_relative_sasa,
            exposed_residues=[],
            patches=[],
            total_epitope_surface_a2=0.0,
        )

    # Only cluster residues that have Cα coordinates
    exposed_with_ca = [r for r in exposed if r in ca_coords]

    # Step 4: Cluster into patches
    clusters = cluster_surface_patches(exposed_with_ca, ca_coords)

    # Step 5: Build SurfacePatch objects and filter by area
    # Build one-letter AA map
    aa_1letter = {}
    for resnum, resname in residue_aa.items():
        for one, three in AA_1TO3.items():
            if three == resname:
                aa_1letter[resnum] = one
                break

    patches = []
    patch_id = 0
    for cluster_residues in clusters:
        total_area = sum(residue_sasa.get(r, 0.0) for r in cluster_residues)

        if total_area < VHH_FOOTPRINT_MIN_A2:
            continue

        # Compute patch properties
        coords = np.array([ca_coords[r] for r in cluster_residues])
        centroid = np.mean(coords, axis=0)

        # Max pairwise dimension
        if len(coords) > 1:
            from scipy.spatial.distance import pdist
            max_dim = np.max(pdist(coords))
        else:
            max_dim = 0.0

        # Average distance from membrane (if available)
        avg_dist = 0.0
        if spatial_filter.residue_distances:
            dists = [spatial_filter.residue_distances.get(r, 0.0) for r in cluster_residues]
            avg_dist = np.mean(dists) if dists else 0.0

        patch = SurfacePatch(
            patch_id=patch_id,
            residue_numbers=sorted(cluster_residues),
            residue_aas={r: aa_1letter.get(r, "X") for r in cluster_residues},
            total_sasa_a2=total_area,
            centroid=centroid,
            max_dimension_a=max_dim,
            avg_distance_from_membrane=avg_dist,
        )
        patches.append(patch)
        patch_id += 1

    total_surface = sum(p.total_sasa_a2 for p in patches)

    logger.info(
        "  %d surface patches (>= %.0f A²), total epitope surface: %.0f A²",
        len(patches), VHH_FOOTPRINT_MIN_A2, total_surface,
    )
    for p in patches:
        logger.info(
            "    Patch %d: %d residues, %.0f A², max dim %.1fA, avg dist %.1fA",
            p.patch_id, len(p.residue_numbers), p.total_sasa_a2,
            p.max_dimension_a, p.avg_distance_from_membrane,
        )

    return SurfaceAnalysis(
        uniprot_id=target.uniprot_id,
        gene_name=target.gene_name,
        residue_sasa=residue_sasa,
        residue_relative_sasa=residue_relative_sasa,
        exposed_residues=exposed,
        patches=patches,
        total_epitope_surface_a2=total_surface,
    )


def cluster_ectodomain_patches(target, structure, membrane, residue_sasa=None,
                               ca_coords=None):
    """
    Cluster ALL surface-exposed ectodomain residues into 3D patches,
    ignoring distance and conservation filters.

    Used for whole-protein specificity screening: each patch >= 600 A²
    is checked against BLAST for off-target cross-reactivity risk.

    Args:
        target: TargetInfo object.
        structure: StructureResult with pdb_path and chain_id.
        membrane: MembraneAnnotation with residue_topology.
        residue_sasa: Optional pre-computed {resnum: sasa_A2}.
        ca_coords: Optional pre-computed {resnum: np.array} Cα coordinates.

    Returns:
        List of SurfacePatch objects (area >= VHH_FOOTPRINT_MIN_A2).
    """
    if residue_sasa is None:
        residue_sasa = _calculate_sasa(structure.pdb_path, structure.chain_id)

    # All ectodomain residues from membrane topology
    ec_residues = set()
    for resnum, topo in membrane.residue_topology.items():
        if topo == "extracellular":
            ec_residues.add(resnum)

    if not ec_residues:
        return []

    # Use passed-in CA coords or extract from structure
    if ca_coords is None:
        ca_coords = extract_ca_coords(structure.pdb_path, structure.chain_id)

    # Parse structure for residue names
    parser = PDBParser(QUIET=True)
    bio_structure = parser.get_structure(target.gene_name, structure.pdb_path)
    model = bio_structure[0]
    chain = get_chain(model, structure.chain_id)

    residue_aa = {}
    for res in chain:
        res_id = res.get_id()
        if res_id[0] != " ":
            continue
        resnum = res_id[1]
        residue_aa[resnum] = res.get_resname()

    # Filter to surface-exposed ectodomain residues (must also have CA coords)
    exposed = []
    for resnum in ec_residues:
        sasa_val = residue_sasa.get(resnum, 0.0)
        resname = residue_aa.get(resnum, "ALA")
        max_sasa = MAX_SASA.get(resname, 150.0)
        rel_sasa = sasa_val / max_sasa if max_sasa > 0 else 0.0
        if rel_sasa > SURFACE_EXPOSURE_THRESHOLD and resnum in ca_coords:
            exposed.append(resnum)

    if len(exposed) < 2:
        return []

    clusters = cluster_surface_patches(exposed, ca_coords)

    # Build SurfacePatch objects for patches >= 600 A²
    aa_1letter = {}
    for resnum, resname in residue_aa.items():
        for one, three in AA_1TO3.items():
            if three == resname:
                aa_1letter[resnum] = one
                break

    patches = []
    patch_id = 0
    for cluster_residues in clusters:
        total_area = sum(residue_sasa.get(r, 0.0) for r in cluster_residues)
        if total_area < VHH_FOOTPRINT_MIN_A2:
            continue

        coords = np.array([ca_coords[r] for r in cluster_residues])
        centroid = np.mean(coords, axis=0)

        if len(coords) > 1:
            from scipy.spatial.distance import pdist
            max_dim = np.max(pdist(coords))
        else:
            max_dim = 0.0

        patch = SurfacePatch(
            patch_id=patch_id,
            residue_numbers=sorted(cluster_residues),
            residue_aas={r: aa_1letter.get(r, "X") for r in cluster_residues},
            total_sasa_a2=total_area,
            centroid=centroid,
            max_dimension_a=max_dim,
            avg_distance_from_membrane=0.0,
        )
        patches.append(patch)
        patch_id += 1

    logger.info(
        "  Ectodomain surface: %d exposed residues -> %d patches >= %.0f A²",
        len(exposed), len(patches), VHH_FOOTPRINT_MIN_A2,
    )

    return patches


# ---------------------------------------------------------------------------
# SASA calculation
# ---------------------------------------------------------------------------

def _calculate_sasa(pdb_path, chain_id):
    """
    Calculate per-residue SASA using FreeSASA.

    Args:
        pdb_path: Path to PDB file.
        chain_id: Target chain ID.

    Returns:
        Dict {residue_number: sasa_A2}.
    """
    import freesasa

    sasa_structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(sasa_structure)

    # residueAreas() returns {chain_id: {resnum_str: ResidueArea}}
    residue_sasa = {}
    all_areas = result.residueAreas()

    # Find our target chain (try exact match, then first available)
    chain_areas = all_areas.get(chain_id)
    if chain_areas is None and len(all_areas) == 1:
        chain_areas = next(iter(all_areas.values()))
    if chain_areas is None:
        logger.warning("  Chain '%s' not found in FreeSASA output (chains: %s)",
                       chain_id, list(all_areas.keys()))
        return residue_sasa

    for resnum_str, area in chain_areas.items():
        try:
            resnum = int(resnum_str)
        except ValueError:
            continue
        residue_sasa[resnum] = area.total

    return residue_sasa


# ---------------------------------------------------------------------------
# Surface patch clustering
# ---------------------------------------------------------------------------

def cluster_surface_patches(exposed_residues, ca_coords, max_distance=None):
    """
    Cluster surface-exposed residues into contiguous patches using
    single-linkage agglomerative clustering on Cα-Cα distances.

    Single-linkage ensures that patches represent connected surface regions:
    two residues are in the same patch if there exists a chain of nearby
    surface-exposed residues connecting them, each pair within max_distance.

    Args:
        exposed_residues: Sorted list of residue numbers.
        ca_coords: Dict {resnum: np.array([x, y, z])}.
        max_distance: Max Cα-Cα distance for connectivity (default: 8.0 A).

    Returns:
        List of lists, each inner list being residue numbers in one cluster.
    """
    if max_distance is None:
        max_distance = PATCH_CLUSTERING_DISTANCE_A

    if len(exposed_residues) < 2:
        return [exposed_residues] if exposed_residues else []

    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import fcluster, linkage

    # Build coordinate matrix in residue order
    residue_list = sorted(exposed_residues)
    coords = np.array([ca_coords[r] for r in residue_list])

    # Pairwise Euclidean distances between Cα atoms
    dist_matrix = pdist(coords)

    # Single-linkage agglomerative clustering
    Z = linkage(dist_matrix, method="single")

    # Cut the dendrogram at max_distance threshold
    labels = fcluster(Z, t=max_distance, criterion="distance")

    # Group residues by cluster label
    clusters = {}
    for res, label in zip(residue_list, labels):
        clusters.setdefault(label, []).append(res)

    return list(clusters.values())


