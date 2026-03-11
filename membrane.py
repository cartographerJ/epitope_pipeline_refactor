"""
Step 3: Membrane Annotation — define the membrane plane and classify residues.

Primary approach: look up the target in the OPM (Orientations of Proteins
in Membranes) database, which provides a pre-oriented PDB with DUM atoms
marking the membrane boundaries.

Fallback: calculate the membrane plane from UniProt transmembrane annotations
mapped onto 3D coordinates using PCA on TM helix Cα atoms.

Every residue is classified as extracellular, transmembrane, or intracellular.
"""

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import requests

from Bio.PDB import PDBParser, PDBIO, Select

from epitope_pipeline.config import (
    CACHE_DIR,
    MEMBRANE_SLAB_HALF_THICKNESS,
    OPM_API,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class MembraneAnnotation:
    """Membrane plane definition and per-residue topology classification."""
    uniprot_id: str
    gene_name: str
    membrane_normal: np.ndarray         # Unit vector perpendicular to membrane
    membrane_center: np.ndarray         # Point on the membrane midplane (3D)
    membrane_half_thickness: float      # Half-thickness in Angstroms
    topology_type: str                  # "single_pass" or "multi_pass"
    tm_segments: list                   # [(start_res, end_res), ...] residue ranges
    residue_topology: dict              # {residue_num: "extracellular"|"transmembrane"|"intracellular"}
    source: str                         # "opm", "uniprot_calculated"


class MembraneAnnotationError(Exception):
    """Could not determine membrane orientation."""
    pass


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def annotate_membrane(target, structure):
    """
    Define the membrane plane and classify residues for a target protein.

    Tries OPM database first for pre-calculated membrane orientation.
    Falls back to computing the membrane plane from UniProt TM annotations
    and 3D coordinates.

    Args:
        target: TargetInfo with features (TM annotations) and sequence.
        structure: StructureResult with pdb_path and chain_id.

    Returns:
        MembraneAnnotation with membrane plane and residue topology.

    Raises:
        MembraneAnnotationError: If no TM annotations are available.
    """
    # Extract TM segments from UniProt features
    tm_segments = _extract_tm_segments(target.features)

    # Check for GPI-anchored protein (no TM but has Lipidation/GPI-anchor)
    if not tm_segments:
        gpi_residue = _extract_gpi_anchor(target.features)
        if gpi_residue is not None:
            return _annotate_gpi_anchored(target, structure, gpi_residue)
        raise MembraneAnnotationError(
            "No transmembrane or GPI-anchor annotations found in UniProt for {}. "
            "Is this a membrane protein?".format(target.gene_name)
        )

    topology_type = "single_pass" if len(tm_segments) == 1 else "multi_pass"
    logger.info(
        "  %s: %s (%d TM segments: %s)",
        target.gene_name, topology_type, len(tm_segments),
        ", ".join("{}-{}".format(s, e) for s, e in tm_segments),
    )

    # Parse the PDB structure
    parser = PDBParser(QUIET=True)
    bio_structure = parser.get_structure(target.gene_name, structure.pdb_path)
    model = bio_structure[0]
    chain = _get_chain(model, structure.chain_id)

    # Try OPM first (only for experimental structures with a PDB ID)
    if structure.source == "pdb_experimental":
        opm_result = _try_opm_lookup(structure.pdb_id, structure.chain_id)
        if opm_result is not None:
            normal, center, half_thick = opm_result
            residue_topo = _classify_residues_from_plane(
                chain, normal, center, half_thick,
                target.features, tm_segments,
            )
            logger.info("  Membrane plane from OPM database")
            return MembraneAnnotation(
                uniprot_id=target.uniprot_id,
                gene_name=target.gene_name,
                membrane_normal=normal,
                membrane_center=center,
                membrane_half_thickness=half_thick,
                topology_type=topology_type,
                tm_segments=tm_segments,
                residue_topology=residue_topo,
                source="opm",
            )

    # Fallback: calculate from UniProt TM annotations + 3D coordinates
    logger.info("  Computing membrane plane from UniProt TM annotations...")
    normal, center = _calculate_membrane_plane(
        chain, tm_segments, topology_type, target.features,
    )
    half_thick = MEMBRANE_SLAB_HALF_THICKNESS

    # Orient normal so extracellular side is positive
    normal = _orient_normal(chain, normal, center, target.features)

    residue_topo = _classify_residues_from_plane(
        chain, normal, center, half_thick,
        target.features, tm_segments,
    )

    return MembraneAnnotation(
        uniprot_id=target.uniprot_id,
        gene_name=target.gene_name,
        membrane_normal=normal,
        membrane_center=center,
        membrane_half_thickness=half_thick,
        topology_type=topology_type,
        tm_segments=tm_segments,
        residue_topology=residue_topo,
        source="uniprot_calculated",
    )


# ---------------------------------------------------------------------------
# TM segment extraction from UniProt features
# ---------------------------------------------------------------------------

def _extract_tm_segments(features):
    """
    Extract transmembrane segment residue ranges from UniProt features.

    Returns list of (start, end) tuples (1-based residue numbers).
    """
    segments = []
    for feat in features:
        if feat["type"] == "Transmembrane":
            start = feat.get("start")
            end = feat.get("end")
            if start is not None and end is not None:
                segments.append((int(start), int(end)))
    segments.sort(key=lambda s: s[0])
    return segments


def _extract_gpi_anchor(features):
    """
    Detect GPI-anchor from UniProt Lipidation features.

    Returns the GPI-anchor residue number (1-based), or None.
    """
    for feat in features:
        if feat["type"] == "Lipidation":
            desc = feat.get("description", "").lower()
            if "gpi" in desc:
                start = feat.get("start")
                if start is not None:
                    return int(start)
    return None


def _annotate_gpi_anchored(target, structure, gpi_residue):
    """
    Define membrane plane for a GPI-anchored protein.

    GPI-anchored proteins have no TM helix. The entire extracellular domain
    is tethered to the membrane via a C-terminal glycolipid anchor. The
    membrane plane is placed at the anchor point, with the normal oriented
    so the bulk of the protein is on the positive (extracellular) side.

    Args:
        target: TargetInfo object.
        structure: StructureResult with pdb_path and chain_id.
        gpi_residue: Residue number of the GPI-anchor attachment.

    Returns:
        MembraneAnnotation for the GPI-anchored protein.
    """
    logger.info("  %s: GPI-anchored protein (anchor at residue %d)",
                target.gene_name, gpi_residue)

    parser = PDBParser(QUIET=True)
    bio_structure = parser.get_structure(target.gene_name, structure.pdb_path)
    model = bio_structure[0]
    chain = _get_chain(model, structure.chain_id)

    # Collect all Cα coordinates
    all_coords = []
    all_resnums = []
    for res in chain:
        res_id = res.get_id()
        if res_id[0] != " ":
            continue
        if "CA" in res:
            all_coords.append(res["CA"].get_vector().get_array())
            all_resnums.append(res_id[1])

    if len(all_coords) < 10:
        raise MembraneAnnotationError(
            "Too few resolved residues ({}) for GPI-anchored protein".format(
                len(all_coords)
            )
        )

    all_coords = np.array(all_coords)
    all_resnums = np.array(all_resnums)

    # PCA: principal axis of the ectodomain = membrane normal
    centroid = np.mean(all_coords, axis=0)
    centered = all_coords - centroid
    cov = np.cov(centered.T)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    normal = eigenvectors[:, -1]  # largest eigenvalue = elongation axis
    normal = normal / np.linalg.norm(normal)

    # Place membrane center at the C-terminal (membrane-proximal) end.
    # Find the most C-terminal resolved residues as the anchor region.
    max_resnum = np.max(all_resnums)
    boundary_mask = all_resnums >= (max_resnum - 20)
    boundary_coords = all_coords[boundary_mask]
    boundary_center = np.mean(boundary_coords, axis=0)

    # Orient normal: bulk of protein (centroid) on the positive side
    if np.dot(centroid - boundary_center, normal) < 0:
        normal = -normal

    # Extrapolate membrane center 15A beyond the boundary (GPI linker + lipid)
    GPI_EXTRAPOLATION_A = 15.0
    membrane_center = boundary_center - normal * GPI_EXTRAPOLATION_A

    logger.info(
        "  Membrane plane estimated from GPI anchor (extrapolated %.0fA beyond C-terminus)",
        GPI_EXTRAPOLATION_A,
    )

    # Classify resolved residues as extracellular (only those N-terminal to GPI anchor)
    # Residues C-terminal to the GPI anchor are cleaved off or buried in membrane
    residue_topo = {}
    for res in chain:
        res_id = res.get_id()
        if res_id[0] != " " or "CA" not in res:
            continue
        resnum = res_id[1]
        if resnum <= gpi_residue:
            residue_topo[resnum] = "extracellular"
        else:
            # Residues after GPI anchor are not accessible (cleaved or membrane-embedded)
            residue_topo[resnum] = "cytoplasmic"

    return MembraneAnnotation(
        uniprot_id=target.uniprot_id,
        gene_name=target.gene_name,
        membrane_normal=normal,
        membrane_center=membrane_center,
        membrane_half_thickness=MEMBRANE_SLAB_HALF_THICKNESS,
        topology_type="gpi_anchored",
        tm_segments=[],
        residue_topology=residue_topo,
        source="gpi_anchor_estimated",
    )


# ---------------------------------------------------------------------------
# OPM database lookup
# ---------------------------------------------------------------------------

def _try_opm_lookup(pdb_id, chain_id):
    """
    Try to retrieve pre-calculated membrane orientation from OPM database.

    OPM provides PDB files oriented so that the membrane plane is the XY
    plane (z=0), with DUM atoms marking the membrane boundaries.

    Args:
        pdb_id: 4-character PDB ID.
        chain_id: Target chain.

    Returns:
        Tuple (normal_vector, center_point, half_thickness) or None.
    """
    cache_path = CACHE_DIR / "opm" / "{}.pdb".format(pdb_id.lower())

    if not cache_path.exists():
        url = "{}/pdb/{}".format(OPM_API, pdb_id.lower())
        try:
            resp = requests.get(url, timeout=15, allow_redirects=True)
            if resp.status_code != 200:
                logger.debug("  OPM lookup for %s: HTTP %d", pdb_id, resp.status_code)
                return None
            # Verify it looks like a PDB file
            content = resp.text
            if "ATOM" not in content:
                logger.debug("  OPM response for %s does not contain ATOM records", pdb_id)
                return None
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_path, "w") as f:
                f.write(content)
        except requests.RequestException as e:
            logger.debug("  OPM lookup failed for %s: %s", pdb_id, e)
            return None

    # Parse DUM atoms to find membrane boundaries
    return _parse_opm_membrane(str(cache_path))


def _parse_opm_membrane(opm_pdb_path):
    """
    Parse an OPM-oriented PDB file for membrane boundary markers.

    OPM convention:
      - Membrane plane is the XY plane (z=0)
      - DUM atoms (HETATM records with residue name DUM) mark the
        upper and lower leaflet boundaries
      - Membrane normal is the Z-axis [0, 0, 1]
      - Center is at the origin [0, 0, 0]
      - Half-thickness = max |z| of DUM atoms

    Returns:
        Tuple (normal, center, half_thickness) or None if no DUM atoms found.
    """
    dum_z_values = []
    with open(opm_pdb_path) as f:
        for line in f:
            if line.startswith("HETATM") and "DUM" in line:
                try:
                    z = float(line[46:54].strip())
                    dum_z_values.append(z)
                except (ValueError, IndexError):
                    continue

    if not dum_z_values:
        logger.debug("  No DUM atoms found in OPM file")
        return None

    half_thickness = max(abs(z) for z in dum_z_values)
    normal = np.array([0.0, 0.0, 1.0])
    center = np.array([0.0, 0.0, 0.0])

    logger.debug("  OPM membrane: half_thickness=%.1fA", half_thickness)
    return normal, center, half_thickness


# ---------------------------------------------------------------------------
# Membrane plane calculation from UniProt TM annotations
# ---------------------------------------------------------------------------

def _calculate_membrane_plane(chain, tm_segments, topology_type, features=None):
    """
    Compute the membrane plane from TM segment Cα coordinates using PCA.

    For single-pass proteins: the TM helix axis (largest PCA eigenvector)
    is approximately the membrane normal.

    For multi-pass proteins: compute the midpoint of each TM helix, then
    fit a plane through those midpoints (smallest PCA eigenvector = normal).

    If TM residues are not present in the structure (common for ectodomain-
    only crystal structures), falls back to estimating the membrane position
    from the terminus closest to the TM domain.

    Args:
        chain: BioPython Chain object.
        tm_segments: List of (start, end) residue number tuples.
        topology_type: "single_pass" or "multi_pass".
        features: UniProt features for topology context.

    Returns:
        Tuple (normal_vector, center_point) as numpy arrays.
    """
    if topology_type == "single_pass":
        ca_coords = _get_ca_coords(chain, tm_segments[0][0], tm_segments[0][1])
        if len(ca_coords) >= 3:
            return _single_pass_plane(chain, tm_segments[0])
        else:
            # TM residues not in structure — estimate from ectodomain boundary
            logger.info("  TM residues not resolved in structure, estimating membrane plane from ectodomain boundary")
            return _estimate_membrane_from_boundary(chain, tm_segments[0], features)
    else:
        # For multi-pass, check if any TM segments are resolved
        has_tm = False
        for start, end in tm_segments:
            ca = _get_ca_coords(chain, start, end)
            if len(ca) >= 3:
                has_tm = True
                break
        if has_tm:
            return _multi_pass_plane(chain, tm_segments)
        else:
            logger.info("  TM residues not resolved in structure, estimating from resolved boundaries")
            return _estimate_membrane_from_boundary(chain, tm_segments[0], features)


def _single_pass_plane(chain, tm_segment):
    """
    For a single TM helix, the helix axis approximates the membrane normal.

    Uses PCA on the TM Cα atoms: the eigenvector with the largest eigenvalue
    corresponds to the helix direction.
    """
    start, end = tm_segment
    ca_coords = _get_ca_coords(chain, start, end)

    if len(ca_coords) < 3:
        raise MembraneAnnotationError(
            "Too few Cα atoms ({}) in TM segment {}-{}".format(len(ca_coords), start, end)
        )

    centroid = np.mean(ca_coords, axis=0)
    centered = ca_coords - centroid

    # PCA: eigendecompose the covariance matrix
    cov = np.cov(centered.T)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)

    # eigh returns sorted ascending — largest eigenvalue is last
    # The helix axis corresponds to the largest eigenvalue (most variance)
    normal = eigenvectors[:, -1]

    # Normalize
    normal = normal / np.linalg.norm(normal)

    return normal, centroid


def _multi_pass_plane(chain, tm_segments):
    """
    For multi-pass proteins, fit a plane through the midpoints of all TM helices.

    The membrane plane passes through these midpoints. The normal is the
    eigenvector with the smallest eigenvalue (least variance direction).
    """
    midpoints = []
    all_ca_coords = []

    for start, end in tm_segments:
        ca_coords = _get_ca_coords(chain, start, end)
        if len(ca_coords) > 0:
            midpoints.append(np.mean(ca_coords, axis=0))
            all_ca_coords.extend(ca_coords.tolist())

    if len(midpoints) < 3:
        # If we have too few TM helices resolved, fall back to using
        # all TM Cα atoms directly
        all_ca = np.array(all_ca_coords)
        if len(all_ca) < 3:
            raise MembraneAnnotationError(
                "Too few TM Cα atoms ({}) to define membrane plane".format(len(all_ca))
            )
        centroid = np.mean(all_ca, axis=0)
        centered = all_ca - centroid
        cov = np.cov(centered.T)
        eigenvalues, eigenvectors = np.linalg.eigh(cov)
        # For a slab of helices, smallest eigenvalue = normal to the slab
        normal = eigenvectors[:, 0]
    else:
        midpoints = np.array(midpoints)
        centroid = np.mean(midpoints, axis=0)
        centered = midpoints - centroid
        cov = np.cov(centered.T)
        eigenvalues, eigenvectors = np.linalg.eigh(cov)
        # Smallest eigenvalue direction = normal to the plane of midpoints
        normal = eigenvectors[:, 0]

    normal = normal / np.linalg.norm(normal)
    return normal, centroid


def _estimate_membrane_from_boundary(chain, tm_segment, features):
    """
    Estimate the membrane plane when TM residues are not in the structure.

    This is common for ectodomain-only crystal structures. The approach:
      1. Collect all resolved Cα coordinates
      2. Find the residues closest to the TM boundary (the membrane-proximal end)
      3. Compute the ectodomain's principal axis via PCA (largest eigenvector)
      4. Place the membrane center at the boundary residues
      5. The principal axis direction serves as the membrane normal

    For type I single-pass (extracellular before TM), the membrane-proximal
    residues are those with the highest residue numbers near the TM start.

    Args:
        chain: BioPython Chain object.
        tm_segment: (start, end) of the TM segment from UniProt.
        features: UniProt features for topology context.

    Returns:
        Tuple (normal_vector, center_point) as numpy arrays.
    """
    tm_start, tm_end = tm_segment

    # Collect all resolved residue Cα coords with their residue numbers
    all_coords = []
    all_resnums = []
    for res in chain:
        res_id = res.get_id()
        if res_id[0] != " ":
            continue
        if "CA" in res:
            all_coords.append(res["CA"].get_vector().get_array())
            all_resnums.append(res_id[1])

    if len(all_coords) < 10:
        raise MembraneAnnotationError(
            "Too few resolved residues ({}) to estimate membrane plane".format(
                len(all_coords)
            )
        )

    all_coords = np.array(all_coords)
    all_resnums = np.array(all_resnums)

    # Determine which end is membrane-proximal by checking if ectodomain
    # is before or after the TM segment
    is_extracellular_before_tm = False
    for feat in (features or []):
        if feat["type"] == "Topological domain" and "Extracellular" in feat.get("description", ""):
            if feat.get("end") and int(feat["end"]) < tm_start:
                is_extracellular_before_tm = True
                break

    if is_extracellular_before_tm:
        # Type I: ectodomain is N-terminal, membrane is at the C-terminal end
        # Find the last ~10 resolved residues as the membrane boundary
        boundary_mask = all_resnums >= (tm_start - 30)
    else:
        # Type II or other: ectodomain is C-terminal
        boundary_mask = all_resnums <= (tm_end + 30)

    # If no boundary residues found, use the last/first 10%
    if not np.any(boundary_mask):
        n = max(10, len(all_resnums) // 10)
        if is_extracellular_before_tm:
            boundary_mask = np.zeros(len(all_resnums), dtype=bool)
            boundary_mask[-n:] = True
        else:
            boundary_mask = np.zeros(len(all_resnums), dtype=bool)
            boundary_mask[:n] = True

    boundary_coords = all_coords[boundary_mask]

    # PCA on all ectodomain Cα atoms — principal axis = ectodomain elongation
    centroid = np.mean(all_coords, axis=0)
    centered = all_coords - centroid
    cov = np.cov(centered.T)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)

    # Largest eigenvalue = principal axis of the ectodomain
    normal = eigenvectors[:, -1]
    normal = normal / np.linalg.norm(normal)

    # Orient: ectodomain centroid should be on the positive side
    boundary_center = np.mean(boundary_coords, axis=0)
    if np.dot(centroid - boundary_center, normal) < 0:
        normal = -normal

    # The membrane center is NOT at the boundary residues — it's BEYOND them.
    # The boundary residues are the last resolved ectodomain residues, but
    # the actual membrane is further out (unresolved linker + TM midpoint).
    # Extrapolate along the negative normal direction past the most
    # membrane-proximal resolved residue.
    #
    # Estimate: ~30A beyond the boundary centroid accounts for:
    #   - Unresolved ectodomain residues (~10-25 residues, ~15A in a coil)
    #   - Half the TM helix (~15A)
    MEMBRANE_EXTRAPOLATION_A = 30.0
    membrane_center = boundary_center - normal * MEMBRANE_EXTRAPOLATION_A

    logger.info(
        "  Estimated membrane center: extrapolated %.0fA beyond boundary residues (n=%d)",
        MEMBRANE_EXTRAPOLATION_A, np.sum(boundary_mask),
    )

    return normal, membrane_center


# ---------------------------------------------------------------------------
# Normal orientation
# ---------------------------------------------------------------------------

def _orient_normal(chain, normal, center, features):
    """
    Orient the membrane normal so that the extracellular side is positive.

    Uses UniProt "Topological domain" annotations to determine which side
    of the membrane is extracellular. If most extracellular-annotated
    residues have negative distance along the normal, flip it.

    Args:
        chain: BioPython Chain object.
        normal: Membrane normal vector.
        center: Membrane center point.
        features: UniProt features list.

    Returns:
        Oriented normal vector (may be flipped).
    """
    extra_distances = []
    for feat in features:
        if feat["type"] == "Topological domain" and "Extracellular" in feat.get("description", ""):
            start = feat.get("start")
            end = feat.get("end")
            if start and end:
                for res in chain:
                    res_id = res.get_id()
                    resnum = res_id[1]
                    if int(start) <= resnum <= int(end) and "CA" in res:
                        ca_pos = res["CA"].get_vector().get_array()
                        d = np.dot(ca_pos - center, normal)
                        extra_distances.append(d)

    if extra_distances:
        mean_d = np.mean(extra_distances)
        if mean_d < 0:
            logger.debug("  Flipping membrane normal (extracellular side was negative)")
            normal = -normal

    return normal


# ---------------------------------------------------------------------------
# Residue classification
# ---------------------------------------------------------------------------

def _classify_residues_from_plane(chain, normal, center, half_thickness, features, tm_segments):
    """
    Classify every residue as extracellular, transmembrane, or intracellular.

    Priority order:
      1. TM residues from UniProt Transmembrane annotations → "transmembrane"
      2. UniProt Topological domain annotations → "extracellular" or "intracellular"
      3. Geometric distance from membrane plane (fallback for unannotated residues)

    UniProt annotations take priority over geometry because AlphaFold
    predictions can place cytoplasmic tails on the wrong side of the
    membrane plane (e.g. disordered tails folding back).

    Args:
        chain: BioPython Chain object.
        normal: Membrane normal (extracellular side = positive).
        center: Membrane center point.
        half_thickness: Half-thickness of the membrane slab.
        features: UniProt features for topology classification.
        tm_segments: TM residue ranges for explicit TM classification.

    Returns:
        Dict mapping residue number -> "extracellular" | "transmembrane" | "intracellular".
    """
    # Build set of TM residue numbers
    tm_residue_set = set()
    for start, end in tm_segments:
        for r in range(start, end + 1):
            tm_residue_set.add(r)

    # Build UniProt topological domain map
    uniprot_topo = {}
    for feat in (features or []):
        if feat["type"] == "Topological domain":
            desc = feat.get("description", "")
            start = feat.get("start")
            end = feat.get("end")
            if start and end:
                if "Extracellular" in desc:
                    topo = "extracellular"
                elif "Cytoplasmic" in desc:
                    topo = "intracellular"
                else:
                    continue
                for r in range(int(start), int(end) + 1):
                    uniprot_topo[r] = topo

    topology = {}
    n_geometric = 0
    for res in chain:
        res_id = res.get_id()
        if res_id[0] != " ":
            continue
        resnum = res_id[1]

        if "CA" not in res:
            continue

        # Priority 1: explicit TM residues
        if resnum in tm_residue_set:
            topology[resnum] = "transmembrane"
            continue

        # Priority 2: UniProt topological domain annotation
        if resnum in uniprot_topo:
            topology[resnum] = uniprot_topo[resnum]
            continue

        # Priority 3: geometric fallback for unannotated residues
        ca_pos = res["CA"].get_vector().get_array()
        d = np.dot(ca_pos - center, normal)

        if abs(d) <= half_thickness:
            topology[resnum] = "transmembrane"
        elif d > half_thickness:
            topology[resnum] = "extracellular"
        else:
            topology[resnum] = "intracellular"
        n_geometric += 1

    if n_geometric > 0:
        logger.debug("  %d residues classified by geometry (no UniProt annotation)", n_geometric)

    return topology


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_chain(model, chain_id):
    """
    Get a chain from a BioPython model by ID.

    Falls back to the first chain if the requested chain_id is not found
    (common with AlphaFold predictions which use chain A).
    """
    if chain_id in model:
        return model[chain_id]
    # Fallback: try first chain
    chains = list(model.get_chains())
    if chains:
        logger.warning(
            "  Chain '%s' not found, using first chain '%s'",
            chain_id, chains[0].get_id(),
        )
        return chains[0]
    raise MembraneAnnotationError("No chains found in structure")


def _get_ca_coords(chain, start_res, end_res):
    """
    Extract Cα coordinates for a residue range.

    Args:
        chain: BioPython Chain object.
        start_res: Start residue number (inclusive).
        end_res: End residue number (inclusive).

    Returns:
        numpy array of shape (N, 3) with Cα coordinates.
    """
    coords = []
    for res in chain:
        res_id = res.get_id()
        if res_id[0] != " ":
            continue
        resnum = res_id[1]
        if start_res <= resnum <= end_res and "CA" in res:
            coords.append(res["CA"].get_vector().get_array())

    return np.array(coords) if coords else np.empty((0, 3))
