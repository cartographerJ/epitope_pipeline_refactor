"""
PDB-file helpers: chain access and Cα-coordinate extraction.

These were previously in utils.py — moved here because they're I/O concerns
(BioPython parsing) rather than general-purpose utilities.
"""

from Bio.PDB import PDBParser


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
