"""
Diagnostic: print per-residue sliding window mismatch density for
ERBB2 proximal rejected Patch 1.

Run: python -m epitope_pipeline.debug_sliding_window
"""

import math
import sys
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser

from epitope_pipeline import config
from epitope_pipeline.target_input import resolve_targets
from epitope_pipeline.structure import acquire_structure
from epitope_pipeline.membrane import annotate_membrane
from epitope_pipeline.spatial import filter_ectodomain
from epitope_pipeline.surface import analyze_surface
from epitope_pipeline.conservation import _align_sequences, _map_alignment_to_residues
from epitope_pipeline.config import VHH_FOOTPRINT_MIN_A2


def main():
    import tempfile

    # Resolve ERBB2
    targets = resolve_targets(["ERBB2"])
    target = targets[0]
    print(f"Target: {target.gene_name} ({target.uniprot_id})")
    print(f"  Sequence length: {len(target.sequence)} aa")
    print(f"  Cyno sequence length: {len(target.cyno_sequence)} aa")

    # Get structure (cached) - use temp dir for structures
    tmp_dir = Path(tempfile.mkdtemp())
    structure = acquire_structure(target, tmp_dir)
    membrane = annotate_membrane(target, structure)

    # Extract CA coords from PDB
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("erbb2", structure.pdb_path)
    model = struct[0]
    chain = model[structure.chain_id]

    ca_coords = {}
    for residue in chain.get_residues():
        if residue.id[0] != " ":
            continue
        resnum = residue.id[1]
        if "CA" in residue:
            ca_coords[resnum] = residue["CA"].get_vector().get_array()

    # Proximal mode: ECD residues <= 40A
    spatial = filter_ectodomain(target, structure, membrane,
                                max_distance=config.PROXIMAL_MAX_DISTANCE_A)
    print(f"\n  Proximal qualifying residues: {len(spatial.qualifying_residues)}")

    # Surface analysis
    surface = analyze_surface(target, structure, spatial)
    print(f"  Surface patches: {len(surface.patches)}")

    # Find the big rejected patch (Patch 1, 129 residues)
    big_patch = None
    for p in surface.patches:
        print(f"    Patch {p.patch_id}: {len(p.residue_numbers)} residues, "
              f"{p.total_sasa_a2:.0f} A²")
        if len(p.residue_numbers) > 50:
            big_patch = p

    if not big_patch:
        print("No large patch found!")
        sys.exit(1)

    print(f"\n{'='*80}")
    print(f"SLIDING WINDOW ANALYSIS — Patch {big_patch.patch_id}")
    print(f"  {len(big_patch.residue_numbers)} residues, {big_patch.total_sasa_a2:.0f} A²")
    print(f"{'='*80}")

    # Align human vs cyno
    aligned_human, aligned_cyno = _align_sequences(target.sequence, target.cyno_sequence)
    residue_conservation = _map_alignment_to_residues(aligned_human, aligned_cyno)

    # Identify mismatches in this patch
    patch_residues = sorted(big_patch.residue_numbers)
    mismatches = [r for r in patch_residues if not residue_conservation.get(r, False)]

    print(f"\n  Total mismatches in patch: {len(mismatches)}")
    print(f"  Mismatch positions: {mismatches}")

    # Show human vs cyno at each mismatch
    print(f"\n  Mismatch details:")
    for r in mismatches:
        human_aa = target.sequence[r - 1] if r <= len(target.sequence) else "?"
        # Get cyno AA from alignment
        cyno_aa = "?"
        human_pos = 0
        for ah, ac in zip(aligned_human, aligned_cyno):
            if ah != "-":
                human_pos += 1
            if human_pos == r and ah != "-":
                cyno_aa = ac
                break
        dist = spatial.residue_distances.get(r, float("nan"))
        print(f"    Res {r}: {human_aa} -> {cyno_aa}  (dist from membrane: {dist:.1f}A)")

    # Run sliding window
    footprint_radius = math.sqrt(VHH_FOOTPRINT_MIN_A2 / math.pi)
    print(f"\n  Footprint radius: {footprint_radius:.1f}A (~600A² circle)")

    res_with_coords = [r for r in patch_residues if r in ca_coords]
    coords = np.array([ca_coords[r] for r in res_with_coords])
    is_mismatch = np.array([
        not residue_conservation.get(r, False) for r in res_with_coords
    ])
    dists = cdist(coords, coords)

    print(f"\n  Per-residue sliding window (residues with ≥1 mismatch in window):")
    print(f"  {'Res':>5}  {'AA':>2}  {'Conserved':>9}  {'Window_Mismatches':>17}  "
          f"{'Neighbors':>9}  {'DistFromMembrane':>16}  {'MismatchResNums'}")
    print(f"  {'-'*5}  {'-'*2}  {'-'*9}  {'-'*17}  {'-'*9}  {'-'*16}  {'-'*20}")

    results = []
    for i, resnum in enumerate(res_with_coords):
        neighbors_mask = dists[i] <= footprint_radius
        local_bad = int(np.sum(is_mismatch[neighbors_mask]))
        n_neighbors = int(np.sum(neighbors_mask))

        # Which mismatch residues are in this window?
        mismatch_in_window = [
            res_with_coords[j] for j in range(len(res_with_coords))
            if neighbors_mask[j] and is_mismatch[j]
        ]

        results.append((resnum, local_bad, n_neighbors, mismatch_in_window))

    # Print all residues that have at least 1 mismatch in their window
    for resnum, local_bad, n_neighbors, mismatch_in_window in results:
        if local_bad == 0:
            continue
        human_aa = target.sequence[resnum - 1] if resnum <= len(target.sequence) else "?"
        conserved = residue_conservation.get(resnum, False)
        dist = spatial.residue_distances.get(resnum, float("nan"))
        mm_str = ",".join(str(r) for r in mismatch_in_window)
        flag = " *** OVER LIMIT" if local_bad > config.MAX_CYNO_MISMATCHES_PER_600A2 else ""
        print(f"  {resnum:>5}  {human_aa:>2}  {'yes' if conserved else 'NO':>9}  "
              f"{local_bad:>17}  {n_neighbors:>9}  {dist:>16.1f}  {mm_str}{flag}")

    # Summary
    worst_pos = max(results, key=lambda x: x[1])
    print(f"\n  WORST: Residue {worst_pos[0]} with {worst_pos[1]} mismatches in window")
    print(f"  Threshold: max {config.MAX_CYNO_MISMATCHES_PER_600A2} per ~600A² window")
    if worst_pos[1] > config.MAX_CYNO_MISMATCHES_PER_600A2:
        print(f"  VERDICT: REJECTED (worst window has {worst_pos[1]} > {config.MAX_CYNO_MISMATCHES_PER_600A2})")
    else:
        print(f"  VERDICT: PASSED")

    # Also print where the mismatches are NOT in the patch (i.e. residues
    # in the 400-700 range that didn't make it into the patch at all)
    print(f"\n{'='*80}")
    print(f"RESIDUES 400-700: Why aren't they all in the patch?")
    print(f"{'='*80}")
    patch_set = set(patch_residues)
    qualifying_set = set(spatial.qualifying_residues)
    exposed_set = set(surface.exposed_residues) if hasattr(surface, 'exposed_residues') else None

    print(f"\n  {'Res':>5}  {'ECD?':>5}  {'<=40A?':>6}  {'Exposed?':>8}  "
          f"{'InPatch?':>8}  {'Dist':>8}  {'relSASA':>8}")
    print(f"  {'-'*5}  {'-'*5}  {'-'*6}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")

    # Get per-residue SASA
    residue_sasa = surface.residue_sasa if hasattr(surface, 'residue_sasa') else {}

    for r in range(400, 701):
        is_ecd = r in set(range(target.ectodomain_start, target.ectodomain_end + 1)) if hasattr(target, 'ectodomain_start') else None
        is_qualifying = r in qualifying_set
        is_in_patch = r in patch_set
        dist = spatial.residue_distances.get(r, float("nan"))
        sasa = residue_sasa.get(r, 0.0)
        # Only print every 5th residue unless it's interesting
        is_interesting = is_qualifying or is_in_patch or r % 10 == 0
        if not is_interesting:
            continue

        # Check if extracellular from membrane annotation
        ecd_str = "?"
        if hasattr(membrane, 'residue_topology'):
            topo = membrane.residue_topology.get(r, "?")
            ecd_str = "ECD" if topo == "Extracellular" else topo[:3]

        qual_str = "yes" if is_qualifying else "no"
        in_patch_str = "yes" if is_in_patch else "no"

        # Check if exposed
        exposed_str = "?"
        if exposed_set is not None:
            exposed_str = "yes" if r in exposed_set else "no"
        elif sasa > 0:
            from epitope_pipeline.config import MAX_SASA, AA_1TO3
            aa = target.sequence[r - 1] if r <= len(target.sequence) else "X"
            max_s = MAX_SASA.get(AA_1TO3.get(aa, "ALA"), 180.0)
            rel = sasa / max_s if max_s > 0 else 0
            exposed_str = f"{rel:.2f}"

        print(f"  {r:>5}  {ecd_str:>5}  {qual_str:>6}  {exposed_str:>8}  "
              f"{in_patch_str:>8}  {dist:>8.1f}  {sasa:>8.1f}")


if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.WARNING)
    main()
