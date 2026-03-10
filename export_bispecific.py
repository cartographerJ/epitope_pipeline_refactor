"""
Bispecific export — CSV, XLSX, zone-specific PDBs, and dual PML scripts.
"""

import csv
import logging
from pathlib import Path
from typing import Dict, List

from epitope_pipeline.config import PALETTE, DUAL_PML_GAP_A
from epitope_pipeline.export import (
    export_annotated_pdb,
    generate_membrane_cgo_pml,
    generate_shared_membrane_cgo_pml,
)

logger = logging.getLogger("epitope_pipeline")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def export_bispecific_all(run_dir, pair_results, zone_results,
                          structures, membranes):
    """
    Generate all bispecific-specific output files.

    Args:
        run_dir: Path to the run directory.
        pair_results: List of BispecificPairResult.
        zone_results: Dict {(uid, zone): TargetZoneResult}.
        structures: Dict {uid: StructureResult}.
        membranes: Dict {uid: MembraneAnnotation}.
    """
    run_dir = Path(run_dir)
    pymol_dir = run_dir / "pymol"
    pymol_dir.mkdir(exist_ok=True)
    zone_sessions_dir = pymol_dir / "zone_sessions"
    zone_sessions_dir.mkdir(exist_ok=True)

    # 1. CSV
    export_bispecific_csv(run_dir, pair_results)

    # 2. XLSX
    export_bispecific_xlsx(run_dir, pair_results)

    # 3. Zone-specific annotated PDBs + PMLs → pymol/zone_sessions/
    exported_zones = set()
    for key, zr in zone_results.items():
        uid, zone = key
        if (uid, zone) in exported_zones:
            continue
        export_zone_annotated_pdb(run_dir, zr)
        exported_zones.add((uid, zone))

    # 4. Dual PML scripts (one per valid orientation) → pymol/
    for pr in pair_results:
        for orientation in (pr.orientation_ab, pr.orientation_ba):
            if not orientation.is_valid:
                continue
            export_bispecific_pml(
                run_dir, pr, orientation, structures, membranes,
            )

    logger.info("  Bispecific outputs written to: %s", run_dir)


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------

def export_bispecific_csv(run_dir, pair_results):
    """Write bispecific_pairs.csv."""
    csv_path = Path(run_dir) / "bispecific_pairs.csv"
    rows = []

    for pr in pair_results:
        for orientation in (pr.orientation_ab, pr.orientation_ba):
            distal_gene = orientation.distal_zone.target.gene_name if orientation.distal_zone else "—"
            proximal_gene = orientation.proximal_zone.target.gene_name if orientation.proximal_zone else "—"

            distal_best = 0.0
            distal_area = 0.0
            distal_patch_id = "—"
            if orientation.distal_zone and orientation.distal_zone.scores:
                s = orientation.distal_zone.scores[0]
                distal_best = s.composite_score
                distal_area = s.patch_area_a2
                distal_patch_id = s.patch_id

            proximal_best = 0.0
            proximal_area = 0.0
            proximal_patch_id = "—"
            if orientation.proximal_zone and orientation.proximal_zone.scores:
                s = orientation.proximal_zone.scores[0]
                proximal_best = s.composite_score
                proximal_area = s.patch_area_a2
                proximal_patch_id = s.patch_id

            rows.append({
                "pair": "{} x {}".format(pr.target_a.gene_name, pr.target_b.gene_name),
                "target_a": pr.target_a.gene_name,
                "target_b": pr.target_b.gene_name,
                "orientation": orientation.label,
                "distal_target": distal_gene,
                "proximal_target": proximal_gene,
                "distal_best_patch": distal_patch_id,
                "distal_best_score": round(distal_best, 4),
                "distal_best_area_A2": round(distal_area, 1),
                "proximal_best_patch": proximal_patch_id,
                "proximal_best_score": round(proximal_best, 4),
                "proximal_best_area_A2": round(proximal_area, 1),
                "orientation_score": round(orientation.orientation_score, 4),
                "is_valid": orientation.is_valid,
                "dual_valid": pr.both_valid,
                "flexibility_bonus": pr.both_valid,
                "final_pair_score": round(pr.final_pair_score, 4),
            })

    if rows:
        fieldnames = list(rows[0].keys())
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        logger.info("  Wrote %s (%d rows)", csv_path.name, len(rows))


# ---------------------------------------------------------------------------
# XLSX
# ---------------------------------------------------------------------------

def export_bispecific_xlsx(run_dir, pair_results):
    """Write bispecific_pairs.xlsx with Cartography styling."""
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
    from openpyxl.utils import get_column_letter

    xlsx_path = Path(run_dir) / "bispecific_pairs.xlsx"
    wb = Workbook()

    # Styles
    header_fill = PatternFill(start_color="28154C", end_color="28154C", fill_type="solid")
    header_font = Font(name="Consolas", bold=True, color="FFFFFF", size=10)
    data_font = Font(name="Consolas", size=10)
    score_font = Font(name="Consolas", size=10, bold=True)
    valid_fill = PatternFill(start_color="6BC291", end_color="6BC291", fill_type="solid")
    invalid_fill = PatternFill(start_color="D3D3D3", end_color="D3D3D3", fill_type="solid")
    dual_fill = PatternFill(start_color="18B5CB", end_color="18B5CB", fill_type="solid")
    thin_border = Border(
        left=Side(style="thin"), right=Side(style="thin"),
        top=Side(style="thin"), bottom=Side(style="thin"),
    )

    # === Sheet 1: Pair Summary ===
    ws = wb.active
    ws.title = "Pair Summary"

    headers = [
        "Pair", "Target A", "Target B", "Final Pair Score",
        "Dual-Valid", "Best Orientation",
        "Distal Target", "Distal Best Score",
        "Proximal Target", "Proximal Best Score",
    ]
    for col, h in enumerate(headers, 1):
        cell = ws.cell(row=1, column=col, value=h)
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center")
        cell.border = thin_border

    for row_idx, pr in enumerate(pair_results, 2):
        best = pr.best_orientation
        distal_gene = best.distal_zone.target.gene_name if best.distal_zone else "—"
        proximal_gene = best.proximal_zone.target.gene_name if best.proximal_zone else "—"
        distal_score = best.distal_zone.best_score if best.distal_zone else 0.0
        proximal_score = best.proximal_zone.best_score if best.proximal_zone else 0.0

        values = [
            "{} x {}".format(pr.target_a.gene_name, pr.target_b.gene_name),
            pr.target_a.gene_name,
            pr.target_b.gene_name,
            round(pr.final_pair_score, 4),
            "Yes" if pr.both_valid else "No",
            best.label,
            distal_gene,
            round(distal_score, 4),
            proximal_gene,
            round(proximal_score, 4),
        ]

        if pr.both_valid:
            fill = dual_fill
        elif pr.final_pair_score > 0:
            fill = valid_fill
        else:
            fill = invalid_fill

        for col, val in enumerate(values, 1):
            cell = ws.cell(row=row_idx, column=col, value=val)
            cell.font = score_font if col == 4 else data_font
            cell.fill = fill
            cell.border = thin_border
            cell.alignment = Alignment(horizontal="center")

    for col in range(1, len(headers) + 1):
        ws.column_dimensions[get_column_letter(col)].width = 20

    # === Sheet 2: All Orientations ===
    ws2 = wb.create_sheet("Orientation Detail")

    orient_headers = [
        "Pair", "Orientation", "Valid",
        "Distal Target", "Distal Score", "Distal Patches", "Distal Area (A2)",
        "Proximal Target", "Proximal Score", "Proximal Patches", "Proximal Area (A2)",
        "Orientation Score",
    ]
    for col, h in enumerate(orient_headers, 1):
        cell = ws2.cell(row=1, column=col, value=h)
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center")
        cell.border = thin_border

    row_idx = 2
    for pr in pair_results:
        for orientation in (pr.orientation_ab, pr.orientation_ba):
            dz = orientation.distal_zone
            pz = orientation.proximal_zone

            values = [
                "{} x {}".format(pr.target_a.gene_name, pr.target_b.gene_name),
                orientation.label,
                "Yes" if orientation.is_valid else "No",
                dz.target.gene_name if dz else "—",
                round(dz.best_score, 4) if dz else 0.0,
                len(dz.scores) if dz else 0,
                round(sum(s.patch_area_a2 for s in dz.scores), 1) if dz and dz.scores else 0.0,
                pz.target.gene_name if pz else "—",
                round(pz.best_score, 4) if pz else 0.0,
                len(pz.scores) if pz else 0,
                round(sum(s.patch_area_a2 for s in pz.scores), 1) if pz and pz.scores else 0.0,
                round(orientation.orientation_score, 4),
            ]

            fill = valid_fill if orientation.is_valid else invalid_fill
            for col, val in enumerate(values, 1):
                cell = ws2.cell(row=row_idx, column=col, value=val)
                cell.font = data_font
                cell.fill = fill
                cell.border = thin_border
                cell.alignment = Alignment(horizontal="center")
            row_idx += 1

    for col in range(1, len(orient_headers) + 1):
        ws2.column_dimensions[get_column_letter(col)].width = 18

    wb.save(str(xlsx_path))
    logger.info("  Wrote %s", xlsx_path.name)


# ---------------------------------------------------------------------------
# Zone-specific annotated PDBs
# ---------------------------------------------------------------------------

def export_zone_annotated_pdb(run_dir, zone_result):
    """
    Write {gene}_{zone}_epitope.pdb with zone-specific B-factor tiers.

    Reuses the existing export_annotated_pdb() from export.py but writes
    to a zone-specific filename.
    """
    target = zone_result.target
    gene = target.gene_name.lower()
    zone = zone_result.zone

    # Temporarily rename gene_name to get zone-specific filenames
    # from export_annotated_pdb (which uses target.gene_name.lower())
    original_name = target.gene_name
    target.gene_name = "{}_{}".format(original_name, zone)

    try:
        export_annotated_pdb(
            run_dir=str(run_dir),
            target=target,
            structure=zone_result.structure,
            membrane=zone_result.membrane,
            spatial_filter=zone_result.spatial_filter,
            surface_analysis=zone_result.surface_analysis,
            conservation_result=zone_result.conservation_result,
            scores=zone_result.scores,
            pdb_subdir="pymol/zone_sessions",
        )
    finally:
        # Restore original name
        target.gene_name = original_name


# ---------------------------------------------------------------------------
# Dual PML script
# ---------------------------------------------------------------------------

def export_bispecific_pml(run_dir, pair_result, orientation,
                          structures, membranes):
    """
    Write a dual-target PyMOL script with both structures side-by-side.
    """
    run_dir = Path(run_dir)
    pml_dir = run_dir / "pymol"
    pml_dir.mkdir(parents=True, exist_ok=True)

    dz = orientation.distal_zone
    pz = orientation.proximal_zone
    distal_gene = dz.target.gene_name.lower()
    proximal_gene = pz.target.gene_name.lower()

    distal_pdb = "{}_distal_epitope".format(distal_gene)
    proximal_pdb = "{}_proximal_epitope".format(proximal_gene)

    pml_path = pml_dir / "{}.pml".format(orientation.label.lower())

    # Build residue sets for each target's patches
    distal_patches = {}
    for s in dz.scores:
        distal_patches[s.patch_id] = sorted(s.patch.residue_numbers)

    proximal_patches = {}
    for s in pz.scores:
        proximal_patches[s.patch_id] = sorted(s.patch.residue_numbers)

    # Helper: format residue numbers for PyMOL resi selection
    def resi_str(resnums):
        if not resnums:
            return "none"
        nums = sorted(resnums)
        ranges = []
        start = nums[0]
        end = nums[0]
        for n in nums[1:]:
            if n == end + 1:
                end = n
            else:
                ranges.append("{}-{}".format(start, end) if end > start else str(start))
                start = end = n
        ranges.append("{}-{}".format(start, end) if end > start else str(start))
        return "+".join(ranges)

    # Determine chains
    distal_chain = dz.structure.chain_id
    proximal_chain = pz.structure.chain_id

    # Cytoplasmic residues to hide
    distal_cyto = set()
    proximal_cyto = set()
    if dz.membrane:
        for resnum, topo in dz.membrane.residue_topology.items():
            if topo == "intracellular":
                distal_cyto.add(resnum)
    if pz.membrane:
        for resnum, topo in pz.membrane.residue_topology.items():
            if topo == "intracellular":
                proximal_cyto.add(resnum)

    gap = DUAL_PML_GAP_A

    # --- Compute membrane-aligned view basis EARLY ---
    # This determines the translation direction so structures appear
    # purely side-by-side on screen (no vertical offset).
    import numpy as np
    normal_d = np.array(dz.membrane.membrane_normal)
    normal_p = np.array(pz.membrane.membrane_normal)
    avg_normal = normal_d + normal_p
    n_len = np.linalg.norm(avg_normal)
    avg_normal = avg_normal / n_len if n_len > 0 else np.array([0., 1., 0.])

    # Orthonormal basis: screen_y = membrane normal, screen_x = side-by-side
    up = avg_normal
    right = np.array([1.0, 0.0, 0.0])
    right = right - np.dot(right, up) * up  # Gram-Schmidt
    if np.linalg.norm(right) < 0.01:
        right = np.array([0.0, 0.0, 1.0])
        right = right - np.dot(right, up) * up
    right = right / np.linalg.norm(right)
    fwd = np.cross(right, up)
    fwd = fwd / np.linalg.norm(fwd)

    # Translate along the right vector so structures separate purely
    # horizontally on screen (no vertical component from tilted normal).
    # Also align membrane planes vertically: shift proximal along normal
    # so both TM regions sit at the same height.
    center_d = np.array(dz.membrane.membrane_center)
    center_p_raw = np.array(pz.membrane.membrane_center)
    vertical_offset = (np.dot(center_d, up) - np.dot(center_p_raw, up)) * up
    translate_vec = right * gap + vertical_offset

    center_p = center_p_raw + translate_vec
    mem_center = (center_d + center_p) / 2.0

    lines = [
        "# ================================================================",
        "# BISPECIFIC VISUALIZATION",
        "# {} (distal) x {} (proximal)".format(
            dz.target.gene_name, pz.target.gene_name),
        "# Pair score: {:.3f}{}".format(
            pair_result.final_pair_score,
            " [DUAL-VALID, {:.0f}% bonus]".format(
                (pair_result.final_pair_score / orientation.orientation_score - 1) * 100
            ) if pair_result.both_valid and orientation.orientation_score > 0 else "",
        ),
        "# ================================================================",
        "",
        "# --- Cartography palette ---",
        "set_color carto_gray,       [0.827, 0.827, 0.827]",
        "set_color carto_palegreen,  [0.878, 0.953, 0.859]",
        "set_color carto_green,      [0.420, 0.761, 0.569]",
        "set_color carto_teal,       [0.094, 0.710, 0.796]",
        "set_color carto_blue,       [0.180, 0.584, 0.824]",
        "set_color carto_purple,     [0.157, 0.082, 0.298]",
        "",
        "# --- Load both structures ---",
        "load zone_sessions/{}.pdb, {}".format(distal_pdb, distal_pdb),
        "load zone_sessions/{}.pdb, {}".format(proximal_pdb, proximal_pdb),
        "",
        "# --- Translate proximal target for side-by-side view ---",
        "translate [{:.2f}, {:.2f}, {:.2f}], {}".format(
            translate_vec[0], translate_vec[1], translate_vec[2], proximal_pdb),
        "",
        "# --- Base display ---",
        "hide everything",
        "show cartoon, {}".format(distal_pdb),
        "show cartoon, {}".format(proximal_pdb),
        "color white, {}".format(distal_pdb),
        "color white, {}".format(proximal_pdb),
        "set cartoon_transparency, 0.3, {}".format(distal_pdb),
        "set cartoon_transparency, 0.3, {}".format(proximal_pdb),
        "",
    ]

    # Hide cytoplasmic for both
    if distal_cyto:
        lines.append("# --- Hide cytoplasmic: {} ---".format(dz.target.gene_name))
        lines.append("select _distal_cyto, {} and chain {} and resi {}".format(
            distal_pdb, distal_chain, resi_str(distal_cyto)))
        lines.append("hide cartoon, _distal_cyto")
        lines.append("disable _distal_cyto")
        lines.append("")

    if proximal_cyto:
        lines.append("# --- Hide cytoplasmic: {} ---".format(pz.target.gene_name))
        lines.append("select _proximal_cyto, {} and chain {} and resi {}".format(
            proximal_pdb, proximal_chain, resi_str(proximal_cyto)))
        lines.append("hide cartoon, _proximal_cyto")
        lines.append("disable _proximal_cyto")
        lines.append("")

    # Distal epitope patches
    lines.append("# =============================================")
    lines.append("# DISTAL TARGET: {} (>= {}A)".format(
        dz.target.gene_name,
        dz.spatial_filter.min_distance_threshold if dz.spatial_filter else "?"))
    lines.append("# =============================================")
    lines.append("")

    for pid, resnums in distal_patches.items():
        sel_name = "Distal_Patch_{}".format(pid)
        lines.append("select {}, {} and chain {} and resi {}".format(
            sel_name, distal_pdb, distal_chain, resi_str(resnums)))
        lines.append("color carto_green, {}".format(sel_name))
        lines.append("show surface, {}".format(sel_name))
        lines.append("set transparency, 0.3, {}".format(sel_name))
        lines.append("set cartoon_transparency, 0.0, {}".format(sel_name))
        lines.append("")

    # Proximal epitope patches
    lines.append("# =============================================")
    lines.append("# PROXIMAL TARGET: {} (<= {}A)".format(
        pz.target.gene_name,
        pz.spatial_filter.max_distance_threshold if pz.spatial_filter else "?"))
    lines.append("# =============================================")
    lines.append("")

    for pid, resnums in proximal_patches.items():
        sel_name = "Proximal_Patch_{}".format(pid)
        lines.append("select {}, {} and chain {} and resi {}".format(
            sel_name, proximal_pdb, proximal_chain, resi_str(resnums)))
        lines.append("color carto_green, {}".format(sel_name))
        lines.append("show surface, {}".format(sel_name))
        lines.append("set transparency, 0.3, {}".format(sel_name))
        lines.append("set cartoon_transparency, 0.0, {}".format(sel_name))
        lines.append("")

    # Shared membrane bilayer CGO spanning both structures
    pdb_dir = run_dir / "pymol" / "zone_sessions"
    if dz.membrane and pz.membrane:
        distal_pdb_path = str(pdb_dir / "{}.pdb".format(distal_pdb))
        proximal_pdb_path = str(pdb_dir / "{}.pdb".format(proximal_pdb))
        lines.extend(generate_shared_membrane_cgo_pml(
            membrane_a=dz.membrane, pdb_path_a=distal_pdb_path, chain_a=distal_chain,
            membrane_b=pz.membrane, pdb_path_b=proximal_pdb_path, chain_b=proximal_chain,
            translate_b=translate_vec,
        ))
    elif dz.membrane:
        d_pdb_path = str(pdb_dir / "{}.pdb".format(distal_pdb))
        lines.extend(generate_membrane_cgo_pml(
            dz.membrane, d_pdb_path, distal_chain,
        ))
    elif pz.membrane:
        p_pdb_path = str(pdb_dir / "{}.pdb".format(proximal_pdb))
        lines.extend(generate_membrane_cgo_pml(
            pz.membrane, p_pdb_path, proximal_chain,
            translate_vec=translate_vec,
        ))

    # Labels
    lines.extend([
        "# --- Labels ---",
        "set label_size, 16",
        "set label_color, black",
        "",
    ])

    # --- Membrane-aligned view (basis computed above) ---
    lines.extend([
        "# --- Membrane-aligned view ---",
        "deselect",
        "set_view (\\",
        "    {:.6f}, {:.6f}, {:.6f},\\".format(right[0], right[1], right[2]),
        "    {:.6f}, {:.6f}, {:.6f},\\".format(up[0], up[1], up[2]),
        "    {:.6f}, {:.6f}, {:.6f},\\".format(fwd[0], fwd[1], fwd[2]),
        "    0.000000, 0.000000, -400.000000,\\",
        "    {:.2f}, {:.2f}, {:.2f},\\".format(
            mem_center[0], mem_center[1], mem_center[2]),
        "    100.000000, 900.000000, -20.000000)",
        "",
        "# Ensure ectodomain faces UP (flip if center of mass is below membrane)",
        "python",
        "from pymol import cmd",
        "v = cmd.get_view()",
        "com = cmd.centerofmass('visible')",
        "mcx, mcy, mcz = {:.2f}, {:.2f}, {:.2f}".format(
            mem_center[0], mem_center[1], mem_center[2]),
        "dx = com[0] - mcx",
        "dy = com[1] - mcy",
        "dz2 = com[2] - mcz",
        "sy_com = v[3]*dx + v[4]*dy + v[5]*dz2",
        "if sy_com < 0:",
        "    cmd.turn('x', 180)",
        "python end",
        "",
        "zoom",
        "",
    ])

    with open(pml_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    logger.info("  Wrote dual PML: %s", pml_path.name)
