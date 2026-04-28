"""
Bispecific Dual-Targeting — evaluate target pairs for complementary epitopes.

For bispecific antibody design, one arm binds a "distal" epitope (>=60A from
membrane) and the other binds a "proximal" epitope (<=40A from membrane).
Both orientations are evaluated for each pair:
  - Target A distal + Target B proximal
  - Target A proximal + Target B distal

Usage:
    from epitope_pipeline.bispecific import run_bispecific
    results = run_bispecific([("ERBB2", "NECTIN4")])

    python -m epitope_pipeline.bispecific ERBB2:NECTIN4
"""

import logging
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from epitope_pipeline import config
from epitope_pipeline.utils import setup_logging
from epitope_pipeline.io.pdb import extract_ca_coords
from epitope_pipeline.io.targets import resolve_targets, TargetResolutionError
from epitope_pipeline.io.structure import acquire_structure, StructureAcquisitionError
from epitope_pipeline.io.membrane import annotate_membrane, MembraneAnnotationError
from epitope_pipeline.spatial import filter_ectodomain
from epitope_pipeline.surface import analyze_surface, cluster_ectodomain_patches
from epitope_pipeline.conservation import (
    analyze_conservation, ConservationError, ConservationResult,
)
from epitope_pipeline.specificity import filter_specificity
from epitope_pipeline.scoring import score_epitopes, compute_target_epitope_metric
from epitope_pipeline.visualize import plot_epitope_map, plot_blast_offtargets


logger = logging.getLogger("epitope_pipeline")


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class TargetZoneResult:
    """Results for one target in one spatial mode (distal or proximal)."""
    target: object                      # TargetInfo
    zone: str                           # "distal" or "proximal"
    structure: object                   # StructureResult
    membrane: object                    # MembraneAnnotation
    spatial_filter: object              # SpatialFilter
    surface_analysis: object            # SurfaceAnalysis
    conservation_result: object         # ConservationResult or None
    specificity_result: object          # SpecificityResult or None
    scores: list = field(default_factory=list)  # [EpitopeScore, ...]
    best_score: float = 0.0            # Best composite score (0.0 if no patches)


@dataclass
class OrientationResult:
    """One orientation of a bispecific pair."""
    distal_zone: TargetZoneResult       # The target playing "distal"
    proximal_zone: TargetZoneResult     # The target playing "proximal"
    label: str                          # e.g. "ERBB2_distal__MSLN_proximal"
    orientation_score: float            # geomean(best_distal, best_proximal)
    is_valid: bool                      # Both sides have >=1 qualifying patch


@dataclass
class BispecificPairResult:
    """Complete evaluation of one target pair."""
    target_a: object                    # TargetInfo
    target_b: object                    # TargetInfo
    orientation_ab: OrientationResult   # A=distal, B=proximal
    orientation_ba: OrientationResult   # A=proximal, B=distal
    both_valid: bool                    # Both orientations work
    best_orientation: OrientationResult
    final_pair_score: float             # best_score * flexibility_bonus if both_valid


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_bispecific(
    pairs,
    run_name=None,
    distal_min_distance_a=None,
    proximal_max_distance_a=None,
    cyno_max_mismatches=None,
    cyno_mismatch_percent=None,
    nonspecific_percent=None,
    force_experimental=False,
    verbose=True,
):
    """
    Evaluate bispecific antibody epitope pairs.

    Args:
        pairs: List of (identifier_a, identifier_b) tuples.
        run_name: Optional custom run directory name.
        distal_min_distance_a: Override distal threshold (default 60A).
        proximal_max_distance_a: Override proximal threshold (default 40A).
        cyno_max_mismatches: Override max cyno mismatches per 600A² (default 2) [DEPRECATED].
        cyno_mismatch_percent: Override cyno mismatch percent threshold (default 15.0).
        nonspecific_percent: Override nonspecific percent threshold (default 15.0).
        force_experimental: Use experimental PDB instead of AlphaFold.
        verbose: Whether to log to console.

    Returns:
        Dict with run_dir, pair_results, zone_results.
    """
    distal_dist = distal_min_distance_a or config.DISTAL_MIN_DISTANCE_A
    proximal_dist = proximal_max_distance_a or config.PROXIMAL_MAX_DISTANCE_A

    # Apply threshold overrides (new whole-patch percent thresholds)
    if cyno_mismatch_percent is not None:
        config.MAX_CYNO_MISMATCH_PERCENT = cyno_mismatch_percent
    if nonspecific_percent is not None:
        config.MAX_NONSPECIFIC_PERCENT = nonspecific_percent

    # Apply legacy threshold overrides (deprecated)
    if cyno_max_mismatches is not None:
        config.MAX_CYNO_MISMATCHES_PER_600A2 = cyno_max_mismatches

    # --- Setup logging ---
    setup_logging(verbose)
    logger.info("=" * 70)
    logger.info("BISPECIFIC EPITOPE PIPELINE v0.1.0")
    logger.info("=" * 70)
    logger.info("Pairs: %s", ", ".join("{} x {}".format(a, b) for a, b in pairs))
    logger.info("Distal threshold: >= %.0fA", distal_dist)
    logger.info("Proximal threshold: <= %.0fA", proximal_dist)

    # --- Create run directory ---
    if run_name is None:
        pair_slug = "_".join(
            "{}_{}".format(a.lower()[:6], b.lower()[:6]) for a, b in pairs[:2]
        )
        if len(pairs) > 2:
            pair_slug += "_etc"
        run_name = "{}_bispecific_{}".format(
            datetime.now().strftime("%y%m%d_%H%M"), pair_slug
        )

    run_dir = config.RUNS_DIR / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    structures_dir = run_dir / ".structures"
    structures_dir.mkdir(exist_ok=True)
    figures_dir = run_dir / "Figures"
    figures_dir.mkdir(exist_ok=True)
    zone_details_dir = figures_dir / "zone_details"
    zone_details_dir.mkdir(exist_ok=True)
    pymol_dir = run_dir / "Structures"
    pymol_dir.mkdir(exist_ok=True)

    # Add file handler
    supp_dir = run_dir / "Supplementary Files"
    supp_dir.mkdir(exist_ok=True)
    logs_dir = supp_dir / "Logs"
    logs_dir.mkdir(exist_ok=True)
    log_path = logs_dir / "log.txt"
    fh = logging.FileHandler(str(log_path))
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(fh)
    logger.info("Run directory: %s", run_dir)

    # =====================================================================
    # Step 1: Deduplicate and resolve all unique targets
    # =====================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 1: Resolving targets")
    logger.info("=" * 70)

    unique_ids = list(dict.fromkeys(
        id_ for pair in pairs for id_ in pair
    ))  # Preserves order, deduplicates

    try:
        all_targets = resolve_targets(unique_ids)
    except TargetResolutionError as e:
        logger.error("Target resolution failed: %s", e)
        logger.removeHandler(fh)
        fh.close()
        return {"run_dir": str(run_dir), "pair_results": [], "zone_results": {}}

    # Build lookup by gene name, uniprot_id, and original input aliases
    target_lookup = {}
    for t in all_targets:
        target_lookup[t.uniprot_id] = t
        target_lookup[t.gene_name.upper()] = t
    # Also map original input IDs (e.g., "HER2" alias → resolved ERBB2 target)
    for orig_id, resolved in zip(unique_ids, all_targets):
        target_lookup[orig_id.upper()] = resolved

    logger.info("Resolved %d unique targets", len(all_targets))

    # =====================================================================
    # Steps 2-3: Acquire structures and annotate membranes (shared)
    # =====================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEPS 2-3: Structure acquisition + membrane annotation")
    logger.info("=" * 70)

    structures = {}
    membranes = {}
    failed_targets = set()

    for target in all_targets:
        uid = target.uniprot_id
        logger.info("\n--- %s (%s) ---", target.gene_name, uid)

        # Structure
        try:
            structure = acquire_structure(target, str(structures_dir),
                                          force_experimental=force_experimental)
            structures[uid] = structure
            logger.info("  Structure: %s (%s, chain %s)",
                        structure.pdb_id, structure.source, structure.chain_id)
        except (StructureAcquisitionError, Exception) as e:
            logger.error("  Structure FAILED for %s: %s", target.gene_name, e)
            failed_targets.add(uid)
            continue

        # Membrane
        try:
            membrane = annotate_membrane(target, structure)
            membranes[uid] = membrane
        except (MembraneAnnotationError, Exception) as e:
            logger.error("  Membrane FAILED for %s: %s", target.gene_name, e)
            failed_targets.add(uid)
            continue

    # =====================================================================
    # Steps 4-8: Analyze each target in both zones
    # =====================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEPS 4-8: Zone analysis (distal + proximal per target)")
    logger.info("=" * 70)

    zone_results = {}  # {(uid, zone): TargetZoneResult}

    for target in all_targets:
        uid = target.uniprot_id
        if uid in failed_targets:
            continue

        for zone in ("distal", "proximal"):
            logger.info("\n--- %s [%s] ---", target.gene_name, zone.upper())

            if zone == "distal":
                result = _analyze_target_zone(
                    target, structures[uid], membranes[uid],
                    zone="distal", min_distance=distal_dist,
                )
            else:
                result = _analyze_target_zone(
                    target, structures[uid], membranes[uid],
                    zone="proximal", max_distance=proximal_dist,
                )

            zone_results[(uid, zone)] = result

            if result.scores:
                logger.info("  %s [%s]: %d patches, best score %.3f",
                            target.gene_name, zone, len(result.scores),
                            result.best_score)
            else:
                logger.info("  %s [%s]: 0 qualifying patches",
                            target.gene_name, zone)

    # =====================================================================
    # Pair evaluation
    # =====================================================================
    logger.info("\n" + "=" * 70)
    logger.info("PAIR EVALUATION")
    logger.info("=" * 70)

    pair_results = []

    for id_a, id_b in pairs:
        target_a = target_lookup.get(id_a.upper(), target_lookup.get(id_a))
        target_b = target_lookup.get(id_b.upper(), target_lookup.get(id_b))

        if target_a is None or target_b is None:
            logger.error("  Cannot find target(s) for pair %s x %s", id_a, id_b)
            continue

        uid_a = target_a.uniprot_id
        uid_b = target_b.uniprot_id

        if uid_a in failed_targets or uid_b in failed_targets:
            logger.error("  Skipping %s x %s — structure/membrane failed",
                         target_a.gene_name, target_b.gene_name)
            continue

        logger.info("\n--- Pair: %s x %s ---", target_a.gene_name, target_b.gene_name)

        # Orientation AB: A=distal, B=proximal
        zone_a_distal = zone_results.get((uid_a, "distal"))
        zone_b_proximal = zone_results.get((uid_b, "proximal"))
        orientation_ab = _score_orientation(zone_a_distal, zone_b_proximal)

        # Orientation BA: A=proximal, B=distal
        zone_a_proximal = zone_results.get((uid_a, "proximal"))
        zone_b_distal = zone_results.get((uid_b, "distal"))
        orientation_ba = _score_orientation(zone_b_distal, zone_a_proximal)

        both_valid = orientation_ab.is_valid and orientation_ba.is_valid

        # Pick best orientation
        if orientation_ab.orientation_score >= orientation_ba.orientation_score:
            best = orientation_ab
        else:
            best = orientation_ba

        # Apply flexibility bonus
        if both_valid:
            final_score = best.orientation_score * config.FLEXIBILITY_BONUS
        else:
            final_score = best.orientation_score

        pair_result = BispecificPairResult(
            target_a=target_a,
            target_b=target_b,
            orientation_ab=orientation_ab,
            orientation_ba=orientation_ba,
            both_valid=both_valid,
            best_orientation=best,
            final_pair_score=final_score,
        )
        pair_results.append(pair_result)

        # Log results
        logger.info("  Orientation: %s(distal) x %s(proximal) — score %.3f %s",
                     target_a.gene_name, target_b.gene_name,
                     orientation_ab.orientation_score,
                     "[VALID]" if orientation_ab.is_valid else "[no patches]")
        logger.info("  Orientation: %s(distal) x %s(proximal) — score %.3f %s",
                     target_b.gene_name, target_a.gene_name,
                     orientation_ba.orientation_score,
                     "[VALID]" if orientation_ba.is_valid else "[no patches]")
        if both_valid:
            logger.info("  DUAL-VALID — flexibility bonus %.0f%% applied",
                         (config.FLEXIBILITY_BONUS - 1) * 100)
        logger.info("  Final pair score: %.3f", final_score)

    # Sort pairs by final score descending
    pair_results.sort(key=lambda p: p.final_pair_score, reverse=True)

    # =====================================================================
    # Build zone-specific data for exports
    # =====================================================================
    parameters = {
        "mode": "bispecific",
        "distal_min_distance_a": distal_dist,
        "proximal_max_distance_a": proximal_dist,
        "cyno_mismatch_percent_base": config.MAX_CYNO_MISMATCH_PERCENT,
        "cyno_mismatch_scaling": "min(base% * sqrt(n_residues/20), 30%)",
        "nonspecific_percent_base": config.MAX_NONSPECIFIC_PERCENT,
        "nonspecific_rule": "worst single paralog match fraction <= nonspecific_percent_base",
        "vhh_footprint_min_a2": config.VHH_FOOTPRINT_MIN_A2,
        "flexibility_bonus": config.FLEXIBILITY_BONUS,
    }

    # =====================================================================
    # Export: bispecific outputs
    # =====================================================================
    logger.info("\n" + "=" * 70)
    logger.info("EXPORT: Bispecific outputs")
    logger.info("=" * 70)

    from epitope_pipeline.io.export_bispecific import export_bispecific_all
    export_bispecific_all(
        run_dir=str(run_dir),
        pair_results=pair_results,
        zone_results=zone_results,
        structures=structures,
        membranes=membranes,
    )

    # Write input manifest
    import json
    manifest = {
        "mode": "bispecific",
        "pairs": ["{} x {}".format(a, b) for a, b in pairs],
        "parameters": parameters,
        "targets": [
            {"gene_name": t.gene_name, "uniprot_id": t.uniprot_id}
            for t in all_targets
        ],
    }
    manifest_path = supp_dir / "Logs" / "input_manifest.json"
    with open(str(manifest_path), "w") as f:
        json.dump(manifest, f, indent=2)
    logger.info("  Wrote input manifest: input_manifest.json")

    # =====================================================================
    # Visualization: zone-specific single-target maps (4 maps)
    # =====================================================================
    logger.info("\n" + "=" * 70)
    logger.info("VISUALIZATION: Zone-specific epitope maps")
    logger.info("=" * 70)

    for (uid, zone_name), zr in zone_results.items():
        if zone_name == "distal":
            zone_label = "DISTAL (>={:.0f}A)".format(distal_dist)
        else:
            zone_label = "PROXIMAL (<={:.0f}A)".format(proximal_dist)

        plot_epitope_map(
            target=zr.target,
            membrane=zr.membrane,
            spatial_filter=zr.spatial_filter,
            surface_analysis=zr.surface_analysis,
            conservation_result=zr.conservation_result,
            specificity_result=zr.specificity_result,
            scores=zr.scores,
            output_path=str(zone_details_dir / "{}_{}_epitope_map.png".format(
                zr.target.gene_name.lower(), zone_name
            )),
            title_suffix=zone_label,
        )

    # BLAST off-target dot plots (one per unique target)
    plotted_blast = set()
    for key, zr in zone_results.items():
        uid = key[0]
        if uid in plotted_blast:
            continue
        if zr.specificity_result:
            blast_fig_dir = figures_dir / "BLAST"
            blast_fig_dir.mkdir(exist_ok=True)
            plot_blast_offtargets(
                target=zr.target,
                specificity_result=zr.specificity_result,
                output_path=str(blast_fig_dir / "{}_blast_offtargets.png".format(
                    zr.target.gene_name.lower()
                )),
            )
            plotted_blast.add(uid)

    # =====================================================================
    # Visualization: bispecific maps
    # =====================================================================
    logger.info("\n" + "=" * 70)
    logger.info("VISUALIZATION: Bispecific epitope maps")
    logger.info("=" * 70)

    from epitope_pipeline.visualize_bispecific import (
        plot_bispecific_epitope_map, plot_bispecific_summary,
        plot_bispecific_combined,
    )

    for pair_result in pair_results:
        for orientation in (pair_result.orientation_ab, pair_result.orientation_ba):
            out_path = str(figures_dir / "bispecific_{}.png".format(
                orientation.label.lower()
            ))
            plot_bispecific_epitope_map(
                orientation=orientation,
                output_path=out_path,
            )

    # Combined orientation figures (vertical layout, domains + epitope only)
    for pair_result in pair_results:
        out_path = str(figures_dir / "bispecific_combined_{}_x_{}.png".format(
            pair_result.target_a.gene_name.lower(),
            pair_result.target_b.gene_name.lower(),
        ))
        plot_bispecific_combined(
            pair_result=pair_result,
            output_path=out_path,
        )

    if pair_results:
        plot_bispecific_summary(
            pair_results=pair_results,
            output_path=str(figures_dir / "bispecific_pair_summary.png"),
        )

    # =====================================================================
    # Summary
    # =====================================================================
    logger.info("\n" + "=" * 70)
    logger.info("BISPECIFIC PIPELINE COMPLETE")
    logger.info("=" * 70)
    logger.info("Run directory: %s", run_dir)
    for pr in pair_results:
        logger.info(
            "  %s x %s: pair score %.3f%s",
            pr.target_a.gene_name, pr.target_b.gene_name,
            pr.final_pair_score,
            " [DUAL-VALID]" if pr.both_valid else "",
        )

    logger.removeHandler(fh)
    fh.close()

    return {
        "run_dir": str(run_dir),
        "pair_results": pair_results,
        "zone_results": zone_results,
    }


# ---------------------------------------------------------------------------
# Internal analysis
# ---------------------------------------------------------------------------

def _analyze_target_zone(target, structure, membrane, zone,
                         min_distance=None, max_distance=None):
    """
    Run Steps 4-8 of the epitope pipeline for one target in one spatial zone.
    """
    # Extract CA coordinates once — passed to all downstream steps
    ca_coords = extract_ca_coords(structure.pdb_path, structure.chain_id)

    # Step 4: Spatial filter
    if zone == "proximal":
        spatial = filter_ectodomain(target, structure, membrane,
                                     max_distance=max_distance,
                                     ca_coords=ca_coords)
    else:
        spatial = filter_ectodomain(target, structure, membrane,
                                     min_distance=min_distance,
                                     ca_coords=ca_coords)

    # Step 5: Surface analysis
    surface = analyze_surface(target, structure, spatial, membrane,
                               ca_coords=ca_coords)

    # Step 6: Conservation
    conservation = None
    if target.cyno_sequence:
        try:
            conservation = analyze_conservation(target, surface, ca_coords)
        except (ConservationError, Exception) as e:
            logger.error("  Conservation FAILED for %s [%s]: %s",
                         target.gene_name, zone, e)

    # Step 7: Specificity
    specificity = None
    try:
        ec_patches = None
        if surface.residue_sasa:
            ec_patches = cluster_ectodomain_patches(
                target, structure, membrane, surface.residue_sasa,
                ca_coords=ca_coords,
            )

        if conservation:
            specificity = filter_specificity(
                target, conservation,
                ectodomain_patches=ec_patches, ca_coords=ca_coords)
        else:
            empty_cons = ConservationResult(
                uniprot_id=target.uniprot_id,
                gene_name=target.gene_name,
                alignment_human="",
                alignment_cyno="",
                overall_identity=0.0,
                residue_conservation={},
                conserved_patches=[],
                rejected_patches=[],
                patch_conservation={},
            )
            specificity = filter_specificity(
                target, empty_cons,
                ectodomain_patches=ec_patches, ca_coords=ca_coords)
    except Exception as e:
        logger.error("  Specificity FAILED for %s [%s]: %s",
                     target.gene_name, zone, e)

    # Step 8: Scoring
    scores = []
    if (spatial.qualifying_residues and surface.patches and
            conservation and conservation.conserved_patches):
        specific_patches = specificity.specific_patches if specificity else []
        if specific_patches:
            scores = score_epitopes(
                target, specific_patches, conservation, specificity,
                spatial, surface, distance_mode=zone,
            )

    best_score = scores[0].composite_score if scores else 0.0

    return TargetZoneResult(
        target=target,
        zone=zone,
        structure=structure,
        membrane=membrane,
        spatial_filter=spatial,
        surface_analysis=surface,
        conservation_result=conservation,
        specificity_result=specificity,
        scores=scores,
        best_score=best_score,
    )


def _score_orientation(distal_zone, proximal_zone):
    """Score one orientation as geometric mean of best scores."""
    if distal_zone is None or proximal_zone is None:
        return OrientationResult(
            distal_zone=distal_zone,
            proximal_zone=proximal_zone,
            label="unknown",
            orientation_score=0.0,
            is_valid=False,
        )

    a = distal_zone.best_score
    b = proximal_zone.best_score
    is_valid = a > 0 and b > 0
    score = (a * b) ** 0.5 if is_valid else 0.0

    label = "{}_distal__{}_proximal".format(
        distal_zone.target.gene_name.lower(),
        proximal_zone.target.gene_name.lower(),
    )

    return OrientationResult(
        distal_zone=distal_zone,
        proximal_zone=proximal_zone,
        label=label,
        orientation_score=score,
        is_valid=is_valid,
    )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Bispecific epitope pipeline")
    parser.add_argument("pairs", nargs="+",
                        help="TARGET_A:TARGET_B pairs (colon-separated)")
    parser.add_argument("--distal", type=float, default=None,
                        help="Distal min distance (A), default 60")
    parser.add_argument("--proximal", type=float, default=None,
                        help="Proximal max distance (A), default 40")
    args = parser.parse_args()

    pairs = []
    for arg in args.pairs:
        if ":" not in arg:
            print("Error: pairs must be colon-separated (e.g. ERBB2:MSLN)")
            sys.exit(1)
        a, b = arg.split(":", 1)
        pairs.append((a.strip(), b.strip()))

    results = run_bispecific(
        pairs,
        distal_min_distance_a=args.distal,
        proximal_max_distance_a=args.proximal,
    )

    print("\nDone: {}".format(results["run_dir"]))
    for pr in results.get("pair_results", []):
        print("  {} x {}: pair score {:.3f}{}".format(
            pr.target_a.gene_name, pr.target_b.gene_name,
            pr.final_pair_score,
            " [DUAL-VALID]" if pr.both_valid else "",
        ))


if __name__ == "__main__":
    main()
