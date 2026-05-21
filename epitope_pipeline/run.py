"""
Epitope Pipeline — Main orchestrator and CLI entry point.

Wires all pipeline modules together into an end-to-end workflow:
  target resolution -> structure -> membrane -> spatial -> surface ->
  conservation -> specificity -> scoring -> export -> visualize

Usage:
    # From Python
    from epitope_pipeline.run import run_pipeline
    results = run_pipeline(["ERBB2", "EGFR"])

    # From command line
    python -m epitope_pipeline.run ERBB2 EGFR

Per-target resilience: if one target fails at any step, the error is
logged and the pipeline continues with remaining targets.
"""

import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from epitope_pipeline import config
from epitope_pipeline.utils import setup_logging, empty_metric
from epitope_pipeline.io.pdb import extract_ca_coords
from epitope_pipeline.io.targets import resolve_targets, TargetResolutionError
from epitope_pipeline.io.structure import (
    acquire_structure, StructureAcquisitionError,
)
from epitope_pipeline.io.membrane import annotate_membrane, MembraneAnnotationError
from epitope_pipeline.compute.spatial import filter_ectodomain
from epitope_pipeline.compute.surface import analyze_surface, cluster_ectodomain_patches
from epitope_pipeline.compute.conservation import analyze_conservation, ConservationError
from epitope_pipeline.compute.specificity import filter_specificity
from epitope_pipeline.compute.scoring import score_epitopes, compute_target_epitope_metric
from epitope_pipeline.io.export import (
    export_all, build_summary_dataframe, write_summary_metrics_csv,
)
from epitope_pipeline.viz.visualize import plot_epitope_map, plot_scoring_summary, plot_blast_offtargets


logger = logging.getLogger("epitope_pipeline")


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def _mode_label(no_distance_filter, max_distance_a):
    """Human-readable spatial-mode label for the run parameters block."""
    if no_distance_filter:
        return "whole ectodomain (no distance filter)"
    if max_distance_a:
        return "proximal (max_distance_a={})".format(max_distance_a)
    return "distal (min_distance_a={})".format(config.ECTODOMAIN_MIN_DISTANCE_A)


def _summary_distance_args(no_distance_filter, max_distance_a):
    """Distance kwargs for plot_scoring_summary, honest about whole-ECD runs."""
    if no_distance_filter:
        return {
            "distance_label": "whole ECD",
            "distance_value": None,
            "distance_mode": "whole_ecd",
        }
    if max_distance_a:
        return {
            "distance_label": "\u2264{:.0f}\u00c5".format(max_distance_a),
            "distance_value": max_distance_a,
            "distance_mode": "proximal",
        }
    return {
        "distance_label": "\u2265{:.0f}\u00c5".format(config.ECTODOMAIN_MIN_DISTANCE_A),
        "distance_value": config.ECTODOMAIN_MIN_DISTANCE_A,
        "distance_mode": "distal",
    }


# ---------------------------------------------------------------------------
# Pipeline orchestrator
# ---------------------------------------------------------------------------

def run_pipeline(
    identifiers,
    run_name=None,
    min_distance_a=None,
    max_distance_a=None,
    no_distance_filter=False,
    cyno_max_mismatches=None,
    cyno_mismatch_percent=None,
    skip_cyno_gate=False,
    nonspecific_percent=None,
    force_experimental=False,
    aliases=None,
    verbose=True,
):
    """
    End-to-end epitope identification pipeline.

    Args:
        identifiers: List of UniProt accession IDs or gene names.
        run_name: Optional custom name for the run directory.
        min_distance_a: Override min ectodomain distance (default 80A, distal mode).
        max_distance_a: Override max ectodomain distance (proximal mode).
            When set, pipeline finds epitopes <= max_distance_a from membrane
            instead of >= min_distance_a.
        cyno_max_mismatches: Override max cyno mismatches per 600A² (default 2).
        cyno_mismatch_percent: Override per-patch cyno mismatch threshold (default 15%).
        skip_cyno_gate: If True, bypass the cyno-conservation short-circuit so
            scoring proceeds for all surface patches regardless of cyno divergence.
            Per-residue cyno identity is still computed and reported. Useful for
            viral or non-mammalian targets, or for exploratory runs.
        force_experimental: Use experimental PDB when available (default: AlphaFold).
        aliases: Optional dict mapping any of the input identifier, the resolved
            UniProt accession, or the resolved gene name → a display alias.
            Matching is case-insensitive. The alias overwrites
            ``target.gene_name``, so it propagates through every output (figure
            file names, CSV rows, log lines, scoring summary chart).
        verbose: Whether to log to console (default True).

    Returns:
        Dict with:
          run_dir: Path to output directory
          targets: List of resolved TargetInfo objects
          scores: Dict {uniprot_id: [EpitopeScore, ...]}
          metrics: Dict {uniprot_id: metric_dict}
    """
    # --- Apply overrides ---
    if min_distance_a is not None:
        config.ECTODOMAIN_MIN_DISTANCE_A = min_distance_a
    if cyno_max_mismatches is not None:
        config.MAX_CYNO_MISMATCHES_PER_600A2 = cyno_max_mismatches
    if cyno_mismatch_percent is not None:
        config.MAX_CYNO_MISMATCH_PERCENT = cyno_mismatch_percent
    if nonspecific_percent is not None:
        config.MAX_NONSPECIFIC_PERCENT = nonspecific_percent

    # --- Setup logging ---
    setup_logging(verbose)
    logger.info("=" * 70)
    logger.info("EPITOPE PIPELINE v%s", "0.1.0")
    logger.info("=" * 70)
    logger.info("Targets: %s", ", ".join(identifiers))

    # --- Create run directory (datetime-stamped, minute resolution) ---
    if run_name is None:
        target_slug = "_".join(i.lower()[:8] for i in identifiers[:3])
        if len(identifiers) > 3:
            target_slug += "_etc"
        run_name = "{}_{}".format(
            datetime.now().strftime("%y%m%d_%H%M"), target_slug
        )

    run_dir = config.RUNS_DIR / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    structures_dir = run_dir / ".structures"
    structures_dir.mkdir(exist_ok=True)
    figures_dir = run_dir / "Figures"
    figures_dir.mkdir(exist_ok=True)

    # Add file handler for logging
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

    # Record parameters used
    parameters = {
        "mode": _mode_label(no_distance_filter, max_distance_a),
        "min_distance_a": config.ECTODOMAIN_MIN_DISTANCE_A,
        "max_distance_a": max_distance_a,
        "cyno_mismatch_percent_base": config.MAX_CYNO_MISMATCH_PERCENT,
        "cyno_mismatch_scaling": "min(base% * sqrt(n_residues/20), 30%)",
        "cyno_gate_skipped": skip_cyno_gate,
        "nonspecific_percent_base": config.MAX_NONSPECIFIC_PERCENT,
        "nonspecific_rule": "worst single paralog match fraction <= nonspecific_percent_base",
        "vhh_footprint_min_a2": config.VHH_FOOTPRINT_MIN_A2,
        "patch_clustering_distance_a": config.PATCH_CLUSTERING_DISTANCE_A,
        "surface_exposure_threshold": config.SURFACE_EXPOSURE_THRESHOLD,
    }

    # =====================================================================
    # Step 1: Resolve targets
    # =====================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 1: Resolving targets")
    logger.info("=" * 70)

    try:
        targets = resolve_targets(identifiers)
    except TargetResolutionError as e:
        logger.error("Target resolution failed: %s", e)
        return {"run_dir": str(run_dir), "targets": [], "scores": {}, "metrics": {}}

    logger.info("Resolved %d targets", len(targets))

    # Apply display aliases: overwrite target.gene_name so the alias flows
    # through every downstream output (figures, CSV rows, file names, logs).
    # Match by input identifier, resolved gene name, or UniProt accession —
    # case-insensitive — so the caller can key on whichever they have handy.
    if aliases:
        alias_lookup = {k.upper(): v for k, v in aliases.items()}
        for ident, target in zip(identifiers, targets):
            alias = (
                alias_lookup.get(ident.upper())
                or alias_lookup.get(target.gene_name.upper())
                or alias_lookup.get(target.uniprot_id.upper())
            )
            if alias and alias != target.gene_name:
                logger.info("  Alias: %s (%s) -> %s",
                            target.gene_name, target.uniprot_id, alias)
                target.gene_name = alias

    # Per-target pipeline state
    all_structures = {}
    all_membranes = {}
    all_spatial = {}
    all_surface = {}
    all_conservation = {}
    all_specificity = {}
    all_scores = {}
    all_metrics = {}

    for target in targets:
        uid = target.uniprot_id
        logger.info("\n" + "-" * 70)
        logger.info("Processing: %s (%s)", target.gene_name, uid)
        logger.info("-" * 70)

        # =================================================================
        # Step 2: Acquire structure
        # =================================================================
        logger.info("\nStep 2: Structure acquisition")
        try:
            structure = acquire_structure(target, str(structures_dir),
                                          force_experimental=force_experimental)
            all_structures[uid] = structure
            logger.info(
                "  Structure: %s (%s, chain %s)",
                structure.pdb_id, structure.source, structure.chain_id,
            )
        except (StructureAcquisitionError, Exception) as e:
            logger.error("  Structure acquisition FAILED for %s: %s", target.gene_name, e)
            all_metrics[uid] = empty_metric()
            continue

        # Extract CA coordinates once — passed to all downstream steps
        ca_coords = extract_ca_coords(structure.pdb_path, structure.chain_id)

        # =================================================================
        # Step 3: Membrane annotation
        # =================================================================
        logger.info("\nStep 3: Membrane annotation")
        try:
            membrane = annotate_membrane(target, structure)
            all_membranes[uid] = membrane
        except (MembraneAnnotationError, Exception) as e:
            logger.error("  Membrane annotation FAILED for %s: %s", target.gene_name, e)
            all_metrics[uid] = empty_metric()
            continue

        # =================================================================
        # Step 4: Ectodomain spatial filter
        # =================================================================
        if no_distance_filter:
            logger.info("\nStep 4: Ectodomain distance filter (OFF — all ECD residues)")
            spatial = filter_ectodomain(target, structure, membrane,
                                         no_distance_filter=True,
                                         ca_coords=ca_coords)
        elif max_distance_a is not None:
            logger.info("\nStep 4: Ectodomain distance filter (<= %.0fA, proximal)", max_distance_a)
            spatial = filter_ectodomain(target, structure, membrane,
                                         max_distance=max_distance_a,
                                         ca_coords=ca_coords)
        else:
            logger.info("\nStep 4: Ectodomain distance filter (>= %.0fA)", config.ECTODOMAIN_MIN_DISTANCE_A)
            spatial = filter_ectodomain(target, structure, membrane,
                                         ca_coords=ca_coords)
        all_spatial[uid] = spatial

        # =================================================================
        # Step 5: Surface analysis
        # =================================================================
        # Always run surface analysis — needed for SASA track in figure
        # even when no residues qualify (returns empty SurfaceAnalysis).
        has_qualifying = bool(spatial.qualifying_residues)
        logger.info("\nStep 5: Surface analysis (SASA + patch clustering)")
        surface = analyze_surface(target, structure, spatial, membrane,
                                   ca_coords=ca_coords)
        all_surface[uid] = surface

        if not has_qualifying:
            logger.info("  No residues qualify — %s has no epitope space at %.0fA",
                        target.gene_name, config.ECTODOMAIN_MIN_DISTANCE_A)

        # =================================================================
        # Step 6: Cyno conservation
        # =================================================================
        # Always run conservation when cyno ortholog exists — the per-residue
        # data is shown in the figure regardless of whether patches qualify.
        if target.cyno_sequence:
            logger.info("\nStep 6: Cynomolgus monkey conservation")
            try:
                conservation = analyze_conservation(target, surface, ca_coords)
                all_conservation[uid] = conservation
            except (ConservationError, Exception) as e:
                logger.error("  Conservation analysis FAILED for %s: %s", target.gene_name, e)
        else:
            logger.info("\nStep 6: No cynomolgus ortholog available for %s", target.gene_name)

        # =================================================================
        # Step 7: Human proteome specificity
        # =================================================================
        # Always run full-sequence BLAST for per-residue specificity in the
        # figure, even if no patches pass conservation.
        logger.info("\nStep 7: Human proteome specificity (BLAST)")
        try:
            conservation_for_spec = all_conservation.get(uid)
            ec_patches = None
            if surface.residue_sasa:
                ec_patches = cluster_ectodomain_patches(
                    target, structure, membrane, surface.residue_sasa,
                    ca_coords=ca_coords,
                )
            from epitope_pipeline.compute.conservation import ConservationResult
            from dataclasses import replace

            if conservation_for_spec is None:
                # No cyno ortholog — synthesize an empty result so BLAST still runs
                conservation_for_spec = ConservationResult(
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

            # With --no-cyno, widen the candidate set to every surface patch so
            # cyno-divergent patches still reach scoring. Per-residue cyno
            # identity (residue_conservation / patch_conservation) is preserved.
            if skip_cyno_gate and surface and surface.patches:
                conservation_for_spec = replace(
                    conservation_for_spec,
                    conserved_patches=list(surface.patches),
                )

            specificity = filter_specificity(
                target, conservation_for_spec,
                ectodomain_patches=ec_patches, ca_coords=ca_coords)
            all_specificity[uid] = specificity
        except Exception as e:
            logger.error("  Specificity screening FAILED for %s: %s", target.gene_name, e)

        # --- Check if we can proceed to scoring ---
        if not has_qualifying:
            all_metrics[uid] = empty_metric()
            continue

        surface = all_surface.get(uid)
        if not surface or not surface.patches:
            if surface and not surface.patches:
                logger.info("  No surface patches >= %.0f A²", config.VHH_FOOTPRINT_MIN_A2)
            all_metrics[uid] = empty_metric()
            continue

        conservation = all_conservation.get(uid)
        if not conservation or not conservation.conserved_patches:
            if not skip_cyno_gate:
                if conservation:
                    logger.info("  No patches pass cyno conservation (max %.1f%% mismatches)",
                                config.MAX_CYNO_MISMATCH_PERCENT)
                all_metrics[uid] = empty_metric()
                continue
            # skip_cyno_gate: fall through with the gate bypassed. Scoring
            # needs a conservation object — synthesize one if absent.
            if conservation is None:
                from epitope_pipeline.compute.conservation import ConservationResult
                conservation = ConservationResult(
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
                all_conservation[uid] = conservation
            logger.info("  --no-cyno: bypassing cyno conservation gate")

        specific_patches = all_specificity[uid].specific_patches if uid in all_specificity else []
        if not specific_patches:
            logger.info("  No patches pass specificity screen")
            all_metrics[uid] = empty_metric()
            continue

        # =================================================================
        # Step 8: Scoring
        # =================================================================
        logger.info("\nStep 8: Epitope scoring")
        scores = score_epitopes(
            target, specific_patches, conservation, specificity, spatial, surface,
        )
        all_scores[uid] = scores

        metric = compute_target_epitope_metric(scores, surface)
        all_metrics[uid] = metric
        logger.info(
            "  %s: %d patches, total epitope area %.0f A², best score %.3f",
            target.gene_name, metric["n_patches"],
            metric["total_epitope_area_a2"], metric["best_score"],
        )

    # =====================================================================
    # Step 9: Export all outputs
    # =====================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 9: Exporting results")
    logger.info("=" * 70)

    # Fill in empty metrics for targets that weren't scored
    for target in targets:
        uid = target.uniprot_id
        if uid not in all_metrics:
            all_metrics[uid] = empty_metric()

    # Build the gene × info summary DataFrame (one row per input target,
    # including failures) and write it as a pandas-readable CSV.
    metrics_df = build_summary_dataframe(
        targets=targets,
        target_metrics=all_metrics,
        surface_analyses=all_surface,
        conservation_results=all_conservation,
        specificity_results=all_specificity,
        epitope_scores=all_scores,
        membranes=all_membranes,
        parameters=parameters,
        run_name=run_name,
    )
    summary_csv_path = write_summary_metrics_csv(supp_dir, metrics_df)

    export_all(
        run_dir=str(run_dir),
        targets=targets,
        epitope_scores=all_scores,
        structures=all_structures,
        membranes=all_membranes,
        spatial_filters=all_spatial,
        surface_analyses=all_surface,
        conservation_results=all_conservation,
        specificity_results=all_specificity,
        target_metrics=all_metrics,
        parameters=parameters,
    )

    # =====================================================================
    # Step 10: Visualization
    # =====================================================================
    logger.info("\n" + "=" * 70)
    logger.info("STEP 10: Generating figures")
    logger.info("=" * 70)

    for target in targets:
        uid = target.uniprot_id
        if uid in all_membranes and uid in all_spatial:
            plot_epitope_map(
                target=target,
                membrane=all_membranes[uid],
                spatial_filter=all_spatial[uid],
                surface_analysis=all_surface.get(uid),
                conservation_result=all_conservation.get(uid),
                specificity_result=all_specificity.get(uid),
                scores=all_scores.get(uid, []),
                output_path=str(figures_dir / "{}_epitope_map.png".format(
                    target.gene_name.lower()
                )),
            )

    # BLAST off-target dot plots
    for target in targets:
        uid = target.uniprot_id
        if all_specificity.get(uid):
            blast_fig_dir = figures_dir / "BLAST"
            blast_fig_dir.mkdir(exist_ok=True)
            plot_blast_offtargets(
                target=target,
                specificity_result=all_specificity[uid],
                output_path=str(blast_fig_dir / "{}_blast_offtargets.png".format(
                    target.gene_name.lower()
                )),
            )

    # Multi-target summary
    if len(targets) > 1:
        plot_scoring_summary(
            all_scores=all_scores,
            target_metrics=all_metrics,
            targets=targets,
            output_path=str(figures_dir / "scoring_summary.png"),
            **_summary_distance_args(no_distance_filter, max_distance_a),
        )

    # =====================================================================
    # Summary
    # =====================================================================
    logger.info("\n" + "=" * 70)
    logger.info("PIPELINE COMPLETE")
    logger.info("=" * 70)
    logger.info("Run directory: %s", run_dir)
    for target in targets:
        uid = target.uniprot_id
        metric = all_metrics.get(uid, {})
        n = metric.get("n_patches", 0)
        area = metric.get("total_epitope_area_a2", 0.0)
        best = metric.get("best_score", 0.0)
        logger.info(
            "  %s: %d epitope patches, %.0f A² total, best score %.3f",
            target.gene_name, n, area, best,
        )

    # Remove file handler to avoid issues on re-runs
    logger.removeHandler(fh)
    fh.close()

    return {
        "run_dir": str(run_dir),
        "targets": targets,
        "scores": all_scores,
        "metrics": all_metrics,
        "metrics_df": metrics_df,
        "summary_csv_path": str(summary_csv_path),
    }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _split_alias(token):
    """Split ``IDENT=ALIAS`` into ``(ident, alias)`` (``alias`` may be None).

    Whitespace around the delimiter is stripped. An empty alias after ``=`` is
    treated as no alias.
    """
    if "=" in token:
        ident, _, alias = token.partition("=")
        ident = ident.strip()
        alias = alias.strip() or None
        return ident, alias
    return token.strip(), None


def _load_targets_file(path):
    """Read one ``IDENT[=ALIAS]`` per line; ignore blank lines and `#` comments."""
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"--targets-file not found: {path}")
    out = []
    for line in p.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line)
    return out


def _dedupe(items):
    """Preserve order, drop later duplicates."""
    seen = set()
    out = []
    for it in items:
        if it not in seen:
            seen.add(it)
            out.append(it)
    return out


def build_argparser():
    import argparse
    parser = argparse.ArgumentParser(
        prog="epitope-pipeline",
        description=(
            "Find druggable VHH epitopes on human membrane proteins. "
            "Targets can be UniProt IDs (P04626) or gene names (ERBB2)."
        ),
    )
    parser.add_argument("targets", nargs="*", metavar="TARGET",
                        help="UniProt accessions or gene names. Append "
                             "`=ALIAS` (e.g. ERBB2=HER2) to override the "
                             "display name used in figures and CSVs.")
    parser.add_argument("--targets-file", metavar="PATH",
                        help="File with one IDENT[=ALIAS] per line; merges with positional args. "
                             "Blank lines and `#` comments are ignored.")
    parser.add_argument("--run-name", metavar="NAME",
                        help="Custom run-directory name under runs/")
    parser.add_argument("--min-distance", type=float, dest="min_distance_a",
                        help="Distal mode: minimum distance from membrane (Å); default 80")
    parser.add_argument("--max-distance", type=float, dest="max_distance_a",
                        help="Proximal mode: maximum distance from membrane (Å)")
    parser.add_argument("--no-distance-filter", action="store_true",
                        help="Skip the ectodomain distance filter entirely")
    parser.add_argument("--cyno-mismatch-percent", type=float,
                        dest="cyno_mismatch_percent",
                        help="Per-patch cyno mismatch tolerance (default 15)")
    parser.add_argument("--no-cyno", action="store_true", dest="skip_cyno_gate",
                        help="Bypass the cyno-conservation gate entirely "
                             "(per-residue cyno identity is still reported)")
    parser.add_argument("--nonspecific-percent", type=float,
                        dest="nonspecific_percent",
                        help="Worst-paralog match fraction allowed (default 85)")
    parser.add_argument("--force-experimental", action="store_true",
                        help="Prefer experimental PDB structures over AlphaFold")
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress console logging (file log still written)")
    return parser


def main():
    parser = build_argparser()
    args = parser.parse_args()

    raw_tokens = list(args.targets)
    if args.targets_file:
        raw_tokens.extend(_load_targets_file(args.targets_file))

    identifiers = []
    aliases = {}
    for token in raw_tokens:
        ident, alias = _split_alias(token)
        if not ident:
            continue
        identifiers.append(ident)
        if alias:
            aliases[ident] = alias
    identifiers = _dedupe(identifiers)

    if not identifiers:
        parser.error("no targets provided (pass positional args or --targets-file)")

    # Default to proximal mode unless the user picks distal explicitly
    max_distance = args.max_distance_a
    if max_distance is None and args.min_distance_a is None and not args.no_distance_filter:
        max_distance = config.PROXIMAL_MAX_DISTANCE_A

    results = run_pipeline(
        identifiers,
        run_name=args.run_name,
        min_distance_a=args.min_distance_a,
        max_distance_a=max_distance,
        no_distance_filter=args.no_distance_filter,
        cyno_mismatch_percent=args.cyno_mismatch_percent,
        skip_cyno_gate=args.skip_cyno_gate,
        nonspecific_percent=args.nonspecific_percent,
        force_experimental=args.force_experimental,
        aliases=aliases or None,
        verbose=not args.quiet,
    )

    print("\nDone! Results in: {}".format(results["run_dir"]))
    for target in results.get("targets", []):
        uid = target.uniprot_id
        metric = results.get("metrics", {}).get(uid, {})
        print("  {}: {} patches, {:.0f} A² epitope space".format(
            target.gene_name,
            metric.get("n_patches", 0),
            metric.get("total_epitope_area_a2", 0.0),
        ))


if __name__ == "__main__":
    main()
