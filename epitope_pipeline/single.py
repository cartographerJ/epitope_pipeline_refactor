"""
Single-protein ECD suitability — assess whether one protein's whole
extracellular domain is a viable VHH-epitope target.

Unlike the bispecific pipeline (which scores each target in a distal and a
proximal membrane zone for complementary antibody arms), this entry point
applies NO distance zone: the entire ectodomain surface is considered.

Usage:
    from epitope_pipeline.single import run_single
    results = run_single(["ERBB2"])

    python -m epitope_pipeline.single ERBB2 EGFR
"""

from datetime import datetime

from epitope_pipeline.run import (
    run_pipeline, _split_alias, _load_targets_file, _dedupe,
)


def run_single(
    identifiers,
    run_name=None,
    cyno_mismatch_percent=None,
    skip_cyno_gate=False,
    nonspecific_percent=None,
    force_experimental=False,
    aliases=None,
    verbose=True,
):
    """Whole-ectodomain suitability assessment for single proteins.

    Delegates to run_pipeline with no membrane-distance filter, so the full
    ectodomain surface is evaluated for druggable, cyno-conserved,
    human-specific epitope patches. No distal/proximal zones are applied.

    Args mirror run_pipeline's non-distance arguments, except the deprecated
    cyno_max_mismatches count is intentionally omitted — use
    cyno_mismatch_percent for cyno-conservation tuning. See run_pipeline for
    the structure of the returned dict.
    """
    if run_name is None:
        slug = "_".join(i.lower()[:8] for i in identifiers[:3])
        if len(identifiers) > 3:
            slug += "_etc"
        run_name = "{}_single_{}".format(
            datetime.now().strftime("%y%m%d_%H%M"), slug
        )

    return run_pipeline(
        identifiers,
        run_name=run_name,
        no_distance_filter=True,
        cyno_mismatch_percent=cyno_mismatch_percent,
        skip_cyno_gate=skip_cyno_gate,
        nonspecific_percent=nonspecific_percent,
        force_experimental=force_experimental,
        aliases=aliases,
        verbose=verbose,
    )


def build_argparser():
    import argparse

    parser = argparse.ArgumentParser(
        prog="epitope-single",
        description=(
            "Assess whether a single protein's whole extracellular domain is a "
            "suitable VHH-epitope target. No distal/proximal zones are applied. "
            "Targets can be UniProt IDs (P04626) or gene names (ERBB2)."
        ),
    )
    parser.add_argument("targets", nargs="*", metavar="TARGET",
                        help="UniProt accessions or gene names. Append "
                             "`=ALIAS` (e.g. ERBB2=HER2) to override the "
                             "display name used in figures and CSVs.")
    parser.add_argument("--targets-file", metavar="PATH",
                        help="File with one IDENT[=ALIAS] per line; merges with "
                             "positional args. Blank lines and `#` comments are ignored.")
    parser.add_argument("--run-name", metavar="NAME",
                        help="Custom run-directory name under runs/")
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

    results = run_single(
        identifiers,
        run_name=args.run_name,
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
