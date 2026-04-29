"""Baseline benchmark for the epitope pipeline.

Times each pipeline stage on a small fixture set and dumps the result to
``benchmarks/baseline_YYYY-MM-DD_<commit>.json``. This is the ground-truth
reference for measuring perf-branch wins — re-run on each perf branch and
diff the JSONs.

Usage:
    .venv/bin/python benchmarks/run_baseline.py
    .venv/bin/python benchmarks/run_baseline.py ERBB2 EGFR

Notes:
    First run on a fresh checkout fetches structures from AlphaFold DB / RCSB
    and queries BLAST. With local BLAST + cached structures a full run is
    typically seconds to a few minutes; without a local BLAST DB or with cold
    caches it can take 5-30 minutes per target.
"""

import json
import platform
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

from epitope_pipeline import config
from epitope_pipeline.io.targets import resolve_targets
from epitope_pipeline.io.structure import acquire_structure
from epitope_pipeline.io.membrane import annotate_membrane
from epitope_pipeline.io.pdb import extract_ca_coords
from epitope_pipeline.compute.spatial import filter_ectodomain
from epitope_pipeline.compute.surface import (
    analyze_surface, cluster_ectodomain_patches,
)
from epitope_pipeline.compute.conservation import (
    analyze_conservation, ConservationResult,
)
from epitope_pipeline.compute.specificity import filter_specificity
from epitope_pipeline.compute.scoring import score_epitopes


DEFAULT_FIXTURES = ["ERBB2", "EGFR"]
REPO = Path(__file__).resolve().parent.parent


def _git_short_sha():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short=7", "HEAD"], cwd=REPO,
        ).decode().strip()
    except Exception:
        return "unknown"


def _time(fn, *args, **kwargs):
    """Run fn and return (result, elapsed_seconds)."""
    t0 = time.perf_counter()
    out = fn(*args, **kwargs)
    return out, time.perf_counter() - t0


def benchmark_target(target, structures_dir):
    """Run the per-target pipeline body and return a stage timing dict."""
    timings = {}

    structure, t = _time(acquire_structure, target, str(structures_dir))
    timings["structure"] = t

    ca_coords, t = _time(extract_ca_coords, structure.pdb_path, structure.chain_id)
    timings["ca_coords"] = t

    membrane, t = _time(annotate_membrane, target, structure)
    timings["membrane"] = t

    spatial, t = _time(
        filter_ectodomain, target, structure, membrane, ca_coords=ca_coords,
    )
    timings["spatial"] = t

    surface, t = _time(
        analyze_surface, target, structure, spatial, membrane,
        ca_coords=ca_coords,
    )
    timings["surface"] = t

    if target.cyno_sequence:
        cons, t = _time(analyze_conservation, target, surface, ca_coords)
        timings["conservation"] = t
    else:
        cons = ConservationResult(
            uniprot_id=target.uniprot_id, gene_name=target.gene_name,
            alignment_human="", alignment_cyno="",
            overall_identity=0.0, residue_conservation={},
            conserved_patches=[], rejected_patches=[], patch_conservation={},
        )
        timings["conservation"] = 0.0

    ec_patches, t = _time(
        cluster_ectodomain_patches, target, structure, membrane,
        surface.residue_sasa, ca_coords=ca_coords,
    )
    timings["ec_patches"] = t

    spec, t = _time(
        filter_specificity, target, cons,
        ectodomain_patches=ec_patches, ca_coords=ca_coords,
    )
    timings["specificity"] = t

    scoring_t = 0.0
    n_scored = 0
    if spec.specific_patches:
        scores, scoring_t = _time(
            score_epitopes, target, spec.specific_patches, cons, spec,
            spatial, surface,
        )
        n_scored = len(scores)
    timings["scoring"] = scoring_t

    timings["total"] = sum(timings.values())

    return {
        "target": target.gene_name,
        "uniprot_id": target.uniprot_id,
        "sequence_length": target.sequence_length,
        "n_qualifying_residues": spatial.total_qualifying,
        "n_surface_patches": len(surface.patches) if surface else 0,
        "n_conserved_patches": len(cons.conserved_patches) if cons else 0,
        "n_specific_patches": len(spec.specific_patches),
        "n_scored": n_scored,
        "stages_s": {k: round(v, 4) for k, v in timings.items()},
    }


def main():
    fixtures = sys.argv[1:] or DEFAULT_FIXTURES

    print("Benchmarking on fixtures:", fixtures)
    print("USE_LOCAL_BLAST:", config.USE_LOCAL_BLAST)

    targets, t_resolve = _time(resolve_targets, fixtures)
    print("resolve_targets:        {:8.3f}s".format(t_resolve))

    structures_dir = REPO / "benchmarks" / ".structures"
    structures_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for target in targets:
        print("\n--- {} ({}) ---".format(target.gene_name, target.uniprot_id))
        try:
            r = benchmark_target(target, structures_dir)
            results.append(r)
            for stage, t in r["stages_s"].items():
                print("  {:18s} {:8.3f}s".format(stage, t))
        except Exception as e:
            print("  FAILED:", e)
            results.append({
                "target": target.gene_name,
                "uniprot_id": target.uniprot_id,
                "error": str(e),
            })

    sha = _git_short_sha()
    output = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "git_commit": sha,
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "fixtures": fixtures,
        "use_local_blast": config.USE_LOCAL_BLAST,
        "resolve_targets_s": round(t_resolve, 4),
        "results": results,
    }

    out_dir = REPO / "benchmarks"
    out_dir.mkdir(exist_ok=True)
    date = datetime.now().strftime("%Y-%m-%d")
    out_path = out_dir / "baseline_{}_{}.json".format(date, sha)
    out_path.write_text(json.dumps(output, indent=2))
    print("\nWrote", out_path)


if __name__ == "__main__":
    main()
