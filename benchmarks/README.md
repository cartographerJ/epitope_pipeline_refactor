# Benchmarks

Baseline timing reference for the epitope pipeline. The harness runs the full
per-target pipeline (structure → membrane → spatial → surface → conservation →
specificity → scoring) on a small fixture set and writes a JSON report.

## Running

```bash
.venv/bin/python benchmarks/run_baseline.py                  # default fixtures (ERBB2, EGFR)
.venv/bin/python benchmarks/run_baseline.py ERBB2 EGFR ITGB6 # custom
```

Output: `benchmarks/baseline_YYYY-MM-DD_<git-sha>.json`

## What's recorded

Per target, per stage, wall-clock seconds (`time.perf_counter`). Plus
sequence length, qualifying-residue / patch counts at each filter stage,
git commit, Python version, OS, and `USE_LOCAL_BLAST` flag.

## Comparing perf branches

```bash
git checkout perf/<branch>
.venv/bin/python benchmarks/run_baseline.py
diff <(jq -S . benchmarks/baseline_*_<old-sha>.json) \
     <(jq -S . benchmarks/baseline_*_<new-sha>.json)
```

Commit the JSON for each perf branch's run alongside the code change so the
delta is reviewable.

## Caveats

- **First run is slow.** Cold caches mean structure downloads (RCSB / AlphaFold
  DB) and BLAST queries against the human proteome. Subsequent runs reuse
  `cache/blast/` and `benchmarks/.structures/`.
- **BLAST dominates.** Without a local BLAST DB (`USE_LOCAL_BLAST=False`) the
  benchmark is mostly measuring NCBI rate-limited remote calls, which is not a
  meaningful signal for code-level perf work. Set up local BLAST via
  `bash scripts/setup_blast_db.sh` before establishing a serious baseline.
- **Network variability.** On a fresh run the structure-acquisition stage
  includes UniProt + RCSB + AlphaFold DB API latency. Re-runs with the
  structure cache populated isolate compute time.
