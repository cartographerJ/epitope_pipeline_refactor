# Single-protein ECD suitability mode — design

**Date:** 2026-05-21
**Status:** Approved (Approach A)

## Problem

The package is perceived as built for *target pairs*. The `bispecific.py` entry
point evaluates each protein in two membrane-distance zones (distal ≥60 Å and
proximal ≤50 Å) so a bispecific antibody can have one arm reaching high and one
low. Even the single-target `run.py` CLI quietly defaults to **proximal** mode
(`max_distance = PROXIMAL_MAX_DISTANCE_A`), which is bispecific-flavored.

We want a clean, first-class way to ask a simpler question: **is the entire
extracellular domain of one protein a suitable VHH-epitope target?** — with no
distal/proximal zone framing at all.

## Decision

Add a dedicated entry point for single-protein, whole-ectodomain assessment,
parallel to `epitope-bispecific`. Implemented as a **thin wrapper** over the
existing per-target orchestrator (`run_pipeline`), not a forked orchestration
loop.

Rationale: the analysis a single protein needs is identical to `run_pipeline`
with the distance filter switched off. Bispecific needed its own loop because it
runs each target through two zones; single mode does not. A wrapper gives the
same "separate entry point" experience (own module, CLI, console script, output
framing) without duplicating the ~200-line structure→membrane→surface→
conservation→specificity→scoring→export→visualize loop, and it inherits future
fixes to the shared engine automatically.

Out of scope (per brainstorming): no new go/no-go "verdict" output. Suitability
is read off the existing patch-count / score / summary-CSV artifacts.

## Changes

### New file: `epitope_pipeline/single.py`

```python
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
    """Whole-ectodomain suitability assessment for one or more single proteins.

    Delegates to run_pipeline(..., no_distance_filter=True). No membrane-distance
    zones are applied — the full ectodomain surface is considered.
    """
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
```

- **No** `min_distance`/`max_distance`/`no_distance_filter` parameters are
  exposed — removing the distal/proximal concept is the entire point.
- `main()` argparse mirrors `run.py`'s, **minus** the three distance flags
  (`--min-distance`, `--max-distance`, `--no-distance-filter`). Retains:
  positional targets with `=ALIAS`, `--targets-file`, `--run-name`, `--no-cyno`,
  `--cyno-mismatch-percent`, `--nonspecific-percent`, `--force-experimental`,
  `--quiet`. Reuses `run.py`'s `_split_alias`, `_load_targets_file`, `_dedupe`
  helpers (import them from `epitope_pipeline.run`).
- Default run-name slug: `"{YYMMDD_HHMM}_single_{slug}"` (e.g.
  `260521_1800_single_erbb2_egfr`).

### Honesty fixes in `epitope_pipeline/run.py`

The existing code mislabels output when the distance filter is off. Fix as part
of this work:

1. `parameters["mode"]`: when `no_distance_filter` is True, label it
   `"whole ectodomain (no distance filter)"` instead of the current
   `"distal (min_distance_a=80)"` (which is computed because the ternary only
   checks `max_distance_a`).
2. Multi-target `plot_scoring_summary` call: when `no_distance_filter` is True,
   pass a "whole ECD" `distance_label`/`distance_mode` instead of the hardcoded
   `"≥80Å"` / distal so the summary figure is not wrong.

### Wiring

- `pyproject.toml` `[project.scripts]`: add
  `epitope-single = "epitope_pipeline.single:main"`.

### Tests

- Add `tests/test_single_cli.py` (mirroring `test_run_cli.py`): assert the
  single-mode argparser parses targets/aliases and does **not** accept the
  distance flags.
- Add a `run_single` smoke assertion (mirroring `test_pipeline_smoke.py`) that it
  invokes `run_pipeline` with `no_distance_filter=True` — can be a fast unit test
  that monkeypatches `run_pipeline` and inspects kwargs, avoiding network/BLAST.

### Docs

- `README.md`: add a "Single-protein ECD suitability" section alongside the
  existing batch / bispecific sections, showing
  `python -m epitope_pipeline.single ERBB2` / `epitope-single ERBB2`.

## Unchanged

- All compute modules, `bispecific.py`, and `run.py`'s own default behavior
  (still defaults to proximal for backward compatibility).
- All output artifact formats (figures, summary CSV, exports).

## Verification

- `pytest` passes (new + existing).
- `epitope-single ERBB2 --help` shows no distance flags.
- A real `run_single(["ERBB2"])` produces a run dir whose
  `parameters["mode"]` reads `"whole ectodomain (no distance filter)"`.
