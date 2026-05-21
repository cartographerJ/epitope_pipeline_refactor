# Tunable cyno conservation mode (faithful old local test as default) — design

**Date:** 2026-05-21
**Status:** Approved

## Problem

The cyno-conservation gate was redefined in commit `ba29865` (2026-03-10,
"Refactor to whole-patch evaluation"). It replaced the original **local
sliding-window** test with a permissive **whole-patch percentage** test:

- **Old (≤ v1.0):** `MAX_CYNO_MISMATCHES_PER_600A2 = 2` — scan every ~600 Å²
  neighborhood (radius ≈ 13.8 Å, ~one VHH footprint) in a patch; if any
  neighborhood exceeds the limit, trim the mismatch-dense residues and
  re-cluster survivors into sub-patches ≥600 Å² (with a contiguity check).
- **Current:** `MAX_CYNO_MISMATCH_PERCENT = 15%`, size-scaled to a 30% cap —
  a whole-patch *average* mismatch tolerance.

The whole-patch mode is much more permissive for low-conservation targets. For
example, CEACAM5 (overall ~80% cyno identity, mismatches spread across the
patch) now yields epitopes that the old local test rejected. We want the
faithful old local behavior back **as the default**, while keeping the
whole-patch mode available and making both tunable.

## Decisions (from brainstorming)

1. **Fidelity:** full old behavior — local-density gate + trim-and-salvage
   (re-cluster survivors into sub-patches) + conserved sub-patch contiguity
   check.
2. **Whole-patch mode:** kept as an opt-in alternative, not removed.
3. **Exposure:** CLI flags + Python params on all three entry points
   (`run.py`, `single.py`, `bispecific.py`).

**Explicit non-goal:** the specificity (paralog) filter also has sliding-window
history; it stays on its current whole-patch percentage logic. This change
touches cyno conservation only.

## Changes

### config.py

```python
# Cyno conservation mode: "local" (sliding-window, default) | "whole_patch"
CYNO_MODE = "local"

# Local mode tunable: max cyno mismatches per ~600 Å² neighborhood
MAX_CYNO_MISMATCHES_PER_600A2 = 2        # (un-deprecated; local-mode threshold)

# Whole-patch mode tunable (unchanged)
MAX_CYNO_MISMATCH_PERCENT = 15.0
```

`MAX_CYNO_MISMATCHES_PER_600A2` loses its "DEPRECATED" comment and becomes the
documented local-mode knob.

### conservation.py

Restore the four helper functions **verbatim** from `ba29865~1:conservation.py`,
with one change: update the import path inside `_recluster_survivors` from
`from epitope_pipeline.surface import cluster_surface_patches, SurfacePatch` to
`from epitope_pipeline.compute.surface import cluster_surface_patches,
SurfacePatch`. The restored functions:

- `check_local_mismatch_density(patch_residues, residue_conservation, ca_coords)`
  → `(worst_resnum, worst_count)`. Uses `footprint_radius = sqrt(VHH_FOOTPRINT_MIN_A2/π)`.
- `identify_failing_residues(patch_residues, boolean_map, ca_coords, max_bad_per_window)`
  → `set` of residues to trim.
- `_recluster_survivors(survivors, ca_coords, residue_sasa, original_patch,
  residue_conservation, next_patch_id)` → `list[SurfacePatch]` (applies the
  Gate-2 contiguity check per sub-cluster; builds new `SurfacePatch` objects).
- `_check_conserved_subpatch(conserved_residues, ca_coords, residue_sasa)` →
  `bool` (largest conserved sub-cluster ≥ `VHH_FOOTPRINT_MIN_A2`).

The current `SurfacePatch` dataclass fields
(`patch_id, residue_numbers, residue_aas, total_sasa_a2, centroid,
max_dimension_a, avg_distance_from_membrane`) match exactly what
`_recluster_survivors` constructs — no field changes needed.

`analyze_conservation` gains mode dispatch:

```python
def analyze_conservation(target, surface_analysis, ca_coords=None,
                         cyno_mode=None, max_mismatches_per_window=None):
    ...
    mode = cyno_mode if cyno_mode is not None else config.CYNO_MODE
    max_bad = (max_mismatches_per_window if max_mismatches_per_window is not None
               else config.MAX_CYNO_MISMATCHES_PER_600A2)
    # (after alignment + per-residue conservation, for each patch:)
    #   mode == "local"       -> two-gate trim-and-salvage path (old logic),
    #                            using max_bad as the per-window threshold
    #   mode == "whole_patch" -> evaluate_patch_conservation (current logic)
```

- The `"local"` path is the old per-patch loop body (Gate 1 sliding window →
  trim → `_recluster_survivors`; Gate 2 contiguity), tracking `next_patch_id`
  from `max(existing patch_ids)+1`.
- Faithful quirk preserved: when `ca_coords` is falsy, the local path cannot run
  the spatial tests and the patch **passes** (matches old behavior).
- `evaluate_patch_conservation` is retained unchanged for `"whole_patch"`.

### Threading + CLI

`run_pipeline` (run.py), `run_single` (single.py), `run_bispecific`
(bispecific.py) each gain `cyno_mode=None` and `cyno_max_window_mismatches=None`
kwargs:
- Applied as config overrides at the top of the run (same pattern as
  `cyno_mismatch_percent`): if provided, set `config.CYNO_MODE` /
  `config.MAX_CYNO_MISMATCHES_PER_600A2`.
- Passed into every `analyze_conservation` call (run.py Step 6; bispecific
  `_analyze_target_zone`).

CLI flags on all three entry points (following the existing override-flag
pattern: argparse `default=None`, config supplies the real default, the flag
only overrides config when the user passes it):
- `--cyno-mode {local,whole-patch}` (argparse default `None`; help text states
  "default: local"; maps `whole-patch` → `"whole_patch"`, `local` → `"local"`).
- `--cyno-max-window-mismatches N` (type int, argparse default `None`; help text
  states "default: 2") — local-mode threshold.

Interactions (documented in `--help`):
- `--no-cyno` / `skip_cyno_gate` still bypasses the gate entirely (unchanged).
- `--cyno-mismatch-percent` only affects `whole_patch` mode.

### Run parameters block (honesty)

The `parameters` dict in run.py / bispecific.py records the active cyno mode and
its threshold:
- always: `"cyno_mode": <mode>`
- local mode: `"cyno_max_mismatches_per_600a2": <N>`
- whole_patch mode: keep existing `cyno_mismatch_percent_base` /
  `cyno_mismatch_scaling`.

## Testing

New `tests/test_compute_conservation_local.py` (Layer-1 unit tests, no network):
- `check_local_mismatch_density`: synthetic patch with a cluster of mismatches
  → worst_count reflects the dense neighborhood; an all-conserved patch → `(0,0)`.
- `identify_failing_residues`: returns the residues whose window exceeds the
  threshold; empty when none exceed.
- `_check_conserved_subpatch`: contiguous conserved residues with area ≥ min →
  True; sparse/insufficient area → False.
- `_recluster_survivors`: survivors forming a ≥600 Å² contiguous cluster yield a
  valid `SurfacePatch` with correct fields; sub-min clusters are dropped.
- Mode dispatch via `analyze_conservation` on a small synthetic target: a
  mismatch-dense patch is **rejected in `local` mode** and **passes in
  `whole_patch` mode**. (Construct minimal `TargetInfo`/`SurfaceAnalysis`
  fixtures with `ca_coords`; reuse existing fixtures in `tests/_fixtures.py`
  where possible.)

CLI tests (extend `tests/test_run_cli.py`, `tests/test_single_cli.py`, and add
`tests/test_bispecific_cli.py` if none exists):
- `--cyno-mode whole-patch` → forwards `cyno_mode="whole_patch"`.
- `--cyno-mode local` → forwards `cyno_mode="local"`.
- no flag → forwards `cyno_mode=None` (run function resolves to config default
  `local`).
- `--cyno-max-window-mismatches 3` → forwards `cyno_max_window_mismatches=3`;
  no flag → `None`.

## Docs

Update README cyno sections (the "Filtering Strategy" / "Default Thresholds" /
"Tuning Guidance" areas, and the single-protein section) to:
- describe the two cyno modes, with `local` as the default,
- document `--cyno-mode` and `--cyno-max-window-mismatches`,
- note that `--cyno-mismatch-percent` applies to `whole_patch` mode only.

## Verification

- `pytest tests/` green (new + existing).
- `epitope-single --help`, `epitope-pipeline --help`, `epitope-bispecific --help`
  all list `--cyno-mode` and `--cyno-max-window-mismatches`.
- Re-running the PFE panel under the default (`local`) yields stricter cyno
  filtering than the prior whole-patch default (CEACAM5 epitopes drop unless a
  clean conserved sub-patch survives).
