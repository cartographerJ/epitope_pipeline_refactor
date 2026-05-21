# Tunable Cyno Conservation Mode (local default) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restore the faithful old local sliding-window cyno-conservation test (with trim-and-salvage) as the default, keep the whole-patch percentage mode as an opt-in, and make both selectable/tunable via Python params and CLI flags on all three entry points.

**Architecture:** Add a `CYNO_MODE` config default (`"local"`). Restore four helper functions in `conservation.py` verbatim from commit `ba29865~1` (one import-path fix) and make `analyze_conservation` dispatch on mode. Thread selection through the entry points using the **existing config-override pattern** (like `cyno_mismatch_percent` / `cyno_max_mismatches` already do): the new Python param is `cyno_mode`; the per-window threshold reuses the existing `cyno_max_mismatches` param. CLI flags `--cyno-mode` and `--cyno-max-window-mismatches` map onto these.

**Tech Stack:** Python 3.10+, numpy, scipy, BioPython, pytest. Spec: `docs/superpowers/specs/2026-05-21-cyno-local-mode-design.md`.

**Test command:** use the project venv: `.venv/bin/pytest`. Run the suite scoped to `tests/` (`.venv/bin/pytest tests/`) — the `archive/` directory has pre-existing collection errors that are out of scope.

---

### Task 1: config — add CYNO_MODE, un-deprecate the local threshold

**Files:**
- Modify: `epitope_pipeline/config.py` (the conservation block ~lines 82-89)
- Test: `tests/test_config_cyno.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_config_cyno.py
"""Cyno-mode config defaults."""

from epitope_pipeline import config


def test_cyno_mode_defaults_to_local():
    assert config.CYNO_MODE == "local"


def test_local_threshold_default():
    assert config.MAX_CYNO_MISMATCHES_PER_600A2 == 2


def test_whole_patch_percent_default_kept():
    assert config.MAX_CYNO_MISMATCH_PERCENT == 15.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/pytest tests/test_config_cyno.py -v`
Expected: FAIL with `AttributeError: module 'epitope_pipeline.config' has no attribute 'CYNO_MODE'`

- [ ] **Step 3: Edit config**

In `epitope_pipeline/config.py`, replace this block:

```python
# ---------------------------------------------------------------------------
# Conservation
# ---------------------------------------------------------------------------
# Whole-patch evaluation (NEW)
MAX_CYNO_MISMATCH_PERCENT = 15.0        # Max % cyno mismatches to accept patch

# DEPRECATED: Legacy sliding window threshold (not used in whole-patch mode)
MAX_CYNO_MISMATCHES_PER_600A2 = 2       # Kept for backward compatibility only
```

with:

```python
# ---------------------------------------------------------------------------
# Conservation
# ---------------------------------------------------------------------------
# Cyno-conservation mode selector:
#   "local"       — sliding-window local-density test with trim-and-salvage
#                   (default; faithful to the pre-2026-03 behavior)
#   "whole_patch" — whole-patch average mismatch percentage (size-scaled)
CYNO_MODE = "local"

# Local mode: max cyno mismatches allowed in any ~600 Å² neighborhood of a patch
MAX_CYNO_MISMATCHES_PER_600A2 = 2

# Whole-patch mode: max % cyno mismatches to accept a patch (size-scaled to 30%)
MAX_CYNO_MISMATCH_PERCENT = 15.0
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/pytest tests/test_config_cyno.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add epitope_pipeline/config.py tests/test_config_cyno.py
git commit -m "Add CYNO_MODE config (default local); un-deprecate local threshold

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: conservation — restore the four local-mode helper functions

Restore `check_local_mismatch_density`, `identify_failing_residues`,
`_recluster_survivors`, `_check_conserved_subpatch` verbatim from
`ba29865~1:conservation.py`, with the import path inside `_recluster_survivors`
updated to the current package layout (`epitope_pipeline.compute.surface`).

**Files:**
- Modify: `epitope_pipeline/compute/conservation.py` (append the four functions after `evaluate_patch_conservation`)
- Test: `tests/test_compute_conservation_local.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_compute_conservation_local.py
"""Unit tests for the restored local sliding-window cyno helpers."""

import numpy as np

from epitope_pipeline.compute.conservation import (
    check_local_mismatch_density,
    identify_failing_residues,
    _check_conserved_subpatch,
    _recluster_survivors,
)
from tests._fixtures import make_patch


def test_density_counts_clustered_mismatches():
    # 5 residues all within ~14 Å of each other; residues 1,2,3 are mismatches
    coords = {1: (0, 0, 0), 2: (3, 0, 0), 3: (6, 0, 0), 4: (9, 0, 0), 5: (12, 0, 0)}
    coords = {r: np.array(c, dtype=float) for r, c in coords.items()}
    cons = {1: False, 2: False, 3: False, 4: True, 5: True}
    worst_res, worst_count = check_local_mismatch_density([1, 2, 3, 4, 5], cons, coords)
    assert worst_count == 3


def test_density_all_conserved_is_zero():
    coords = {r: np.array((r * 3.0, 0, 0)) for r in range(1, 6)}
    cons = {r: True for r in range(1, 6)}
    assert check_local_mismatch_density([1, 2, 3, 4, 5], cons, coords) == (0, 0)


def test_identify_failing_residues_local_cluster():
    # Cluster A (1-3) mismatches close together; cluster B (10-12) conserved + far
    coords = {
        1: (0, 0, 0), 2: (2, 0, 0), 3: (4, 0, 0),
        10: (100, 0, 0), 11: (102, 0, 0), 12: (104, 0, 0),
    }
    coords = {r: np.array(c, dtype=float) for r, c in coords.items()}
    cons = {1: False, 2: False, 3: False, 10: True, 11: True, 12: True}
    failing = identify_failing_residues([1, 2, 3, 10, 11, 12], cons, coords, 2)
    assert failing == {1, 2, 3}


def test_identify_failing_none_when_threshold_high():
    coords = {1: (0, 0, 0), 2: (2, 0, 0), 3: (4, 0, 0)}
    coords = {r: np.array(c, dtype=float) for r, c in coords.items()}
    cons = {1: False, 2: False, 3: False}
    assert identify_failing_residues([1, 2, 3], cons, coords, 5) == set()


def test_conserved_subpatch_contiguous_passes():
    # 7 conserved residues, 5 Å apart (contiguous), 100 Å² each = 700 ≥ 600
    coords = {r: np.array((i * 5.0, 0, 0)) for i, r in enumerate(range(1, 8))}
    sasa = {r: 100.0 for r in range(1, 8)}
    assert _check_conserved_subpatch(list(range(1, 8)), coords, sasa) is True


def test_conserved_subpatch_sparse_fails():
    # 2 residues far apart, 100 Å² each → no single cluster ≥ 600
    coords = {1: np.array((0.0, 0, 0)), 2: np.array((100.0, 0, 0))}
    sasa = {1: 100.0, 2: 100.0}
    assert _check_conserved_subpatch([1, 2], coords, sasa) is False


def test_recluster_survivors_builds_subpatch():
    # 7 contiguous conserved survivors, 100 Å² each = 700 ≥ 600 → one sub-patch
    survivors = list(range(1, 8))
    coords = {r: np.array((i * 5.0, 0, 0)) for i, r in enumerate(survivors)}
    sasa = {r: 100.0 for r in survivors}
    cons = {r: True for r in survivors}
    original = make_patch(patch_id=0, residue_numbers=survivors,
                          residue_aas={r: "A" for r in survivors})
    subs = _recluster_survivors(survivors, coords, sasa, original, cons, next_patch_id=1)
    assert len(subs) == 1
    assert subs[0].patch_id == 1
    assert sorted(subs[0].residue_numbers) == survivors
    assert subs[0].total_sasa_a2 == 700.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/pytest tests/test_compute_conservation_local.py -v`
Expected: FAIL with `ImportError: cannot import name 'check_local_mismatch_density'`

- [ ] **Step 3: Append the four helpers to conservation.py**

Append to `epitope_pipeline/compute/conservation.py` (after `evaluate_patch_conservation`, end of file). Note the import path in `_recluster_survivors` is `epitope_pipeline.compute.surface` (changed from the old `epitope_pipeline.surface`):

```python
# ---------------------------------------------------------------------------
# Local sliding-window evaluation (faithful to pre-2026-03 behavior)
# ---------------------------------------------------------------------------

def check_local_mismatch_density(patch_residues, residue_conservation, ca_coords):
    """
    Sliding window: for every residue in the patch, examine its ~600A²
    neighborhood (radius = sqrt(600/π) ≈ 13.8A) and count mismatches.

    Returns (worst_resnum, worst_mismatch_count) — the residue whose
    neighborhood has the most mismatches. (0, 0) if no coords available.
    """
    from scipy.spatial.distance import cdist

    footprint_radius = math.sqrt(VHH_FOOTPRINT_MIN_A2 / math.pi)

    res_with_coords = [r for r in patch_residues if r in ca_coords]
    if not res_with_coords:
        return (0, 0)

    coords = np.array([ca_coords[r] for r in res_with_coords])
    is_bad = np.array([
        not residue_conservation.get(r, False) for r in res_with_coords
    ])

    dists = cdist(coords, coords)

    worst_pos = 0
    worst_count = 0
    for i, resnum in enumerate(res_with_coords):
        neighbors = dists[i] <= footprint_radius
        local_bad = int(np.sum(is_bad[neighbors]))
        if local_bad > worst_count:
            worst_count = local_bad
            worst_pos = resnum

    return (worst_pos, worst_count)


# Backward-compatible alias for the old private name
_check_local_mismatch_density = check_local_mismatch_density


def identify_failing_residues(patch_residues, boolean_map, ca_coords,
                              max_bad_per_window):
    """
    Identify all residues whose ~600A² neighborhood contains more than
    max_bad_per_window mismatches. These are the residues to trim.

    Returns a set of residue numbers that should be trimmed.
    """
    from scipy.spatial.distance import cdist

    footprint_radius = math.sqrt(VHH_FOOTPRINT_MIN_A2 / math.pi)

    res_with_coords = [r for r in patch_residues if r in ca_coords]
    if not res_with_coords:
        return set()

    coords = np.array([ca_coords[r] for r in res_with_coords])
    is_bad = np.array([
        not boolean_map.get(r, False) for r in res_with_coords
    ])

    dists = cdist(coords, coords)

    failing = set()
    for i, resnum in enumerate(res_with_coords):
        neighbors = dists[i] <= footprint_radius
        local_bad = int(np.sum(is_bad[neighbors]))
        if local_bad > max_bad_per_window:
            failing.add(resnum)

    return failing


def _check_conserved_subpatch(conserved_residues, ca_coords, residue_sasa):
    """
    Verify that conserved residues within a patch form a contiguous 3D region
    large enough for a VHH CDR footprint (largest sub-cluster >= VHH min).
    """
    if len(conserved_residues) < 2:
        area = sum(residue_sasa.get(r, 0.0) for r in conserved_residues)
        return area >= VHH_FOOTPRINT_MIN_A2

    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import fcluster, linkage

    res_list = sorted(conserved_residues)
    coords = np.array([ca_coords[r] for r in res_list])

    dist_matrix = pdist(coords)
    Z = linkage(dist_matrix, method="single")
    labels = fcluster(Z, t=PATCH_CLUSTERING_DISTANCE_A, criterion="distance")

    clusters = {}
    for res, label in zip(res_list, labels):
        clusters.setdefault(label, []).append(res)

    max_area = 0.0
    for cluster_residues in clusters.values():
        area = sum(residue_sasa.get(r, 0.0) for r in cluster_residues)
        if area > max_area:
            max_area = area

    return max_area >= VHH_FOOTPRINT_MIN_A2


def _recluster_survivors(survivors, ca_coords, residue_sasa, original_patch,
                         residue_conservation, next_patch_id):
    """
    Re-cluster survivor residues after trimming and build new SurfacePatch
    objects for sub-clusters >= VHH min that also pass the conserved sub-patch
    contiguity check.
    """
    from epitope_pipeline.compute.surface import cluster_surface_patches, SurfacePatch

    survivors_with_ca = [r for r in survivors if r in ca_coords]
    if len(survivors_with_ca) < 2:
        return []

    clusters = cluster_surface_patches(survivors_with_ca, ca_coords)

    result = []
    pid = next_patch_id
    for cluster_residues in clusters:
        total_area = sum(residue_sasa.get(r, 0.0) for r in cluster_residues)
        if total_area < VHH_FOOTPRINT_MIN_A2:
            continue

        conserved_in_cluster = [
            r for r in cluster_residues
            if residue_conservation.get(r, False) and r in ca_coords
        ]
        if not _check_conserved_subpatch(conserved_in_cluster, ca_coords, residue_sasa):
            logger.info(
                "      Sub-patch (%d residues, %.0f A²): conserved residues "
                "don't form contiguous >= %.0f A², skipping",
                len(cluster_residues), total_area, VHH_FOOTPRINT_MIN_A2,
            )
            continue

        coords = np.array([ca_coords[r] for r in cluster_residues])
        centroid = np.mean(coords, axis=0)

        if len(coords) > 1:
            from scipy.spatial.distance import pdist
            max_dim = float(np.max(pdist(coords)))
        else:
            max_dim = 0.0

        sp = SurfacePatch(
            patch_id=pid,
            residue_numbers=sorted(cluster_residues),
            residue_aas={r: original_patch.residue_aas.get(r, "X")
                         for r in cluster_residues},
            total_sasa_a2=total_area,
            centroid=centroid,
            max_dimension_a=max_dim,
            avg_distance_from_membrane=original_patch.avg_distance_from_membrane,
        )
        result.append(sp)
        pid += 1

    return result
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/pytest tests/test_compute_conservation_local.py -v`
Expected: PASS (7 tests)

- [ ] **Step 5: Commit**

```bash
git add epitope_pipeline/compute/conservation.py tests/test_compute_conservation_local.py
git commit -m "Restore local sliding-window cyno helpers (verbatim + import fix)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: conservation — mode dispatch in analyze_conservation

Add `cyno_mode` / `max_mismatches_per_window` params and branch the per-patch
evaluation between the local trim-and-salvage path and the existing whole-patch
path.

**Files:**
- Modify: `epitope_pipeline/compute/conservation.py` (`analyze_conservation`, ~lines 54-163)
- Test: `tests/test_compute_conservation_dispatch.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_compute_conservation_dispatch.py
"""analyze_conservation dispatches between local and whole_patch modes."""

from types import SimpleNamespace

import numpy as np

from epitope_pipeline.compute.conservation import analyze_conservation
from tests._fixtures import make_patch


def _build_target_and_surface():
    """30-residue patch: residues 1-4 are clustered mismatches, 5-30 conserved.

    Sequences are equal length and differ only at positions 1-4 (1:1 alignment),
    so residue_conservation maps cleanly: 1-4 False, 5-30 True.
    """
    human = "A" * 30
    cyno = "C" * 4 + "A" * 26  # positions 1-4 differ
    target = SimpleNamespace(
        uniprot_id="TEST", gene_name="TEST",
        sequence=human, cyno_sequence=cyno,
    )

    residues = list(range(1, 31))
    patch = make_patch(patch_id=0, residue_numbers=residues,
                       residue_aas={r: "A" for r in residues})
    surface = SimpleNamespace(
        patches=[patch],
        residue_sasa={r: 30.0 for r in residues},  # 30 * 30 = 900 Å² total
    )

    # Coords: residues 1-4 tightly clustered at the origin; 5-30 a far, contiguous
    # chain spaced 5 Å (so they re-cluster into one sub-patch).
    ca = {1: (0, 0, 0), 2: (2, 0, 0), 3: (4, 0, 0), 4: (6, 0, 0)}
    for i, r in enumerate(range(5, 31)):
        ca[r] = (100.0 + i * 5.0, 0, 0)
    ca = {r: np.array(c, dtype=float) for r, c in ca.items()}
    return target, surface, ca


def test_whole_patch_keeps_full_patch():
    target, surface, ca = _build_target_and_surface()
    res = analyze_conservation(target, surface, ca, cyno_mode="whole_patch")
    # 30 residues, 4 mismatches = 13.3% <= scaled 18.4% threshold → full patch kept
    assert len(res.conserved_patches) == 1
    assert len(res.conserved_patches[0].residue_numbers) == 30


def test_local_trims_mismatch_dense_residues():
    target, surface, ca = _build_target_and_surface()
    res = analyze_conservation(target, surface, ca, cyno_mode="local",
                               max_mismatches_per_window=2)
    # Local gate trims clustered mismatches 1-4, salvages residues 5-30
    assert len(res.conserved_patches) == 1
    assert set(res.conserved_patches[0].residue_numbers) == set(range(5, 31))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/pytest tests/test_compute_conservation_dispatch.py -v`
Expected: FAIL — `test_local_trims_mismatch_dense_residues` fails (current code only does whole-patch, so it keeps all 30) and/or `analyze_conservation` rejects the unexpected `cyno_mode` kwarg with `TypeError`.

- [ ] **Step 3: Add mode dispatch**

In `epitope_pipeline/compute/conservation.py`, change the `analyze_conservation`
signature from:

```python
def analyze_conservation(target, surface_analysis, ca_coords=None):
```

to:

```python
def analyze_conservation(target, surface_analysis, ca_coords=None,
                         cyno_mode=None, max_mismatches_per_window=None):
```

Then replace the per-patch evaluation block — everything from the comment
`# Step 3: Evaluate each patch (whole-patch mode, no trimming)` down to (but not
including) the final `logger.info("  %s conservation: %d patches pass...` summary
line — with this dispatch:

```python
    # Resolve mode + threshold (params override config defaults)
    mode = cyno_mode if cyno_mode is not None else config.CYNO_MODE
    max_bad = (max_mismatches_per_window if max_mismatches_per_window is not None
               else config.MAX_CYNO_MISMATCHES_PER_600A2)

    conserved_patches = []
    rejected_patches = []
    patch_conservation = {}

    if mode == "local":
        next_patch_id = max(p.patch_id for p in surface_analysis.patches) + 1
        for patch in surface_analysis.patches:
            n_residues = len(patch.residue_numbers)
            n_conserved = sum(
                1 for r in patch.residue_numbers
                if residue_conservation.get(r, False)
            )
            n_mismatches = n_residues - n_conserved
            identity = n_conserved / n_residues if n_residues > 0 else 0.0
            patch_conservation[patch.patch_id] = identity
            logger.info(
                "    Patch %d: %d residues (%.0f A²), cyno identity %.1f%% "
                "(%d/%d conserved, %d mismatches)",
                patch.patch_id, n_residues, patch.total_sasa_a2,
                identity * 100, n_conserved, n_residues, n_mismatches,
            )

            # Gate 1: sliding window — trim mismatch-dense residues, salvage survivors
            if ca_coords and n_mismatches > 0:
                worst_pos, worst_count = check_local_mismatch_density(
                    patch.residue_numbers, residue_conservation, ca_coords,
                )
                if worst_count > max_bad:
                    failing = identify_failing_residues(
                        patch.residue_numbers, residue_conservation, ca_coords, max_bad,
                    )
                    survivors = [r for r in patch.residue_numbers if r not in failing]
                    logger.info(
                        "      Sliding window: %d residues in mismatch-dense zones "
                        "(worst: res %d with %d), trimming...",
                        len(failing), worst_pos, worst_count,
                    )
                    sub_patches = _recluster_survivors(
                        survivors, ca_coords, surface_analysis.residue_sasa,
                        patch, residue_conservation, next_patch_id,
                    )
                    if sub_patches:
                        for sp in sub_patches:
                            sp_n = len(sp.residue_numbers)
                            sp_conserved = sum(
                                1 for r in sp.residue_numbers
                                if residue_conservation.get(r, False)
                            )
                            sp_identity = sp_conserved / sp_n if sp_n > 0 else 0.0
                            patch_conservation[sp.patch_id] = sp_identity
                            conserved_patches.append(sp)
                            logger.info(
                                "      -> Sub-patch %d: %d residues, %.0f A², "
                                "cyno %.1f%% -> PASSED",
                                sp.patch_id, sp_n, sp.total_sasa_a2, sp_identity * 100,
                            )
                            next_patch_id = sp.patch_id + 1
                    else:
                        rejected_patches.append((patch, identity))
                        logger.info(
                            "      -> REJECTED (no surviving sub-patches >= %.0f A²)",
                            VHH_FOOTPRINT_MIN_A2,
                        )
                    continue
                logger.info(
                    "      Sliding window OK (worst: %d mismatches near res %d)",
                    worst_count, worst_pos,
                )

            # Gate 2: conserved residues must form a contiguous >= VHH-min patch
            if ca_coords:
                conserved_residues = [
                    r for r in patch.residue_numbers
                    if residue_conservation.get(r, False) and r in ca_coords
                ]
                if not _check_conserved_subpatch(
                        conserved_residues, ca_coords, surface_analysis.residue_sasa):
                    rejected_patches.append((patch, identity))
                    logger.info(
                        "      -> REJECTED (conserved residues don't form contiguous "
                        "patch >= %.0f A²)", VHH_FOOTPRINT_MIN_A2,
                    )
                    continue

            conserved_patches.append(patch)
            logger.info("      -> PASSED")
    else:
        # Whole-patch mode: average mismatch percentage (threshold scales with size)
        for patch in surface_analysis.patches:
            n_residues = len(patch.residue_numbers)
            n_conserved = sum(
                1 for r in patch.residue_numbers
                if residue_conservation.get(r, False)
            )
            n_mismatches = n_residues - n_conserved
            identity_fraction = n_conserved / n_residues if n_residues > 0 else 0.0
            mismatch_percent = (n_mismatches / n_residues) * 100.0
            patch_conservation[patch.patch_id] = identity_fraction
            logger.info(
                "    Patch %d: %d residues (%.0f A²), cyno identity %.1f%% "
                "(%d conserved, %d mismatches = %.1f%%)",
                patch.patch_id, n_residues, patch.total_sasa_a2,
                identity_fraction * 100, n_conserved, n_mismatches, mismatch_percent,
            )
            passes, _, effective_thresh = evaluate_patch_conservation(
                patch, residue_conservation)
            if passes:
                conserved_patches.append(patch)
                logger.info("      -> PASSED (%.1f%% mismatches <= %.1f%% threshold, %d residues)",
                            mismatch_percent, effective_thresh, n_residues)
            else:
                rejected_patches.append((patch, identity_fraction))
                logger.info("      -> REJECTED (%.1f%% mismatches > %.1f%% threshold, %d residues)",
                            mismatch_percent, effective_thresh, n_residues)
```

(Leave the existing `logger.info("  %s conservation: %d patches pass, %d rejected", ...)`
summary and the `return ConservationResult(...)` block unchanged — they already
read `conserved_patches`, `rejected_patches`, `patch_conservation`.)

Also update the module docstring and `analyze_conservation` docstring step list
to mention both modes (the local sliding-window default and the whole-patch
alternative).

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/pytest tests/test_compute_conservation_dispatch.py tests/test_compute_conservation.py tests/test_compute_conservation_local.py -v`
Expected: PASS (dispatch 2 + existing whole-patch tests + local helper 7)

- [ ] **Step 5: Commit**

```bash
git add epitope_pipeline/compute/conservation.py tests/test_compute_conservation_dispatch.py
git commit -m "Dispatch analyze_conservation between local and whole_patch modes

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: run.py — params, config override, parameters block, CLI flags

Reuse the existing `cyno_max_mismatches` param for the per-window threshold; add
a new `cyno_mode` param. Wire CLI flags `--cyno-mode` and
`--cyno-max-window-mismatches`.

**Files:**
- Modify: `epitope_pipeline/run.py` (signature ~88-101; override block ~133-141; parameters dict ~180-189; argparser ~618-627; main `run_pipeline(...)` call ~664-672)
- Test: `tests/test_run_cli.py` (extend)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_run_cli.py`:

```python
def test_cyno_mode_default_is_none(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2"])
    assert captured[0]["cyno_mode"] is None


def test_cyno_mode_whole_patch(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--cyno-mode", "whole-patch"])
    assert captured[0]["cyno_mode"] == "whole_patch"


def test_cyno_mode_local(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--cyno-mode", "local"])
    assert captured[0]["cyno_mode"] == "local"


def test_cyno_max_window_mismatches(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--cyno-max-window-mismatches", "3"])
    assert captured[0]["cyno_max_mismatches"] == 3
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/pytest tests/test_run_cli.py -k cyno_mode -v`
Expected: FAIL — `KeyError: 'cyno_mode'` (stub didn't receive it) / argparse error on unknown `--cyno-mode`.

- [ ] **Step 3a: Add `cyno_mode` to the signature**

In `epitope_pipeline/run.py`, change the `run_pipeline` signature line:

```python
    cyno_max_mismatches=None,
    cyno_mismatch_percent=None,
```

to:

```python
    cyno_max_mismatches=None,
    cyno_mismatch_percent=None,
    cyno_mode=None,
```

- [ ] **Step 3b: Apply the config override**

In the `# --- Apply overrides ---` block, after the `cyno_max_mismatches` line:

```python
    if cyno_max_mismatches is not None:
        config.MAX_CYNO_MISMATCHES_PER_600A2 = cyno_max_mismatches
```

add:

```python
    if cyno_mode is not None:
        config.CYNO_MODE = cyno_mode
```

- [ ] **Step 3c: Record mode in the parameters dict**

In the `parameters = {` dict, change:

```python
        "cyno_mismatch_percent_base": config.MAX_CYNO_MISMATCH_PERCENT,
        "cyno_mismatch_scaling": "min(base% * sqrt(n_residues/20), 30%)",
```

to:

```python
        "cyno_mode": config.CYNO_MODE,
        "cyno_max_mismatches_per_600a2": config.MAX_CYNO_MISMATCHES_PER_600A2,
        "cyno_mismatch_percent_base": config.MAX_CYNO_MISMATCH_PERCENT,
        "cyno_mismatch_scaling": "min(base% * sqrt(n_residues/20), 30%)",
```

- [ ] **Step 3d: Add CLI flags**

In `build_argparser()`, after the `--cyno-mismatch-percent` argument:

```python
    parser.add_argument("--cyno-mode", choices=["local", "whole-patch"],
                        default=None, dest="cyno_mode",
                        help="Cyno conservation mode (default: local). "
                             "'local' = sliding-window local test; "
                             "'whole-patch' = average %% mismatch.")
    parser.add_argument("--cyno-max-window-mismatches", type=int, default=None,
                        dest="cyno_max_window_mismatches",
                        help="Local mode: max cyno mismatches per ~600 Å² window "
                             "(default: 2)")
```

- [ ] **Step 3e: Forward from main()**

In `main()`, change the `run_pipeline(` call to add (after `cyno_mismatch_percent=args.cyno_mismatch_percent,`):

```python
        cyno_mismatch_percent=args.cyno_mismatch_percent,
        cyno_mode=(args.cyno_mode.replace("-", "_") if args.cyno_mode else None),
        cyno_max_mismatches=args.cyno_max_window_mismatches,
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/pytest tests/test_run_cli.py -v`
Expected: PASS (existing + 4 new)

- [ ] **Step 5: Commit**

```bash
git add epitope_pipeline/run.py tests/test_run_cli.py
git commit -m "run.py: --cyno-mode / --cyno-max-window-mismatches + config wiring

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: single.py — pass cyno_mode + threshold through

**Files:**
- Modify: `epitope_pipeline/single.py` (run_single signature/docstring/call ~23-62; argparser ~85-90; main call ~123-132)
- Test: `tests/test_single_cli.py` (extend)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_single_cli.py`:

```python
def test_single_cyno_mode_whole_patch(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--cyno-mode", "whole-patch"])
    assert captured[0]["cyno_mode"] == "whole_patch"


def test_single_cyno_mode_default_none(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2"])
    assert captured[0]["cyno_mode"] is None


def test_single_cyno_max_window_mismatches(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--cyno-max-window-mismatches", "4"])
    assert captured[0]["cyno_max_mismatches"] == 4
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/pytest tests/test_single_cli.py -k cyno -v`
Expected: FAIL — argparse rejects `--cyno-mode` / stub missing keys.

- [ ] **Step 3a: Update run_single signature, docstring, and delegation**

In `epitope_pipeline/single.py`, change the `run_single` signature:

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
```

to:

```python
def run_single(
    identifiers,
    run_name=None,
    cyno_mismatch_percent=None,
    cyno_mode=None,
    cyno_max_mismatches=None,
    skip_cyno_gate=False,
    nonspecific_percent=None,
    force_experimental=False,
    aliases=None,
    verbose=True,
):
```

Replace the docstring's "except the deprecated cyno_max_mismatches count is
intentionally omitted — use cyno_mismatch_percent for cyno-conservation tuning."
sentence with:

```
    Cyno conservation runs in config.CYNO_MODE by default ("local"); pass
    cyno_mode="whole_patch" for the average-percentage test. cyno_max_mismatches
    is the local-mode per-~600 Å² window threshold; cyno_mismatch_percent tunes
    whole_patch mode.
```

In the `return run_pipeline(` call, after `cyno_mismatch_percent=cyno_mismatch_percent,` add:

```python
        cyno_mismatch_percent=cyno_mismatch_percent,
        cyno_mode=cyno_mode,
        cyno_max_mismatches=cyno_max_mismatches,
```

- [ ] **Step 3b: Add CLI flags**

In `build_argparser()`, after the `--cyno-mismatch-percent` argument:

```python
    parser.add_argument("--cyno-mode", choices=["local", "whole-patch"],
                        default=None, dest="cyno_mode",
                        help="Cyno conservation mode (default: local). "
                             "'local' = sliding-window local test; "
                             "'whole-patch' = average %% mismatch.")
    parser.add_argument("--cyno-max-window-mismatches", type=int, default=None,
                        dest="cyno_max_window_mismatches",
                        help="Local mode: max cyno mismatches per ~600 Å² window "
                             "(default: 2)")
```

- [ ] **Step 3c: Forward from main()**

In `main()`, in the `run_single(` call, after `cyno_mismatch_percent=args.cyno_mismatch_percent,` add:

```python
        cyno_mismatch_percent=args.cyno_mismatch_percent,
        cyno_mode=(args.cyno_mode.replace("-", "_") if args.cyno_mode else None),
        cyno_max_mismatches=args.cyno_max_window_mismatches,
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/pytest tests/test_single_cli.py tests/test_single.py -v`
Expected: PASS (existing + 3 new)

- [ ] **Step 5: Commit**

```bash
git add epitope_pipeline/single.py tests/test_single_cli.py
git commit -m "single.py: forward cyno_mode + window threshold

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: bispecific.py — params, config override, parameters block, CLI flags

**Files:**
- Modify: `epitope_pipeline/bispecific.py` (run_bispecific signature ~88-98; override block ~119-127; parameters dict ~358-368; main argparser + call ~640-664)
- Test: `tests/test_bispecific_cli.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_bispecific_cli.py
"""Argparse / CLI tests for epitope_pipeline.bispecific.main.

run_bispecific is monkeypatched to a recording stub so flags are exercised
without running the actual pipeline.
"""

import sys

import pytest

from epitope_pipeline import bispecific as bispecific_module


@pytest.fixture
def captured(monkeypatch):
    calls = []

    def stub(pairs, **kwargs):
        calls.append({"pairs": list(pairs), **kwargs})
        return {"run_dir": "/tmp/fake", "pair_results": [], "zone_results": {}}

    monkeypatch.setattr(bispecific_module, "run_bispecific", stub)
    return calls


def _invoke(monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", ["epitope-bispecific", *argv])
    bispecific_module.main()


def test_pair_parsing(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2:NECTIN4"])
    assert captured[0]["pairs"] == [("ERBB2", "NECTIN4")]


def test_cyno_mode_whole_patch(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2:NECTIN4", "--cyno-mode", "whole-patch"])
    assert captured[0]["cyno_mode"] == "whole_patch"


def test_cyno_mode_default_none(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2:NECTIN4"])
    assert captured[0]["cyno_mode"] is None


def test_cyno_max_window_mismatches(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2:NECTIN4", "--cyno-max-window-mismatches", "1"])
    assert captured[0]["cyno_max_mismatches"] == 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/pytest tests/test_bispecific_cli.py -v`
Expected: FAIL — argparse rejects `--cyno-mode` / stub missing keys.

- [ ] **Step 3a: Add `cyno_mode` to run_bispecific signature**

In `epitope_pipeline/bispecific.py`, change:

```python
    cyno_max_mismatches=None,
    cyno_mismatch_percent=None,
    nonspecific_percent=None,
```

to:

```python
    cyno_max_mismatches=None,
    cyno_mismatch_percent=None,
    cyno_mode=None,
    nonspecific_percent=None,
```

- [ ] **Step 3b: Apply config override**

After the existing legacy override:

```python
    # Apply legacy threshold overrides (deprecated)
    if cyno_max_mismatches is not None:
        config.MAX_CYNO_MISMATCHES_PER_600A2 = cyno_max_mismatches
```

add:

```python
    if cyno_mode is not None:
        config.CYNO_MODE = cyno_mode
```

- [ ] **Step 3c: Record mode in the parameters dict**

In the `parameters = {` dict, change:

```python
        "cyno_mismatch_percent_base": config.MAX_CYNO_MISMATCH_PERCENT,
        "cyno_mismatch_scaling": "min(base% * sqrt(n_residues/20), 30%)",
```

to:

```python
        "cyno_mode": config.CYNO_MODE,
        "cyno_max_mismatches_per_600a2": config.MAX_CYNO_MISMATCHES_PER_600A2,
        "cyno_mismatch_percent_base": config.MAX_CYNO_MISMATCH_PERCENT,
        "cyno_mismatch_scaling": "min(base% * sqrt(n_residues/20), 30%)",
```

- [ ] **Step 3d: Add CLI flags + forward from main()**

In `main()`'s argparser, after the `--proximal` argument:

```python
    parser.add_argument("--cyno-mode", choices=["local", "whole-patch"],
                        default=None, dest="cyno_mode",
                        help="Cyno conservation mode (default: local). "
                             "'local' = sliding-window local test; "
                             "'whole-patch' = average %% mismatch.")
    parser.add_argument("--cyno-max-window-mismatches", type=int, default=None,
                        dest="cyno_max_window_mismatches",
                        help="Local mode: max cyno mismatches per ~600 Å² window "
                             "(default: 2)")
```

Change the `run_bispecific(` call:

```python
    results = run_bispecific(
        pairs,
        distal_min_distance_a=args.distal,
        proximal_max_distance_a=args.proximal,
    )
```

to:

```python
    results = run_bispecific(
        pairs,
        distal_min_distance_a=args.distal,
        proximal_max_distance_a=args.proximal,
        cyno_mode=(args.cyno_mode.replace("-", "_") if args.cyno_mode else None),
        cyno_max_mismatches=args.cyno_max_window_mismatches,
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/pytest tests/test_bispecific_cli.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add epitope_pipeline/bispecific.py tests/test_bispecific_cli.py
git commit -m "bispecific.py: --cyno-mode / --cyno-max-window-mismatches + wiring

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: README docs + full-suite verification

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Run the full suite**

Run: `.venv/bin/pytest tests/ -q`
Expected: all pass (existing + Tasks 1-6 additions).

- [ ] **Step 2: Verify all three CLIs expose the flags**

Run: `.venv/bin/epitope-single --help` then `.venv/bin/epitope-pipeline --help` then `.venv/bin/epitope-bispecific --help`
Expected: each lists `--cyno-mode` and `--cyno-max-window-mismatches`.

- [ ] **Step 3: Update README**

In `README.md`, in the conservation/tuning area (the "Default Thresholds" /
"Tuning Guidance" sections around the cyno discussion), add a "Cyno conservation
modes" subsection. Insert this block (adjust heading level to match neighbors):

```markdown
### Cyno conservation modes

Cyno-conservation filtering runs in one of two modes (default **local**):

- **`local`** (default) — a sliding-window local-density test. Every ~600 Å²
  neighborhood (~one VHH footprint) of a patch may contain at most
  `--cyno-max-window-mismatches` cyno mismatches (default **2**); residues in
  mismatch-dense zones are trimmed and the survivors re-clustered into
  sub-patches ≥600 Å². This is strict about *local* conservation and is the
  faithful pre-2026-03 behavior.
- **`whole-patch`** — a whole-patch *average* mismatch percentage, size-scaled:
  `min(--cyno-mismatch-percent × √(n/20), 30%)`. Permissive for large patches.

```bash
# Default (local), tighten the per-window limit
epitope-single ERBB2 --cyno-max-window-mismatches 1

# Opt into the permissive whole-patch average instead
epitope-single ERBB2 --cyno-mode whole-patch --cyno-mismatch-percent 20
```

`--cyno-mismatch-percent` applies only to `whole-patch` mode;
`--cyno-max-window-mismatches` applies only to `local` mode. `--no-cyno` still
bypasses the gate entirely. All three entry points (`epitope-pipeline`,
`epitope-single`, `epitope-bispecific`) accept these flags.
```

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "Document local vs whole-patch cyno conservation modes

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Final verification

- [ ] `.venv/bin/pytest tests/ -q` — all green.
- [ ] `epitope-single ERBB2 --help` lists `--cyno-mode` and `--cyno-max-window-mismatches`.
- [ ] (Optional, network) Re-run a known case under the new default and confirm `parameters` records `"cyno_mode": "local"`:
  `epitope-single CEACAM5=PFE-003 --run-name cyno_local_check` → inspect the run's `summary_metrics.csv` / params for `cyno_mode = local`, and confirm stricter cyno filtering than the prior whole-patch default.
