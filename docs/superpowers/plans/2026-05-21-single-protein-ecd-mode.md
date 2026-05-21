# Single-protein ECD Suitability Mode Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a clean, first-class `epitope-single` entry point that assesses whether a single protein's whole extracellular domain is a suitable VHH-epitope target, with no distal/proximal zone framing.

**Architecture:** Thin wrapper over the existing per-target `run_pipeline`, calling it with `no_distance_filter=True`. No orchestration loop is duplicated. Three small "honesty" fixes correct the mode/label text in `run.py` and `visualize.py` that is currently wrong when the distance filter is off; each fix is extracted into a pure helper so it can be unit-tested without running the network/BLAST pipeline.

**Tech Stack:** Python 3.10+, argparse, pytest. Spec: `docs/superpowers/specs/2026-05-21-single-protein-ecd-mode-design.md`.

---

### Task 1: Mode-label helper in run.py

When `no_distance_filter=True`, `parameters["mode"]` currently reads `"distal (min_distance_a=80)"` because the ternary only checks `max_distance_a`. Extract a pure helper and use it.

**Files:**
- Modify: `epitope_pipeline/run.py` (add helper near top of module; replace the `parameters["mode"]` expression ~line 147)
- Test: `tests/test_run_mode_labels.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_run_mode_labels.py
"""Unit tests for the pure mode/label helpers in epitope_pipeline.run."""

from epitope_pipeline.run import _mode_label
from epitope_pipeline import config


def test_mode_label_whole_ecd_takes_priority():
    # no_distance_filter wins even if a max_distance is somehow present
    assert _mode_label(no_distance_filter=True, max_distance_a=50.0) == (
        "whole ectodomain (no distance filter)"
    )
    assert _mode_label(no_distance_filter=True, max_distance_a=None) == (
        "whole ectodomain (no distance filter)"
    )


def test_mode_label_proximal():
    assert _mode_label(no_distance_filter=False, max_distance_a=50.0) == (
        "proximal (max_distance_a=50.0)"
    )


def test_mode_label_distal_default():
    assert _mode_label(no_distance_filter=False, max_distance_a=None) == (
        "distal (min_distance_a={})".format(config.ECTODOMAIN_MIN_DISTANCE_A)
    )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_run_mode_labels.py -v`
Expected: FAIL with `ImportError: cannot import name '_mode_label'`

- [ ] **Step 3: Add the helper**

Add this function in `epitope_pipeline/run.py` immediately above `def run_pipeline(` (after the `logger = logging.getLogger(...)` line and section comment):

```python
def _mode_label(no_distance_filter, max_distance_a):
    """Human-readable spatial-mode label for the run parameters block."""
    if no_distance_filter:
        return "whole ectodomain (no distance filter)"
    if max_distance_a:
        return "proximal (max_distance_a={})".format(max_distance_a)
    return "distal (min_distance_a={})".format(config.ECTODOMAIN_MIN_DISTANCE_A)
```

- [ ] **Step 4: Use the helper in run_pipeline**

In `epitope_pipeline/run.py`, replace this line inside the `parameters = {` dict (~line 147):

```python
        "mode": "proximal (max_distance_a={})".format(max_distance_a) if max_distance_a else "distal (min_distance_a={})".format(config.ECTODOMAIN_MIN_DISTANCE_A),
```

with:

```python
        "mode": _mode_label(no_distance_filter, max_distance_a),
```

- [ ] **Step 5: Run test to verify it passes**

Run: `pytest tests/test_run_mode_labels.py -v`
Expected: PASS (3 tests)

- [ ] **Step 6: Commit**

```bash
git add epitope_pipeline/run.py tests/test_run_mode_labels.py
git commit -m "Extract _mode_label helper; label whole-ECD runs honestly

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Scoring-summary distance-args helper in run.py

The multi-target `plot_scoring_summary` call hardcodes `"≥80Å"`/distal, which is wrong when the distance filter is off. Extract a pure helper that returns the right `distance_label`/`distance_value`/`distance_mode`.

**Files:**
- Modify: `epitope_pipeline/run.py` (add helper near `_mode_label`; replace the `plot_scoring_summary(...)` kwargs ~lines 477-486)
- Test: `tests/test_run_mode_labels.py` (extend)

- [ ] **Step 1: Write the failing test**

Append to `tests/test_run_mode_labels.py`:

```python
from epitope_pipeline.run import _summary_distance_args


def test_summary_distance_args_whole_ecd():
    args = _summary_distance_args(no_distance_filter=True, max_distance_a=None)
    assert args == {
        "distance_label": "whole ECD",
        "distance_value": None,
        "distance_mode": "whole_ecd",
    }


def test_summary_distance_args_proximal():
    args = _summary_distance_args(no_distance_filter=False, max_distance_a=50.0)
    assert args["distance_mode"] == "proximal"
    assert args["distance_value"] == 50.0
    assert args["distance_label"] == "≤50Å"


def test_summary_distance_args_distal():
    args = _summary_distance_args(no_distance_filter=False, max_distance_a=None)
    assert args["distance_mode"] == "distal"
    assert args["distance_value"] == config.ECTODOMAIN_MIN_DISTANCE_A
    assert args["distance_label"] == "≥80Å"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_run_mode_labels.py -v`
Expected: FAIL with `ImportError: cannot import name '_summary_distance_args'`

- [ ] **Step 3: Add the helper**

Add in `epitope_pipeline/run.py` directly below `_mode_label`:

```python
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
            "distance_label": "≤{:.0f}Å".format(max_distance_a),
            "distance_value": max_distance_a,
            "distance_mode": "proximal",
        }
    return {
        "distance_label": "≥{:.0f}Å".format(config.ECTODOMAIN_MIN_DISTANCE_A),
        "distance_value": config.ECTODOMAIN_MIN_DISTANCE_A,
        "distance_mode": "distal",
    }
```

- [ ] **Step 4: Use the helper in run_pipeline**

In `epitope_pipeline/run.py`, replace the multi-target summary call (~lines 477-486). The current source uses `\u` escapes for the ≤/≥/Å glyphs — match it exactly:

```python
    if len(targets) > 1:
        plot_scoring_summary(
            all_scores=all_scores,
            target_metrics=all_metrics,
            targets=targets,
            output_path=str(figures_dir / "scoring_summary.png"),
            distance_label="\u226440\u00c5" if max_distance_a else "\u226580\u00c5",
            distance_value=max_distance_a or config.ECTODOMAIN_MIN_DISTANCE_A,
            distance_mode="proximal" if max_distance_a else "distal",
        )
```

with:

```python
    if len(targets) > 1:
        plot_scoring_summary(
            all_scores=all_scores,
            target_metrics=all_metrics,
            targets=targets,
            output_path=str(figures_dir / "scoring_summary.png"),
            **_summary_distance_args(no_distance_filter, max_distance_a),
        )
```

- [ ] **Step 5: Run test to verify it passes**

Run: `pytest tests/test_run_mode_labels.py -v`
Expected: PASS (6 tests total)

- [ ] **Step 6: Commit**

```bash
git add epitope_pipeline/run.py tests/test_run_mode_labels.py
git commit -m "Extract _summary_distance_args; fix whole-ECD summary label

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: whole_ecd footer branch in visualize.py

`plot_scoring_summary` builds a footer `dist_str` from `distance_mode`, with only `proximal` vs else (distal). Add a `whole_ecd` branch via a pure helper so the new `distance_mode="whole_ecd"` renders correctly.

**Files:**
- Modify: `epitope_pipeline/viz/visualize.py` (add helper near top of `plot_scoring_summary`; replace the `dist_str` if/else ~lines 650-653)
- Test: `tests/test_viz_footer.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_viz_footer.py
"""Unit test for the distance-footer string helper used by plot_scoring_summary."""

from epitope_pipeline.viz.visualize import _distance_footer_str


def test_footer_whole_ecd():
    assert _distance_footer_str("whole_ecd", None) == (
        "whole ectodomain (no distance filter)"
    )


def test_footer_proximal():
    assert _distance_footer_str("proximal", 50.0) == "≤50Å from membrane"


def test_footer_distal():
    assert _distance_footer_str("distal", 80.0) == "≥80Å from membrane"


def test_footer_proximal_default_value():
    assert _distance_footer_str("proximal", None) == "≤40Å from membrane"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_viz_footer.py -v`
Expected: FAIL with `ImportError: cannot import name '_distance_footer_str'`

- [ ] **Step 3: Add the helper**

Add this module-level function in `epitope_pipeline/viz/visualize.py` immediately above `def plot_scoring_summary(`:

```python
def _distance_footer_str(distance_mode, distance_value):
    """Footer fragment describing the spatial filter for the summary chart."""
    if distance_mode == "whole_ecd":
        return "whole ectodomain (no distance filter)"
    if distance_mode == "proximal":
        return "≤{:.0f}Å from membrane".format(distance_value or 40)
    return "≥{:.0f}Å from membrane".format(distance_value or 80)
```

- [ ] **Step 4: Use the helper in plot_scoring_summary**

In `epitope_pipeline/viz/visualize.py`, replace the footer if/else (~lines 650-653):

```python
    if distance_mode == "proximal":
        dist_str = "\u2264{:.0f}\u00c5 from membrane".format(distance_value or 40)
    else:
        dist_str = "\u2265{:.0f}\u00c5 from membrane".format(distance_value or 80)
```

with:

```python
    dist_str = _distance_footer_str(distance_mode, distance_value)
```

- [ ] **Step 5: Run test to verify it passes**

Run: `pytest tests/test_viz_footer.py -v`
Expected: PASS (4 tests)

- [ ] **Step 6: Commit**

```bash
git add epitope_pipeline/viz/visualize.py tests/test_viz_footer.py
git commit -m "Add whole_ecd footer branch to scoring-summary chart

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: single.py module — run_single wrapper

The core new feature: a thin wrapper that runs the whole-ECD assessment for single proteins.

**Files:**
- Create: `epitope_pipeline/single.py`
- Test: `tests/test_single.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_single.py
"""Tests for the single-protein whole-ECD wrapper run_single."""

from epitope_pipeline import single as single_module


def test_run_single_forces_no_distance_filter(monkeypatch):
    """run_single must delegate to run_pipeline with no_distance_filter=True
    and must NOT pass any min/max distance."""
    captured = {}

    def stub(identifiers, **kwargs):
        captured["identifiers"] = list(identifiers)
        captured.update(kwargs)
        return {"run_dir": "/tmp/fake", "targets": [], "metrics": {}}

    monkeypatch.setattr(single_module, "run_pipeline", stub)

    single_module.run_single(["ERBB2"], cyno_mismatch_percent=20.0,
                             skip_cyno_gate=True)

    assert captured["identifiers"] == ["ERBB2"]
    assert captured["no_distance_filter"] is True
    assert captured["cyno_mismatch_percent"] == 20.0
    assert captured["skip_cyno_gate"] is True
    # The distal/proximal knobs must not be forwarded by the single-protein API
    assert "min_distance_a" not in captured
    assert "max_distance_a" not in captured


def test_run_single_default_run_name_marks_single(monkeypatch):
    """When no run_name is given, the auto slug identifies this as a single run."""
    captured = {}

    def stub(identifiers, **kwargs):
        captured.update(kwargs)
        return {"run_dir": "/tmp/fake", "targets": [], "metrics": {}}

    monkeypatch.setattr(single_module, "run_pipeline", stub)
    single_module.run_single(["ERBB2"])
    assert "_single_erbb2" in captured["run_name"]


def test_run_single_explicit_run_name_passes_through(monkeypatch):
    captured = {}

    def stub(identifiers, **kwargs):
        captured.update(kwargs)
        return {"run_dir": "/tmp/fake", "targets": [], "metrics": {}}

    monkeypatch.setattr(single_module, "run_pipeline", stub)
    single_module.run_single(["ERBB2"], run_name="my_run")
    assert captured["run_name"] == "my_run"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_single.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'epitope_pipeline.single'`

- [ ] **Step 3: Create the module**

Create `epitope_pipeline/single.py`:

```python
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

import sys

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

    Args mirror run_pipeline's non-distance arguments. See run_pipeline for the
    structure of the returned dict.
    """
    if run_name is None:
        from datetime import datetime
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
    from datetime import datetime  # noqa: F401  (kept parallel to run.py)

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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_single.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add epitope_pipeline/single.py tests/test_single.py
git commit -m "Add run_single: whole-ECD suitability wrapper over run_pipeline

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: single.py CLI argparse tests

Verify the CLI parses targets/aliases/flags and rejects the distance flags (the whole point of the separate entry point).

**Files:**
- Test: `tests/test_single_cli.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_single_cli.py
"""Argparse / CLI tests for epitope_pipeline.single.main.

run_single is monkeypatched to a recording stub so flags are exercised
without running the actual pipeline.
"""

import sys

import pytest

from epitope_pipeline import single as single_module


@pytest.fixture
def captured(monkeypatch):
    calls = []

    def stub(identifiers, **kwargs):
        calls.append({"identifiers": list(identifiers), **kwargs})
        return {"run_dir": "/tmp/fake", "targets": [], "metrics": {}}

    monkeypatch.setattr(single_module, "run_single", stub)
    return calls


def _invoke(monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", ["epitope-single", *argv])
    single_module.main()


def test_positional_targets(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "EGFR"])
    assert captured[0]["identifiers"] == ["ERBB2", "EGFR"]
    assert captured[0]["skip_cyno_gate"] is False
    assert captured[0]["cyno_mismatch_percent"] is None


def test_aliases_and_flags(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2=HER2", "--no-cyno", "--nonspecific-percent", "70"])
    c = captured[0]
    assert c["identifiers"] == ["ERBB2"]
    assert c["aliases"] == {"ERBB2": "HER2"}
    assert c["skip_cyno_gate"] is True
    assert c["nonspecific_percent"] == 70.0


def test_targets_file_merges_and_dedupes(captured, monkeypatch, tmp_path):
    f = tmp_path / "list.txt"
    f.write_text("EGFR\nMSLN\n")
    _invoke(monkeypatch, ["ERBB2", "EGFR", "--targets-file", str(f)])
    assert captured[0]["identifiers"] == ["ERBB2", "EGFR", "MSLN"]


def test_no_targets_errors(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["epitope-single"])
    with pytest.raises(SystemExit):
        single_module.main()
    assert "no targets" in capsys.readouterr().err.lower()


@pytest.mark.parametrize("flag", ["--min-distance", "--max-distance", "--no-distance-filter"])
def test_distance_flags_rejected(monkeypatch, capsys, flag):
    """The single-protein CLI must NOT accept distal/proximal distance flags."""
    argv = ["epitope-single", "ERBB2", flag]
    if flag != "--no-distance-filter":
        argv.append("50")
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SystemExit):
        single_module.main()
    assert "unrecognized arguments" in capsys.readouterr().err.lower()
```

- [ ] **Step 2: Run test to verify it fails or passes**

Run: `pytest tests/test_single_cli.py -v`
Expected: PASS (the module from Task 4 already provides `main`). If any test fails, fix `single.py` to match.

- [ ] **Step 3: Commit**

```bash
git add tests/test_single_cli.py
git commit -m "Add CLI tests for epitope-single (distance flags rejected)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Register the console script

**Files:**
- Modify: `pyproject.toml:27-29`

- [ ] **Step 1: Add the entry point**

In `pyproject.toml`, under `[project.scripts]`, replace:

```toml
[project.scripts]
epitope-pipeline = "epitope_pipeline.run:main"
epitope-bispecific = "epitope_pipeline.bispecific:main"
```

with:

```toml
[project.scripts]
epitope-pipeline = "epitope_pipeline.run:main"
epitope-bispecific = "epitope_pipeline.bispecific:main"
epitope-single = "epitope_pipeline.single:main"
```

- [ ] **Step 2: Reinstall so the console script registers**

Run: `.venv/bin/pip install -e . --quiet && .venv/bin/epitope-single --help`
Expected: help text prints; the options list shows NO `--min-distance`, `--max-distance`, or `--no-distance-filter`.

- [ ] **Step 3: Commit**

```bash
git add pyproject.toml
git commit -m "Register epitope-single console script

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: Full suite + README docs

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Run the whole non-slow test suite**

Run: `.venv/bin/pytest -q`
Expected: all tests pass (existing + the new mode-label, viz-footer, single, single-CLI tests).

- [ ] **Step 2: Add README section**

In `README.md`, add a "Single-protein ECD suitability" section near the batch/bispecific usage sections. Use this content (adapt heading level to match surrounding sections):

```markdown
## Single-protein ECD suitability

To assess whether one protein's **whole extracellular domain** is a suitable
VHH-epitope target — with no distal/proximal zone framing — use the
`epitope-single` entry point:

```bash
# One target
epitope-single ERBB2

# Several targets, with display aliases and the cyno gate relaxed
epitope-single ERBB2=HER2 EGFR MSLN --cyno-mismatch-percent 25

# From a file (one IDENT[=ALIAS] per line)
epitope-single --targets-file targets.txt

# Or as a module
python -m epitope_pipeline.single ERBB2
```

This evaluates the entire ectodomain surface for druggable, cyno-conserved,
human-specific epitope patches and writes the same run artifacts (figures,
summary CSV, exports) as the standard pipeline. Suitability is read off the
per-target patch counts and scores in the summary outputs.

Unlike `epitope-bispecific` (which scores each target in a distal ≥60 Å and a
proximal ≤50 Å zone for complementary antibody arms), `epitope-single` applies
no distance filter — the run parameters record the mode as
`whole ectodomain (no distance filter)`.
```

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "Document single-protein ECD suitability mode

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Final verification

- [ ] `.venv/bin/pytest -q` — all green.
- [ ] `.venv/bin/epitope-single --help` — no distance flags listed.
- [ ] (Optional, network) `.venv/bin/epitope-single ERBB2 --run-name single_smoke` then confirm `runs/single_smoke/Supplementary Files/.../` parameters record `"mode": "whole ectodomain (no distance filter)"`.
