# Epitope Pipeline — Package Reorg Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reorganize the flat top-level pipeline into a pip-installable package (`epitope_pipeline`) with `compute/`, `io/`, and `viz/` subpackages, drop the Tamarind structure-prediction fallback, and add a tests/ harness — without changing any pipeline logic.

**Architecture:** Repo root becomes the project root (holds `pyproject.toml`, `README.md`, `tests/`, `scripts/`, `diagnostics/`). All Python source moves into a nested `epitope_pipeline/` directory and is grouped into three subpackages by responsibility: `compute/` (numerical pipeline stages), `io/` (PDB parsing + network + file exports), `viz/` (figures). `run.py`, `bispecific.py`, `app.py`, `config.py`, `utils.py`, and `__init__.py` stay at package root. Tamarind-related code is deleted (structure prediction fallback will be reintroduced later via Boltz-2).

**Tech Stack:** Python ≥3.10, setuptools/pyproject.toml, pytest. Existing deps unchanged: numpy, scipy, biopython, freesasa, requests, matplotlib, openpyxl. Drop: `python-dotenv` (was Tamarind-only).

---

## Final layout

```
epitope_pipeline/                      (repo root, becomes project root)
├── pyproject.toml                     NEW
├── README.md                          (small import-path updates)
├── .gitignore                         (unchanged)
├── epitope_pipeline/                  NEW package directory (everything moves into here)
│   ├── __init__.py
│   ├── config.py
│   ├── utils.py                       (logging + empty_metric only)
│   ├── run.py                         (run_pipeline + CLI)
│   ├── bispecific.py
│   ├── app.py                         (streamlit)
│   ├── compute/
│   │   ├── __init__.py
│   │   ├── spatial.py
│   │   ├── surface.py
│   │   ├── conservation.py
│   │   ├── specificity.py
│   │   └── scoring.py
│   ├── io/
│   │   ├── __init__.py
│   │   ├── pdb.py                     NEW (extract_ca_coords + get_chain, lifted from utils.py)
│   │   ├── targets.py                 (was target_input.py)
│   │   ├── structure.py               (Tamarind blocks removed)
│   │   ├── membrane.py
│   │   ├── export.py
│   │   └── export_bispecific.py
│   └── viz/
│       ├── __init__.py
│       ├── visualize.py
│       └── visualize_bispecific.py
├── tests/                             NEW
│   ├── __init__.py
│   ├── conftest.py
│   └── test_imports.py
├── scripts/                           (unchanged location, internal imports updated)
├── diagnostics/                       (unchanged location, internal imports updated)
└── archive/                           (unchanged, no import updates — pre-archived)
```

**External-callsite-impact summary** (clean break, no shim):
- `app.py`: imports `from epitope_pipeline.run import run_pipeline` → unchanged (still at package root after move).
- `diagnostics/*.py`: must update any deep imports like `from epitope_pipeline.specificity import …` → `from epitope_pipeline.compute.specificity import …`.
- `scripts/*.py`: same.
- `target_input` is renamed to `targets`; any imports of it must be updated.
- `extract_ca_coords` and `get_chain` move from `utils` to `io.pdb`; ~10 callsites to update.

---

## Task 1: Create branch

**Files:** none

- [ ] **Step 1: Create the feature branch**

```bash
git -C /Users/jverboon/Analyses/epitope_pipeline checkout -b refactor/package-layout
```

- [ ] **Step 2: Verify clean state**

```bash
git -C /Users/jverboon/Analyses/epitope_pipeline status
```

Expected: "On branch refactor/package-layout" and "nothing to commit, working tree clean".

---

## Task 2: Add pyproject.toml at repo root

**Files:**
- Create: `pyproject.toml`
- Create: `MANIFEST.in` (only if needed for non-Python data; skip otherwise)

This task pins the package shape *before* the file moves. After this commit `pip install -e .` won't work yet (because the package is still at the repo root, not nested). That's expected — Task 3 fixes it.

- [ ] **Step 1: Create pyproject.toml**

```toml
# pyproject.toml
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "epitope-pipeline"
version = "0.2.0"
description = "Structural bioinformatics pipeline for VHH epitope discovery on human membrane proteins."
readme = "README.md"
requires-python = ">=3.10"
authors = [{name = "Jeff Verboon"}]
dependencies = [
    "numpy>=1.24,<2.0",
    "scipy>=1.10",
    "biopython>=1.81",
    "freesasa>=2.2.0",
    "requests>=2.28",
    "matplotlib>=3.7.1,<3.8",
    "openpyxl>=3.0",
]

[project.optional-dependencies]
app = ["streamlit>=1.28"]
dev = ["pytest>=8.0"]

[project.scripts]
epitope-pipeline = "epitope_pipeline.run:main"
epitope-bispecific = "epitope_pipeline.bispecific:main"

[tool.setuptools.packages.find]
include = ["epitope_pipeline*"]
exclude = ["tests*", "diagnostics*", "scripts*", "archive*"]
```

- [ ] **Step 2: Stage and commit**

```bash
git -C /Users/jverboon/Analyses/epitope_pipeline add pyproject.toml
git -C /Users/jverboon/Analyses/epitope_pipeline commit -m "Add pyproject.toml for installable package"
```

Note: `[project.scripts]` references `epitope_pipeline.run:main` and `epitope_pipeline.bispecific:main`. Those `main()` functions don't exist yet — Tasks 3 and 9 wrap the existing `__main__` blocks in `def main()` so the entry points work. Until then, `epitope-pipeline` CLI won't run; `python -m epitope_pipeline.run` still does.

---

## Task 3: Move Python sources into nested `epitope_pipeline/` directory

**Files:** every `.py` file currently at repo root, plus `__init__.py`.

Goal: `git mv` each file from `<repo>/` into `<repo>/epitope_pipeline/`. After this commit imports still resolve to `epitope_pipeline.config`, `epitope_pipeline.run`, etc. — but only if the package is installed in editable mode (Task 4).

- [ ] **Step 1: Create the nested package directory**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
mkdir epitope_pipeline
```

- [ ] **Step 2: Move every `.py` file at repo root into the nested dir**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
git mv __init__.py epitope_pipeline/__init__.py
git mv app.py epitope_pipeline/app.py
git mv bispecific.py epitope_pipeline/bispecific.py
git mv config.py epitope_pipeline/config.py
git mv conservation.py epitope_pipeline/conservation.py
git mv export.py epitope_pipeline/export.py
git mv export_bispecific.py epitope_pipeline/export_bispecific.py
git mv membrane.py epitope_pipeline/membrane.py
git mv run.py epitope_pipeline/run.py
git mv scoring.py epitope_pipeline/scoring.py
git mv spatial.py epitope_pipeline/spatial.py
git mv specificity.py epitope_pipeline/specificity.py
git mv structure.py epitope_pipeline/structure.py
git mv surface.py epitope_pipeline/surface.py
git mv target_input.py epitope_pipeline/target_input.py
git mv utils.py epitope_pipeline/utils.py
git mv visualize.py epitope_pipeline/visualize.py
git mv visualize_bispecific.py epitope_pipeline/visualize_bispecific.py
git mv requirements.txt epitope_pipeline/requirements.txt  # will be deleted in Task 5
```

- [ ] **Step 3: Wrap CLI entry points in `main()` functions for `[project.scripts]` to work**

Edit `epitope_pipeline/run.py:453-474` (the `if __name__ == "__main__":` block). Change:

```python
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python -m epitope_pipeline.run TARGET1 [TARGET2 ...]")
        ...
        sys.exit(1)

    identifiers = sys.argv[1:]
    results = run_pipeline(identifiers, max_distance_a=config.PROXIMAL_MAX_DISTANCE_A)
    ...
```

to:

```python
def main():
    if len(sys.argv) < 2:
        print("Usage: python -m epitope_pipeline.run TARGET1 [TARGET2 ...]")
        print("")
        print("Targets can be UniProt IDs (e.g. P04626) or gene names (e.g. ERBB2)")
        print("")
        print("Example:")
        print("  python -m epitope_pipeline.run ERBB2 EGFR")
        sys.exit(1)

    identifiers = sys.argv[1:]
    results = run_pipeline(identifiers, max_distance_a=config.PROXIMAL_MAX_DISTANCE_A)

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

- [ ] **Step 4: Do the same for `bispecific.py`**

Locate the `if __name__ == "__main__":` block at the bottom of `epitope_pipeline/bispecific.py`. Lift its body into a `def main():` and call `main()` from the guard. Use:

```bash
grep -n "__main__" /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/bispecific.py
```

to find the line, then apply the same transform as Step 3.

- [ ] **Step 5: Verify the file moves**

```bash
ls /Users/jverboon/Analyses/epitope_pipeline/*.py 2>/dev/null
```

Expected: no output (all `.py` files moved out of repo root).

```bash
ls /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/*.py | wc -l
```

Expected: 18 (16 original modules + `__init__.py` + nothing else).

- [ ] **Step 6: Commit**

```bash
git -C /Users/jverboon/Analyses/epitope_pipeline add -A
git -C /Users/jverboon/Analyses/epitope_pipeline commit -m "Move package sources into nested epitope_pipeline/ directory"
```

---

## Task 4: Editable install + import smoke test

**Files:**
- Create: `tests/__init__.py` (empty)
- Create: `tests/conftest.py`
- Create: `tests/test_imports.py`

- [ ] **Step 1: Install in editable mode**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
pip install -e ".[dev]"
```

Expected: "Successfully installed epitope-pipeline-0.2.0" plus pytest deps.

- [ ] **Step 2: Create tests directory**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
mkdir -p tests
touch tests/__init__.py
```

- [ ] **Step 3: Write conftest.py**

```python
# tests/conftest.py
"""Pytest configuration for the epitope_pipeline test suite."""
```

- [ ] **Step 4: Write test_imports.py**

This is the cheap "did the moves break anything" guard. It will be expanded after each subsequent task to cover new module paths.

```python
# tests/test_imports.py
"""Smoke test: every public module imports cleanly.

Catches broken imports introduced by file moves and refactors. Cheap to run
(~0.5 s) — run after every reorg commit.
"""
import importlib

PUBLIC_MODULES = [
    "epitope_pipeline",
    "epitope_pipeline.config",
    "epitope_pipeline.utils",
    "epitope_pipeline.run",
    "epitope_pipeline.bispecific",
    "epitope_pipeline.spatial",
    "epitope_pipeline.surface",
    "epitope_pipeline.conservation",
    "epitope_pipeline.specificity",
    "epitope_pipeline.scoring",
    "epitope_pipeline.target_input",
    "epitope_pipeline.structure",
    "epitope_pipeline.membrane",
    "epitope_pipeline.export",
    "epitope_pipeline.export_bispecific",
    "epitope_pipeline.visualize",
    "epitope_pipeline.visualize_bispecific",
]


def test_all_modules_import():
    for name in PUBLIC_MODULES:
        importlib.import_module(name)
```

- [ ] **Step 5: Run the test**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
pytest tests/test_imports.py -v
```

Expected: PASS for `test_all_modules_import` (provided `tamarind/` is reachable; if it isn't, this test will fail on `epitope_pipeline.structure` — proceed to Task 5 to fix permanently). If the test fails for any other reason, stop and investigate before continuing.

- [ ] **Step 6: Commit**

```bash
git -C /Users/jverboon/Analyses/epitope_pipeline add tests/
git -C /Users/jverboon/Analyses/epitope_pipeline commit -m "Add pytest smoke test for module imports"
```

---

## Task 5: Remove Tamarind structure-prediction fallback

**Files:**
- Modify: `epitope_pipeline/structure.py`
- Modify: `epitope_pipeline/config.py`
- Modify: `epitope_pipeline/run.py`
- Delete: `epitope_pipeline/requirements.txt` (replaced by `pyproject.toml`)

Drops a hard sibling-directory dependency. After this commit AlphaFold DB is the only structure source for proteins without good experimental PDB; targets that miss both raise `StructureAcquisitionError`. Boltz-2 sequence-based folding can be reintroduced later behind a real package boundary.

- [ ] **Step 1: Strip Tamarind from `structure.py`**

In `epitope_pipeline/structure.py`:

a) Delete lines 44-46 (the sys.path hack and Tamarind import):

```python
# Import TamarindClient from sibling project
sys.path.insert(0, str(Path(__file__).parent.parent))
from tamarind.client import TamarindClient, TamarindError
```

b) Delete the `sys` and unused `Path` parent reference if no other use remains. Run:

```bash
grep -n "^import sys\|sys\." /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/structure.py
```

Keep `import sys` only if other code in the file still uses it. Otherwise remove.

c) Replace the call site at lines 144-146:

```python
    # Fallback: Tamarind AlphaFold prediction (de novo)
    logger.info("  Folding %s via Tamarind %s...", target.gene_name, FOLDING_TOOL)
    return _fold_with_tamarind(target, structures_dir)
```

with:

```python
    raise StructureAcquisitionError(
        "No structure available for {} ({}): not in AlphaFold DB and no "
        "experimental PDB. Sequence-based folding has been removed; see "
        "the planned Boltz-2 integration.".format(target.gene_name, target.uniprot_id)
    )
```

d) Delete the entire "Tamarind AlphaFold fallback" section starting at the comment block around line 626 through the end of `_download_fold_result` (around line 771). This removes:
- `_fold_with_tamarind`
- `_download_fold_result`
- Their imports of `zipfile`, `io`, `shutil` if not used elsewhere — verify with `grep -n "zipfile\|^import io\|^import shutil" structure.py` and remove unused imports only.

e) Remove the `FOLDING_NUM_MODELS, FOLDING_NUM_RECYCLES, FOLDING_TOOL, FOLDING_TIMEOUT, RAW_JOBS_DIR` references from the imports at the top of structure.py (lines 26-42). Verify with grep; only delete those that are no longer referenced after the function deletions.

- [ ] **Step 2: Strip Tamarind from `config.py`**

In `epitope_pipeline/config.py`:

a) Delete lines 5 and 19-21 (the docstring mention and the dotenv loading):

```python
Tamarind API key is loaded from the existing tamarind/.env file.
```
```python
# Load Tamarind API key from existing .env
load_dotenv(TAMARIND_ROOT / ".env")
TAMARIND_API_KEY = os.environ.get("TAMARIND_API_KEY", "")
```

b) Delete lines 17 (`TAMARIND_ROOT = ...`).

c) Delete lines 56-59:

```python
FOLDING_TOOL = "alphafold"          # Tamarind folding service name
FOLDING_NUM_MODELS = 1
FOLDING_NUM_RECYCLES = 3
FOLDING_TIMEOUT = 7200              # 2 hours max wait
```

d) Delete the `from dotenv import load_dotenv` import on line 10.

e) Delete the `RAW_JOBS_DIR = PIPELINE_ROOT / "raw_jobs"` line (line 25 — was Tamarind-specific scratch space). Confirm it's unused: `grep -rn "RAW_JOBS_DIR" /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/`.

f) The `PROJECT_ROOT = PIPELINE_ROOT.parent` line (line 16) was only used by `TAMARIND_ROOT`; delete it.

g) Verify config.py still imports cleanly:

```bash
python -c "from epitope_pipeline import config; print('ok')"
```

Expected: `ok`.

- [ ] **Step 3: Strip Tamarind from `run.py`**

In `epitope_pipeline/run.py`:

a) Remove `folding_tool=None,` from the `run_pipeline` signature (around line 60).

b) Remove the corresponding docstring line (around line 75): `folding_tool: Override Tamarind folding tool (default "alphafold").`

c) Remove the override block (around lines 95-96):

```python
if folding_tool is not None:
    config.FOLDING_TOOL = folding_tool
```

d) Remove the `parameters["folding_tool"]` entry (around line 143).

- [ ] **Step 4: Delete the now-redundant `requirements.txt`**

```bash
git -C /Users/jverboon/Analyses/epitope_pipeline rm epitope_pipeline/requirements.txt
```

- [ ] **Step 5: Run the import smoke test**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
pytest tests/test_imports.py -v
```

Expected: PASS. If it fails because of a leftover Tamarind reference, grep for it:

```bash
grep -rn -i "tamarind\|TAMARIND\|FOLDING_" /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/
```

Fix every hit.

- [ ] **Step 6: Update README**

Search the README for Tamarind/folding mentions:

```bash
grep -n -i "tamarind\|folding" /Users/jverboon/Analyses/epitope_pipeline/README.md
```

Update or delete those lines (the structure-acquisition section near line 64 mentions "Tamarind AlphaFold (fallback)" — replace with "AlphaFold DB" or "experimental PDB" only).

- [ ] **Step 7: Commit**

```bash
git -C /Users/jverboon/Analyses/epitope_pipeline add -A
git -C /Users/jverboon/Analyses/epitope_pipeline commit -m "Remove Tamarind structure-prediction fallback"
```

---

## Task 6: Carve out `io/` subpackage

**Files:**
- Create: `epitope_pipeline/io/__init__.py`
- Create: `epitope_pipeline/io/pdb.py`
- Move: `epitope_pipeline/target_input.py` → `epitope_pipeline/io/targets.py`
- Move: `epitope_pipeline/structure.py` → `epitope_pipeline/io/structure.py`
- Move: `epitope_pipeline/membrane.py` → `epitope_pipeline/io/membrane.py`
- Move: `epitope_pipeline/export.py` → `epitope_pipeline/io/export.py`
- Move: `epitope_pipeline/export_bispecific.py` → `epitope_pipeline/io/export_bispecific.py`
- Modify: `epitope_pipeline/utils.py` (carve out `extract_ca_coords` + `get_chain`)
- Modify: every internal callsite of moved modules and of `extract_ca_coords` / `get_chain`

The `io/` subpackage owns everything that touches the filesystem or the network: PDB parsing, UniProt/RCSB/AlphaFold-DB/OPM HTTP, file exports.

- [ ] **Step 1: Create the io subpackage**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
mkdir -p epitope_pipeline/io
```

```python
# epitope_pipeline/io/__init__.py
"""I/O subpackage: PDB parsing, network clients, file exports."""
```

- [ ] **Step 2: Carve `extract_ca_coords` + `get_chain` out of `utils.py` into `io/pdb.py`**

Create `epitope_pipeline/io/pdb.py`:

```python
"""
PDB-file helpers: chain access and Cα-coordinate extraction.

These were previously in utils.py — moved here because they're I/O concerns
(BioPython parsing) rather than general-purpose utilities.
"""

import numpy as np
from Bio.PDB import PDBParser


def get_chain(model, chain_id):
    """Get chain by ID with fallback to first chain."""
    if chain_id in model:
        return model[chain_id]
    chains = list(model.get_chains())
    if chains:
        return chains[0]
    raise ValueError("No chains found in structure")


def extract_ca_coords(pdb_path, chain_id):
    """
    Extract Cα coordinates from a PDB file.

    Args:
        pdb_path: Path to PDB file.
        chain_id: Target chain ID.

    Returns:
        Dict {residue_number: np.array([x, y, z])}.
    """
    parser = PDBParser(QUIET=True)
    bio_struct = parser.get_structure("target", pdb_path)
    model = bio_struct[0]
    chain = get_chain(model, chain_id)

    ca_coords = {}
    for res in chain:
        res_id = res.get_id()
        if res_id[0] != " ":
            continue
        if "CA" in res:
            ca_coords[res_id[1]] = res["CA"].get_vector().get_array()

    return ca_coords
```

Then strip those two functions out of `epitope_pipeline/utils.py`. After the strip, `utils.py` should contain only `setup_logging` and `empty_metric` (plus the module docstring and logger setup). Remove `numpy as np` and `from Bio.PDB import PDBParser` from utils.py if no longer used.

- [ ] **Step 3: Move modules into `io/`**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
git mv epitope_pipeline/target_input.py epitope_pipeline/io/targets.py
git mv epitope_pipeline/structure.py epitope_pipeline/io/structure.py
git mv epitope_pipeline/membrane.py epitope_pipeline/io/membrane.py
git mv epitope_pipeline/export.py epitope_pipeline/io/export.py
git mv epitope_pipeline/export_bispecific.py epitope_pipeline/io/export_bispecific.py
```

- [ ] **Step 4: Update intra-`io/` imports**

Inside the moved files, any `from epitope_pipeline.target_input import …` becomes `from epitope_pipeline.io.targets import …`, and similar. Run:

```bash
grep -rn "from epitope_pipeline\.\(target_input\|structure\|membrane\|export\|export_bispecific\) " /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/io/
```

For each hit, rewrite the import to use the `epitope_pipeline.io.<module>` path.

Same for `from epitope_pipeline.utils import extract_ca_coords` and `from epitope_pipeline.utils import get_chain` — both move to `epitope_pipeline.io.pdb`. Find them:

```bash
grep -rn "extract_ca_coords\|get_chain" /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/
```

Update each callsite. Pay special attention to the deferred imports (`from epitope_pipeline.utils import extract_ca_coords` inside function bodies) — there are several in `export.py` and `spatial.py`.

- [ ] **Step 5: Update callsites in `run.py` and `bispecific.py` (still at package root)**

```bash
grep -n "from epitope_pipeline\." /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/run.py /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/bispecific.py
```

Rewrite each import to point to the new locations:
- `epitope_pipeline.target_input` → `epitope_pipeline.io.targets`
- `epitope_pipeline.structure` → `epitope_pipeline.io.structure`
- `epitope_pipeline.membrane` → `epitope_pipeline.io.membrane`
- `epitope_pipeline.export` → `epitope_pipeline.io.export`
- `epitope_pipeline.utils.extract_ca_coords` → `epitope_pipeline.io.pdb.extract_ca_coords`

- [ ] **Step 6: Update test_imports.py to cover new paths**

In `tests/test_imports.py`, replace:
```python
    "epitope_pipeline.target_input",
    "epitope_pipeline.structure",
    "epitope_pipeline.membrane",
    "epitope_pipeline.export",
    "epitope_pipeline.export_bispecific",
```
with:
```python
    "epitope_pipeline.io",
    "epitope_pipeline.io.pdb",
    "epitope_pipeline.io.targets",
    "epitope_pipeline.io.structure",
    "epitope_pipeline.io.membrane",
    "epitope_pipeline.io.export",
    "epitope_pipeline.io.export_bispecific",
```

- [ ] **Step 7: Run tests**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
pytest tests/test_imports.py -v
```

Expected: PASS. If any module fails to import with a `ModuleNotFoundError`, the failing line tells you which import was missed; fix and re-run.

- [ ] **Step 8: Commit**

```bash
git -C /Users/jverboon/Analyses/epitope_pipeline add -A
git -C /Users/jverboon/Analyses/epitope_pipeline commit -m "Carve I/O modules into epitope_pipeline.io subpackage"
```

---

## Task 7: Carve out `compute/` subpackage

**Files:**
- Create: `epitope_pipeline/compute/__init__.py`
- Move: `epitope_pipeline/spatial.py` → `epitope_pipeline/compute/spatial.py`
- Move: `epitope_pipeline/surface.py` → `epitope_pipeline/compute/surface.py`
- Move: `epitope_pipeline/conservation.py` → `epitope_pipeline/compute/conservation.py`
- Move: `epitope_pipeline/specificity.py` → `epitope_pipeline/compute/specificity.py`
- Move: `epitope_pipeline/scoring.py` → `epitope_pipeline/compute/scoring.py`
- Modify: every callsite of those modules

- [ ] **Step 1: Create the compute subpackage**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
mkdir -p epitope_pipeline/compute
```

```python
# epitope_pipeline/compute/__init__.py
"""Compute subpackage: numerical pipeline stages (spatial, surface, conservation, specificity, scoring)."""
```

- [ ] **Step 2: Move modules**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
git mv epitope_pipeline/spatial.py epitope_pipeline/compute/spatial.py
git mv epitope_pipeline/surface.py epitope_pipeline/compute/surface.py
git mv epitope_pipeline/conservation.py epitope_pipeline/compute/conservation.py
git mv epitope_pipeline/specificity.py epitope_pipeline/compute/specificity.py
git mv epitope_pipeline/scoring.py epitope_pipeline/compute/scoring.py
```

- [ ] **Step 3: Update intra-`compute/` and cross-package imports**

```bash
grep -rn "from epitope_pipeline\.\(spatial\|surface\|conservation\|specificity\|scoring\) " /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/
```

Each hit becomes `from epitope_pipeline.compute.<module>`. Note that some modules import from siblings within compute (e.g., `surface.py` may reference `spatial.py`); update those too.

Also check for relative imports inside the moved files: `from .surface import …` style — those still work after the move (siblings in the same package). Make sure none of them reach back up — search:

```bash
grep -rn "^from \." /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/compute/
```

If any of those references something now in `io/` or at the package root, change to an absolute import.

- [ ] **Step 4: Update test_imports.py**

Replace:
```python
    "epitope_pipeline.spatial",
    "epitope_pipeline.surface",
    "epitope_pipeline.conservation",
    "epitope_pipeline.specificity",
    "epitope_pipeline.scoring",
```
with:
```python
    "epitope_pipeline.compute",
    "epitope_pipeline.compute.spatial",
    "epitope_pipeline.compute.surface",
    "epitope_pipeline.compute.conservation",
    "epitope_pipeline.compute.specificity",
    "epitope_pipeline.compute.scoring",
```

- [ ] **Step 5: Run tests**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
pytest tests/test_imports.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git -C /Users/jverboon/Analyses/epitope_pipeline add -A
git -C /Users/jverboon/Analyses/epitope_pipeline commit -m "Carve compute modules into epitope_pipeline.compute subpackage"
```

---

## Task 8: Carve out `viz/` subpackage

**Files:**
- Create: `epitope_pipeline/viz/__init__.py`
- Move: `epitope_pipeline/visualize.py` → `epitope_pipeline/viz/visualize.py`
- Move: `epitope_pipeline/visualize_bispecific.py` → `epitope_pipeline/viz/visualize_bispecific.py`

- [ ] **Step 1: Create the viz subpackage**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
mkdir -p epitope_pipeline/viz
```

```python
# epitope_pipeline/viz/__init__.py
"""Viz subpackage: matplotlib figure generation."""
```

- [ ] **Step 2: Move modules**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
git mv epitope_pipeline/visualize.py epitope_pipeline/viz/visualize.py
git mv epitope_pipeline/visualize_bispecific.py epitope_pipeline/viz/visualize_bispecific.py
```

- [ ] **Step 3: Update callsites**

```bash
grep -rn "from epitope_pipeline\.\(visualize\|visualize_bispecific\) \|import epitope_pipeline\.\(visualize\|visualize_bispecific\)" /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/
```

Each becomes `from epitope_pipeline.viz.<module>`.

- [ ] **Step 4: Update test_imports.py**

Replace:
```python
    "epitope_pipeline.visualize",
    "epitope_pipeline.visualize_bispecific",
```
with:
```python
    "epitope_pipeline.viz",
    "epitope_pipeline.viz.visualize",
    "epitope_pipeline.viz.visualize_bispecific",
```

- [ ] **Step 5: Run tests**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
pytest tests/test_imports.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git -C /Users/jverboon/Analyses/epitope_pipeline add -A
git -C /Users/jverboon/Analyses/epitope_pipeline commit -m "Carve viz modules into epitope_pipeline.viz subpackage"
```

---

## Task 9: Update external callsites — diagnostics, scripts, app.py, README

**Files:**
- Modify: every `.py` under `diagnostics/`
- Modify: every `.py` under `scripts/`
- Modify: `epitope_pipeline/app.py`
- Modify: `README.md`

- [ ] **Step 1: Inventory external import usages**

```bash
grep -rn "from epitope_pipeline\.\|import epitope_pipeline\." \
  /Users/jverboon/Analyses/epitope_pipeline/diagnostics/ \
  /Users/jverboon/Analyses/epitope_pipeline/scripts/ \
  /Users/jverboon/Analyses/epitope_pipeline/epitope_pipeline/app.py
```

Save the output — every line is a callsite to update.

- [ ] **Step 2: Apply the rewrites**

For each line in the inventory, apply the mapping:

| Old import path                              | New import path                              |
|----------------------------------------------|----------------------------------------------|
| `epitope_pipeline.target_input`              | `epitope_pipeline.io.targets`                |
| `epitope_pipeline.structure`                 | `epitope_pipeline.io.structure`              |
| `epitope_pipeline.membrane`                  | `epitope_pipeline.io.membrane`               |
| `epitope_pipeline.export`                    | `epitope_pipeline.io.export`                 |
| `epitope_pipeline.export_bispecific`         | `epitope_pipeline.io.export_bispecific`      |
| `epitope_pipeline.spatial`                   | `epitope_pipeline.compute.spatial`           |
| `epitope_pipeline.surface`                   | `epitope_pipeline.compute.surface`           |
| `epitope_pipeline.conservation`              | `epitope_pipeline.compute.conservation`      |
| `epitope_pipeline.specificity`               | `epitope_pipeline.compute.specificity`       |
| `epitope_pipeline.scoring`                   | `epitope_pipeline.compute.scoring`           |
| `epitope_pipeline.visualize`                 | `epitope_pipeline.viz.visualize`             |
| `epitope_pipeline.visualize_bispecific`      | `epitope_pipeline.viz.visualize_bispecific`  |
| `epitope_pipeline.utils.extract_ca_coords`   | `epitope_pipeline.io.pdb.extract_ca_coords`  |
| `epitope_pipeline.utils.get_chain`           | `epitope_pipeline.io.pdb.get_chain`          |

`epitope_pipeline.config`, `epitope_pipeline.utils.setup_logging`, `epitope_pipeline.utils.empty_metric`, `epitope_pipeline.run`, `epitope_pipeline.bispecific` are unchanged.

- [ ] **Step 3: Update README import examples**

```bash
grep -n "from epitope_pipeline\.\|import epitope_pipeline\." /Users/jverboon/Analyses/epitope_pipeline/README.md
```

Apply the same mapping to each example. Also update the "Quick Start" section: the install instruction should now include `pip install -e .` instead of `pip install -r requirements.txt`.

- [ ] **Step 4: Verify each external file still imports**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
python -c "import epitope_pipeline.app"
for f in diagnostics/*.py scripts/*.py; do
  python -c "import ast; ast.parse(open('$f').read())" || echo "PARSE FAIL: $f"
done
```

Expected: no `PARSE FAIL` lines and no ImportError. (We don't actually run the scripts — they have side effects — just confirm they parse.)

- [ ] **Step 5: Run smoke tests**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
pytest tests/test_imports.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git -C /Users/jverboon/Analyses/epitope_pipeline add -A
git -C /Users/jverboon/Analyses/epitope_pipeline commit -m "Update external callsites for new package layout"
```

---

## Task 10: End-to-end smoke test on a cached target

**Files:**
- Create: `tests/test_pipeline_smoke.py`

Confirms the full pipeline still runs. Uses a target whose BLAST + structure are likely already cached (`P04626` / ERBB2 was the canonical example in the README). If the caches aren't there this test will hit the network — that's intentional once but you can skip it on subsequent runs with `pytest -k "not pipeline_smoke"`.

- [ ] **Step 1: Write the smoke test**

```python
# tests/test_pipeline_smoke.py
"""End-to-end smoke test on a small fixture target.

Runs the full pipeline (target → structure → membrane → spatial → surface →
conservation → specificity → scoring → export → viz) on one cached UniProt ID
and asserts the run produced *some* output. Network-dependent on first run;
cached afterwards.
"""
import pytest

from epitope_pipeline.run import run_pipeline


@pytest.mark.slow
def test_pipeline_runs_end_to_end(tmp_path, monkeypatch):
    """Smoke: ERBB2 produces a populated metrics dict."""
    # Redirect runs to the temporary directory so we don't pollute repo state
    from epitope_pipeline import config as cfg
    monkeypatch.setattr(cfg, "RUNS_DIR", tmp_path / "runs")

    result = run_pipeline(["P04626"], verbose=False)
    assert "targets" in result
    assert len(result["targets"]) == 1
    metrics = result["metrics"][result["targets"][0].uniprot_id]
    # Run completed; we don't care what the score is, just that we got one
    assert "n_patches" in metrics
    assert "best_score" in metrics
```

- [ ] **Step 2: Add the slow marker to pyproject.toml**

In `pyproject.toml`, append:

```toml
[tool.pytest.ini_options]
markers = [
    "slow: end-to-end tests that may hit network on first run",
]
```

- [ ] **Step 3: Run the fast tests (default)**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
pytest tests/ -v -m "not slow"
```

Expected: only `test_imports.py` runs; PASS.

- [ ] **Step 4: Run the slow test once to confirm it works**

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
pytest tests/test_pipeline_smoke.py -v
```

Expected: PASS, possibly slow (1–10 min on first uncached run; seconds on cached). If it fails, the failure points to a regression introduced by the reorg — investigate before declaring the reorg done.

- [ ] **Step 5: Commit**

```bash
git -C /Users/jverboon/Analyses/epitope_pipeline add -A
git -C /Users/jverboon/Analyses/epitope_pipeline commit -m "Add end-to-end pipeline smoke test"
```

---

## Final state check

After Task 10:

```bash
cd /Users/jverboon/Analyses/epitope_pipeline
git log --oneline main..HEAD
```

Expected output (10 commits):
```
... Add end-to-end pipeline smoke test
... Update external callsites for new package layout
... Carve viz modules into epitope_pipeline.viz subpackage
... Carve compute modules into epitope_pipeline.compute subpackage
... Carve I/O modules into epitope_pipeline.io subpackage
... Remove Tamarind structure-prediction fallback
... Add pytest smoke test for module imports
... Move package sources into nested epitope_pipeline/ directory
... Add pyproject.toml for installable package
```

```bash
pytest tests/ -v
```
Expected: every test passes.

```bash
pip show epitope-pipeline | grep -E "^(Name|Version|Location)"
```
Expected: package is installed in editable mode and resolvable.

The branch is then ready to either:
1. Open a PR against `main` for review and merge.
2. Become the base for the perf branches (`perf/structure-context-cache`, `perf/local-blast-threads`, `perf/parallel-targets`, etc.) outlined in the original analysis.
