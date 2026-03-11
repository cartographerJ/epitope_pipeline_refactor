# Epitope Pipeline — Change History

## 2026-03-03: Distance measurement — bilayer surface reference

**File**: `spatial.py` (lines 108-116)

**Before** (surface_offset from most proximal ECD residue):
```python
# Re-zero distances to the membrane surface: the most proximal ECD
# residue becomes 0A. This measures "distance from the cell surface"
# rather than from the membrane center (~15-20A inside the bilayer).
extracellular_set = set(extracellular_residues)
ecd_dists_raw = [residue_distances[r] for r in extracellular_set
                 if r in residue_distances]
surface_offset = min(ecd_dists_raw) if ecd_dists_raw else 0.0
for r in residue_distances:
    residue_distances[r] = max(0.0, residue_distances[r] - surface_offset)
```

**After** (surface_offset from membrane half_thickness):
```python
# Re-zero distances to the membrane surface (top of bilayer).
# raw distance is from membrane center; subtract half_thickness to get
# distance from the bilayer surface. This is more robust than using the
# most proximal ECD residue, which can be artifactually buried in
# AlphaFold models that lack membrane context (e.g., ERBB2, EGFR).
extracellular_set = set(extracellular_residues)
surface_offset = membrane.membrane_half_thickness
for r in residue_distances:
    residue_distances[r] = max(0.0, residue_distances[r] - surface_offset)
```

**Why**: AlphaFold predicts without membrane context. For ERBB2 and EGFR, juxtamembrane ECD residues collapse to the membrane midplane in the AF model, giving `min(ecd_dists) ≈ 0A`. This inflated all distances by ~15A compared to well-behaved targets where `min(ecd_dists) ≈ half_thickness`. Using `half_thickness` (a known geometric constant) as the offset is deterministic and correct regardless of AF artifacts.

**Impact**:
- Well-behaved targets (NECTIN4, ITGAV, GPNMB): ≤4A shift, no change in results
- AF-artifact targets (ERBB2, EGFR): ~15A correction downward
- ERBB2 max distance: 95.5A → 80.5A (still above 80A threshold)

**To revert**: Replace `surface_offset = membrane.membrane_half_thickness` with the old 4-line block above.

---

## 2026-03-03: Membrane CGO — lipid sphere model

**File**: `export.py`, function `generate_membrane_cgo_pml()`

**Before**: Two flat TRIANGLE_FAN CGO disks (top/bottom of bilayer), warm golden color, 70% transparent.

**After**: Scattered CGO SPHERE primitives mimicking a lipid bilayer. ~350 lipid positions (fixed seed=42), each with headgroup spheres at bilayer surfaces + tail spheres filling hydrophobic core. 10A exclusion zone around TM helix. Mostly gray with sparse red/orange headgroup accents.

**To revert**: Replace the function body with the TRIANGLE_FAN disk version (search git-independent backups or session notes from earlier 2026-03-03).

---

## 2026-03-03: Patch trimming (conservation + specificity)

**Files**: `conservation.py`, `specificity.py`

**Before**: Sliding window finds mismatch-dense zone → reject entire patch.

**After**: Sliding window finds mismatch-dense zone → `identify_failing_residues()` → trim bad residues → `_recluster_survivors()` re-clusters → keep sub-patches ≥600A².

**Impact**: ERBB2 proximal patch 1 (129 res) had 17 cyno mismatches clustered at juxtamembrane. Old: entire patch rejected. New: trims 32 failing residues, salvages 83-res sub-patch (6,882 A²) + 9-res sub-patch (1,018 A²).

---

## 2026-03-03: CA coord consolidation + shared utils

**Files**: `utils.py` (new), `run.py`, `bispecific.py`, `spatial.py`, `surface.py`

**Before**: PDB parsed 4-5x per target. Helper functions duplicated across files.

**After**: `extract_ca_coords()`, `get_chain()`, `setup_logging()`, `empty_metric()` in `utils.py`. CA coords extracted once in orchestrator, passed via `ca_coords=` parameter to all downstream modules.
