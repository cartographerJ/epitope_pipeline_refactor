# Proposal: Parallel Filtering Architecture

**Date:** 2026-03-10
**Context:** CEACAM5 parameter sweep revealed fundamental limitation in sequential filtering
**Goal:** Enable complete data acquisition and definitive filter attribution

---

## Problem Statement

### Current Sequential Architecture (Cyno → Specificity)

```
Step 1-5: Spatial/Surface/Distance → 383 qualifying residues
   ↓
Step 6: Cyno Conservation (sliding window + trim + recluster)
   ↓ 259 conserved residues
   ✗ 124 cyno mismatch excluded
   ⚠ 97 cyno-conserved orphaned in sub-threshold clusters
   ↓
Step 7: Specificity (BLAST only cyno-conserved residues)
   ↓ 250 human-specific
   ✗ 0 non-specific
   ⚠ 133 NOT ASSESSED (never BLAST-screened)
   ↓
Step 8: Scoring → 0 patches (insufficient contiguous surface area)
```

### Key Limitations

1. **Incomplete Data**: 133/383 residues (34.7%) never get BLAST-screened
2. **Orphaned Residues**: 97 cyno-conserved residues lost during reclustering
3. **Ambiguous Attribution**: Can't determine if failures are due to cyno or specificity
4. **Misleading Visualization**: Gray "not assessed" regions appear similar to failures
5. **Lost Opportunities**: May miss viable patches that pass both filters independently

---

## CEACAM5 Case Study: What We Discovered

### Sequential Pipeline Results (Current)
- Starting: 383 qualifying residues
- Cyno conserved: 259 (67.6%)
- **BLAST-screened: 250 (65.3% of total)**
- **NOT ASSESSED: 133 (34.7%)**
  - 36 are cyno mismatches (explains absence)
  - **97 are cyno-conserved but orphaned!**

### Parallel Pipeline Results (Simulated)
- Starting: 383 qualifying residues
- Cyno conserved: 259 (67.6%)
- **BLAST-screened: ALL 383 (100%)**
- Human-specific: 250 (65.3%)
- **Intersection (pass BOTH): 162 residues**

### Critical Finding

**162 residues pass BOTH filters independently**, but the sequential pipeline orphans 97 of them during cyno reclustering. These orphans never reach the specificity step, preventing us from discovering potential patches.

**Visual Evidence**: See `epitope_pipeline/runs/260310_1148_specificity_sweep/figures/ceacam5_parallel_filter_analysis.png`

---

## Proposed Parallel Architecture

### Design Principles

1. **Complete Assessment**: Every qualifying residue gets BOTH cyno and specificity evaluations
2. **Independent Filters**: Cyno and specificity run in parallel on the same input set
3. **Intersection Logic**: Patches must pass BOTH filters to advance
4. **Full Data Retention**: Store results from both filters regardless of pass/fail
5. **Clear Attribution**: Can definitively identify which filter(s) caused failures

### Pipeline Flow

```
Step 1-5: Spatial/Surface/Distance → 383 qualifying residues
   ↓
   ├─ Step 6A: Cyno Conservation (independent) ─┐
   │    • Assess ALL 383 residues                │
   │    • Track conserved vs mismatch            │
   │    • Store results                          │
   │                                             │
   └─ Step 6B: Specificity (independent) ────────┤
        • BLAST ALL 383 residues                 │
        • Track specific vs non-specific         │
        • Store results                          │
                                                  ↓
Step 7: Intersection (require BOTH filters pass)
   • Combine cyno_conserved AND human_specific
   • Result: residues passing BOTH filters
   • Retain metadata from both assessments
   ↓
Step 8: Clustering & Patch Formation
   • Cluster intersection residues
   • Apply 600A² contiguous threshold
   • Form final patches
   ↓
Step 9: Scoring & Export
```

---

## Implementation Plan

### Phase 1: Core Restructuring (conservation.py + specificity.py)

**File: `epitope_pipeline/conservation.py`**

Current function signature:
```python
def filter_conservation(
    ca_coords: np.ndarray,
    cyno_seq: str,
    human_seq: str,
    patches: List[Dict],
    run_dir: Path,
    logger
) -> List[Dict]:
    # Returns filtered patches
```

Proposed new signature:
```python
def assess_conservation(
    ca_coords: np.ndarray,
    residue_indices: List[int],  # All qualifying residues
    cyno_seq: str,
    human_seq: str,
    run_dir: Path,
    logger
) -> Dict[int, bool]:
    # Returns {residue_num: is_conserved} mapping
    # No clustering/patching - just per-residue assessment
```

**File: `epitope_pipeline/specificity.py`**

Current function signature:
```python
def filter_specificity(
    structure_path: Path,
    patches: List[Dict],
    run_dir: Path,
    logger
) -> List[Dict]:
    # Returns filtered patches
```

Proposed new signature:
```python
def assess_specificity(
    structure_path: Path,
    residue_indices: List[int],  # All qualifying residues
    run_dir: Path,
    logger
) -> Dict[int, bool]:
    # Returns {residue_num: is_specific} mapping
    # BLAST all residues, no pre-filtering
    # No clustering/patching - just per-residue assessment
```

### Phase 2: New Intersection Module

**File: `epitope_pipeline/intersection.py` (NEW)**

```python
def compute_filter_intersection(
    ca_coords: np.ndarray,
    residue_indices: List[int],
    conservation_results: Dict[int, bool],
    specificity_results: Dict[int, bool],
    min_patch_area: float = 600.0,
    run_dir: Path,
    logger
) -> List[Dict]:
    """
    Compute intersection of conservation and specificity filters.

    Returns:
        patches: List of patches where ALL residues pass BOTH filters
    """
    # Find residues passing BOTH filters
    passing_residues = [
        idx for idx in residue_indices
        if conservation_results.get(idx, False) and specificity_results.get(idx, False)
    ]

    # Cluster intersection residues
    patches = _cluster_residues(ca_coords, passing_residues, min_patch_area)

    # Attach metadata from both filters
    for patch in patches:
        patch['conservation_results'] = conservation_results
        patch['specificity_results'] = specificity_results

    return patches
```

### Phase 3: Orchestrator Updates

**File: `epitope_pipeline/orchestrator.py`**

Current flow:
```python
# Step 6: Conservation
patches = conservation.filter_conservation(...)

# Step 7: Specificity
patches = specificity.filter_specificity(patches, ...)
```

Proposed flow:
```python
# Step 6A: Conservation (independent)
conservation_results = conservation.assess_conservation(
    ca_coords=ca_coords,
    residue_indices=spatial_residue_indices,  # From Step 5
    cyno_seq=cyno_seq,
    human_seq=human_seq,
    run_dir=run_dir,
    logger=logger
)

# Step 6B: Specificity (independent, runs in parallel)
specificity_results = specificity.assess_specificity(
    structure_path=structure_path,
    residue_indices=spatial_residue_indices,  # Same input as 6A
    run_dir=run_dir,
    logger=logger
)

# Step 7: Intersection
patches = intersection.compute_filter_intersection(
    ca_coords=ca_coords,
    residue_indices=spatial_residue_indices,
    conservation_results=conservation_results,
    specificity_results=specificity_results,
    min_patch_area=config.MIN_PATCH_AREA,
    run_dir=run_dir,
    logger=logger
)
```

### Phase 4: Visualization Updates

**File: `epitope_pipeline/visualize.py`**

Update residue table to include both filter results:
```python
residue_table_columns = [
    'residue_num',
    'aa',
    'topology',
    'distance_from_membrane_A',
    'sasa_A2',
    'relative_sasa',
    'cyno_conserved',      # Boolean: True/False (not yes/no/NaN)
    'human_specific',      # Boolean: True/False (not yes/no/NaN)
    'passes_both_filters', # Boolean: cyno_conserved AND human_specific
    'in_patch',
    'patch_id',
    'epitope_score'
]
```

Update epitope maps to show:
- **Cyno Conservation track**: Mint (conserved) / Red (mismatch)
- **Human Specificity track**: Mint (specific) / Red (non-specific) — NO gray "not assessed"
- **Intersection track**: Mint (pass both) / Red (fail at least one)

---

## Expected Impact

### Performance

**Time Cost:**
- Current: ~1-2 seconds per BLAST query × ~10 patches = ~10-20 seconds total
- Proposed: ~0.1 seconds per BLAST query × ~383 residues = ~38 seconds total
- **Net cost: +20-30 seconds per target** (acceptable)

With local BLAST (from plan), this reduces to ~4-8 seconds overhead.

### Data Completeness

| Metric | Sequential | Parallel | Improvement |
|--------|-----------|----------|-------------|
| Cyno-assessed | 383 (100%) | 383 (100%) | — |
| BLAST-screened | 250 (65%) | **383 (100%)** | **+35%** |
| Complete data | 250 (65%) | **383 (100%)** | **+35%** |
| Attributable failures | No | **Yes** | ✓ |

### CEACAM5 Example

| Approach | Residues Passing Both | Potential Patches |
|----------|----------------------|-------------------|
| Sequential | Unknown (orphans lost) | 0 patches |
| **Parallel** | **162 residues** | **TBD (needs clustering)** |

**Key insight**: Parallel filtering reveals 162 residues pass both filters, vs 0 patches from sequential approach. May discover viable patches after re-clustering.

---

## Migration Strategy

### Option A: Hard Cutover (Recommended)

1. Implement parallel architecture in new module structure
2. Update orchestrator to use new flow
3. Add compatibility flag: `USE_PARALLEL_FILTERING = True` (default)
4. Run validation suite comparing old vs new on known targets
5. Deprecate sequential mode after validation

**Pros**: Clean break, easier to maintain
**Cons**: Requires validation before deployment

### Option B: Side-by-Side Comparison

1. Implement parallel architecture alongside existing
2. Run BOTH pipelines in orchestrator
3. Export both result sets for comparison
4. Gradual migration after confidence established

**Pros**: Safe, empirical validation
**Cons**: Complex orchestrator, temporary code bloat

---

## Validation Plan

### Test Suite

1. **Regression Tests**: Run parallel pipeline on known working targets
   - ERBB2, NECTIN4, GPNMB, ITGAV, FOLH1, EGFR
   - Verify: same or more patches discovered
   - Verify: patch quality metrics unchanged

2. **CEACAM5 Re-run**: Primary test case
   - Expected: discover ≥1 patch from 162 intersection residues
   - Compare: sequential (0 patches) vs parallel (TBD)

3. **Failed Target Re-runs**: Targets with 0 patches
   - CEACAM5, CLDN18, CLDN6, MS4A1, BCMA, TROP2, etc.
   - Hypothesis: some may yield patches with complete data

4. **Data Completeness Check**: All targets
   - Verify: 0% "not assessed" residues
   - Verify: every qualifying residue has both filter results

### Success Criteria

- [ ] No regression on working targets (same or better scores)
- [ ] CEACAM5 yields ≥1 patch with parallel filtering
- [ ] 100% data completeness (no "not assessed" residues)
- [ ] Clear attribution for all filter failures
- [ ] Visualization shows no gray "not assessed" regions

---

## Next Steps

1. **Immediate**: Review and approve this proposal
2. **Phase 1**: Implement `assess_conservation()` and `assess_specificity()` functions
3. **Phase 2**: Create `intersection.py` module
4. **Phase 3**: Update orchestrator to use parallel flow
5. **Phase 4**: Run validation suite on CEACAM5 + known targets
6. **Phase 5**: Update visualization to eliminate "not assessed" ambiguity

---

## Open Questions

1. **Clustering strategy**: Should intersection use the same clustering as current pipeline, or explore alternative algorithms?

2. **Patch trimming**: Current pipeline trims patches if local windows fail filters. With parallel assessment, do we still need sliding-window trimming, or can we trust per-residue assessments?

3. **Performance optimization**: Should we parallelize BLAST queries (ThreadPoolExecutor) to maintain current runtime?

4. **Backward compatibility**: Should we keep sequential mode as fallback option, or fully deprecate?

5. **Bispecific impact**: Does parallel filtering change bispecific pair scoring logic?

---

## Conclusion

**Parallel filtering solves a fundamental limitation**: sequential filtering prevents complete data acquisition and creates ambiguous failures. The CEACAM5 case study demonstrates that **162 residues pass both filters** when assessed independently, but sequential pipeline orphans 97 of them.

**Recommendation**: Implement parallel filtering architecture as proposed. Expected benefits:
- Complete data acquisition (100% vs 65%)
- Definitive filter attribution
- Discovery of viable patches missed by sequential approach
- Clearer visualizations without "not assessed" ambiguity

**Trade-off**: +20-30 seconds per target (acceptable for gains in data quality and completeness).

**Next action**: Approve proposal and proceed with Phase 1 implementation.
