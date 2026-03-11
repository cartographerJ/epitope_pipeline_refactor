# Archive

Historical files preserved for reference. Not used in current production pipeline.

## Contents

### HISTORY.md
Early development history log (pre-March 2026). Superseded by memory/epitope_pipeline_notes.md for detailed session context.

### debug_sliding_window.py
Debug script for troubleshooting sliding window trim-and-recluster approach. Replaced by whole-patch evaluation in March 2026 refactor.

### test_specificity_params.py
Early testing script for specificity threshold tuning. Replaced by comprehensive 2D parameter sweeps.

### proposals/
**PARALLEL_FILTERING_PROPOSAL.md** - Initial proposal for parallel conservation and specificity filtering. Implemented differently using whole-patch evaluation instead of proposed sliding window enhancements.

### diagnostics_deprecated/
Deprecated diagnostic scripts replaced by production versions:

- **bispecific_parameter_sweep.py** - Old 1D parameter sweep (replaced by bispecific_parameter_sweep_2d.py)
- **parameter_sweep.py** - Single-target parameter sweep (replaced by bispecific sweeps)
- **generate_ceacam5_diagnostics.py** - CEACAM5-specific diagnostic (temporary, used for initial validation)
- **generate_ceacam5_patch_diagnostics.py** - CEACAM5 patch analysis (temporary, used for parameter sweep validation)
- **plot_specificity_sweep.py** - Old plotting script (replaced by integrated heatmaps in 2D sweep)
- **test_specificity_sweep.py** - Old test script (replaced by comprehensive validation)
- **visualize_parallel_filters.py** - Parallel filter visualization (deprecated with whole-patch refactor)

## Why These Were Archived

**Sliding window approach** (debug_sliding_window.py, PARALLEL_FILTERING_PROPOSAL.md):
- Complex: ~350 lines across 4 recursive functions
- Biologically inaccurate: trimming changed epitope surfaces
- Hard to interpret: "trimmed 149→19 residues"
- Replaced by: Whole-patch evaluation (150 lines, percentage-based thresholds)

**Old parameter sweeps** (bispecific_parameter_sweep.py, parameter_sweep.py):
- Only tested one threshold at a time (1D sweeps)
- Didn't reveal full parameter space behavior
- Replaced by: 2D matrix sweeps testing all (cyno, spec) combinations

**Temporary diagnostics** (generate_ceacam5_*.py):
- Target-specific scripts used for initial validation
- Not generalizable to other targets
- Replaced by: Generic specificity_diagnostics.py module

## Recovery

If you need to reference old approaches or recover archived scripts:
1. Check this directory for the original files
2. See memory/epitope_pipeline_notes.md for implementation context
3. Git history preserves all changes: `git log --all -- <filename>`
