# Diagnostics Scripts

Production-ready diagnostic and validation tools for the epitope pipeline.

## Current Tools

### bispecific_parameter_sweep_2d.py
**Purpose:** Comprehensive 2D parameter sweep for validating bispecific pair filtering logic

**Usage:**
```bash
python -m epitope_pipeline.diagnostics.bispecific_parameter_sweep_2d
```

**What it does:**
- Tests all combinations of cyno conservation (0-100%) and specificity (0-100%) thresholds
- Default: 6×6 matrix (0%, 20%, 40%, 60%, 80%, 100%)
- Generates pair score heatmaps and patch count matrices
- All results contained in single timestamped run folder

**Outputs:**
- `runs/YYMMDD_HHMM_param_sweep_2d_<pair>/summary/` - heatmaps and text summaries
- `runs/YYMMDD_HHMM_param_sweep_2d_<pair>/threshold_combinations/` - individual bispecific runs

### run_adam8_itgb6_sweep.py
**Purpose:** Example runner script for ADAM8:ITGB6 parameter sweep

**Usage:**
```bash
python -m epitope_pipeline.diagnostics.run_adam8_itgb6_sweep
```

Use this as a template for sweeping other bispecific pairs.

### specificity_diagnostics.py
**Purpose:** Per-residue specificity analysis and visualization

**Usage:**
```python
from epitope_pipeline.diagnostics.specificity_diagnostics import generate_specificity_report

generate_specificity_report(
    target_gene="CEACAM5",
    output_dir="diagnostics_output/"
)
```

**What it does:**
- Detailed BLAST analysis per patch
- Per-residue specificity visualization
- Off-target hit tables with identity scores

## Archive

Deprecated scripts and old approaches are in `../archive/diagnostics_deprecated/`:
- `bispecific_parameter_sweep.py` - old 1D sweep (replaced by 2D version)
- `parameter_sweep.py` - old single-target sweep
- `generate_ceacam5_diagnostics.py` - CEACAM5-specific temporary script
- `generate_ceacam5_patch_diagnostics.py` - CEACAM5 patch-specific temporary script
- `plot_specificity_sweep.py` - old plotting script
- `test_specificity_sweep.py` - old test script
- `visualize_parallel_filters.py` - old parallel filter visualization

## Adding New Diagnostics

When creating new diagnostic scripts:
1. Follow the naming pattern: `<purpose>_<target/analysis>.py`
2. Add a `main()` function for CLI execution
3. Make it runnable as a module: `python -m epitope_pipeline.diagnostics.<script_name>`
4. Document usage in this README
5. Output to timestamped subdirectories in `runs/`
