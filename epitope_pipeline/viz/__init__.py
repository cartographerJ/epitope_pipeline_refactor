"""Viz subpackage: matplotlib figure generation.

All figures for a pipeline run are written under ``runs/<run_name>/Figures/``:

  ``runs/<run_name>/Figures/``
  ├── ``{gene}_epitope_map.png``       — per-target 6-track linear map
  ├── ``scoring_summary.png``          — multi-target druggability bar chart (only
  │                                      when more than one target was run)
  └── ``BLAST/``
      └── ``{gene}_blast_offtargets.png`` — per-target paralog dot plot

Cross-target metrics for batch runs are at
``runs/<run_name>/Supplementary Files/summary_metrics.csv``
(one row per input target, pandas-readable).
"""
