#!/usr/bin/env python3
"""Run bispecific analysis on a batch of target pairs."""

import logging
import sys
from datetime import datetime
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from epitope_pipeline.utils import setup_logging
from epitope_pipeline.bispecific import run_bispecific


def main():
    """Run bispecific analysis on multiple target pairs."""
    timestamp = datetime.now().strftime("%y%m%d_%H%M")

    # Define all target pairs
    pairs = [
        ("ITGB6", "ADAM8"),
        ("MSLN", "ADAM8"),
        ("MSLN", "MELTF"),
        ("CEACAM5", "CLDN18.2"),  # Isoform 2
        ("ROR1", "EGFR"),
        ("ERBB2", "EGFR"),
        ("MUCL3", "ABCC11"),
        ("LY6K", "LY6G6D"),
    ]

    # Create run name from pairs
    run_name = f"{timestamp}_bispecific_batch"

    # Setup logging
    setup_logging(verbose=True)
    logger = logging.getLogger(__name__)

    logger.info("="*80)
    logger.info(f"BISPECIFIC BATCH RUN: {len(pairs)} pairs")
    logger.info("="*80)
    for gene_a, gene_b in pairs:
        logger.info(f"  - {gene_a} x {gene_b}")
    logger.info("="*80)

    # Run all pairs in a single bispecific analysis
    # This creates one unified output directory with all pairs together
    run_bispecific(
        pairs,
        run_name=run_name,
        cyno_mismatch_percent=15.0,
        nonspecific_percent=15.0,
        verbose=True,
    )

    logger.info("")
    logger.info("="*80)
    logger.info("BATCH COMPLETE")
    logger.info("="*80)
    logger.info(f"Results in: epitope_pipeline/runs/{run_name}/")
    logger.info("  - bispecific_pairs.csv/xlsx - Summary table")
    logger.info("  - figures/bispecific_pair_summary.png - All pairs ranked")
    logger.info("  - pymol/ - PyMOL sessions for all pairs")
    logger.info("="*80)


if __name__ == '__main__':
    main()
