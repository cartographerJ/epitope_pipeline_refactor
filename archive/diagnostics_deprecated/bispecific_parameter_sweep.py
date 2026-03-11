"""
Bispecific parameter sweep for CEACAM5:ERBB2 pair.

Tests the bispecific pair at multiple cyno conservation and human specificity
threshold values to determine optimal defaults. All results contained in a
single dated run folder.

This is a temporary analysis script - results inform config.py defaults but
the script itself doesn't modify any production code.
"""

import logging
import os
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from epitope_pipeline import config
from epitope_pipeline.bispecific import run_bispecific
from epitope_pipeline.utils import setup_logging

# Cartography color palette
CARTO_GRAY = '#D3D3D3'
CARTO_MINT = '#E0F3DB'
CARTO_GREEN = '#6BC291'
CARTO_TEAL = '#18B5CB'
CARTO_BLUE = '#2E95D2'
CARTO_PURPLE = '#28154C'


def run_bispecific_threshold_sweep(pair_name, run_dir):
    """
    Run parameter sweep on bispecific pair with varying thresholds.

    Two separate sweeps:
    1. Cyno conservation: 10%, 15%, 20%, 25%, 30% (specificity fixed at 15%)
    2. Human specificity: 10%, 15%, 20%, 25%, 30% (conservation fixed at 15%)

    Args:
        pair_name: Pair in format "GENE_A:GENE_B"
        run_dir: Output directory for all results (single parent folder)
    """
    logger = logging.getLogger(__name__)

    # Parse pair
    gene_a, gene_b = pair_name.split(':')

    # Create subdirectories within the single run folder
    cyno_sweep_dir = run_dir / "cyno_sweep"
    spec_sweep_dir = run_dir / "specificity_sweep"
    summary_dir = run_dir / "summary"
    cyno_sweep_dir.mkdir(parents=True, exist_ok=True)
    spec_sweep_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    thresholds = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]

    # Results storage: {threshold: {pair_score, validity, orientation_scores}}
    cyno_results = {}
    spec_results = {}

    logger.info("=" * 80)
    logger.info(f"BISPECIFIC PARAMETER SWEEP: {pair_name}")
    logger.info("=" * 80)
    logger.info(f"Output directory: {run_dir}")
    logger.info("=" * 80)

    # Sweep 1: Vary cyno conservation (specificity fixed at 15%)
    logger.info("\n" + "=" * 80)
    logger.info("SWEEP 1: Cyno Conservation Thresholds")
    logger.info("(Human Specificity fixed at 15%)")
    logger.info("=" * 80)

    for threshold in thresholds:
        logger.info("\n" + "-" * 80)
        logger.info(f"Testing cyno threshold: {threshold}% (specificity=15%)")
        logger.info("-" * 80 + "\n")

        # Temporarily modify config
        original_cyno = config.MAX_CYNO_MISMATCH_PERCENT
        config.MAX_CYNO_MISMATCH_PERCENT = threshold

        try:
            # Run bispecific analysis
            # Parse pair name (e.g., "CEACAM5:ERBB2") into tuple
            gene_a, gene_b = pair_name.split(":")

            # Create run name to nest inside parameter sweep folder
            custom_run_name = f"{run_dir.name}/cyno_sweep/cyno_{int(threshold)}pct"
            output = run_bispecific([(gene_a, gene_b)], run_name=custom_run_name)
            pair_results = output.get("pair_results", [])

            if pair_results and len(pair_results) > 0:
                result = pair_results[0]

                # Extract orientation scores
                ab_score = result.orientation_ab.orientation_score if result.orientation_ab else 0.0
                ba_score = result.orientation_ba.orientation_score if result.orientation_ba else 0.0

                # Extract patch counts (gene_a distal = ab orientation, gene_a proximal = ba orientation)
                gene_a_distal = len(result.orientation_ab.distal_zone.specificity_result.specific_patches) if result.orientation_ab and result.orientation_ab.distal_zone.specificity_result else 0
                gene_b_proximal = len(result.orientation_ab.proximal_zone.specificity_result.specific_patches) if result.orientation_ab and result.orientation_ab.proximal_zone.specificity_result else 0
                gene_a_proximal = len(result.orientation_ba.proximal_zone.specificity_result.specific_patches) if result.orientation_ba and result.orientation_ba.proximal_zone.specificity_result else 0
                gene_b_distal = len(result.orientation_ba.distal_zone.specificity_result.specific_patches) if result.orientation_ba and result.orientation_ba.distal_zone.specificity_result else 0

                cyno_results[threshold] = {
                    'pair_score': result.final_pair_score,
                    'both_valid': result.both_valid,
                    'best_orientation': result.best_orientation.label if result.best_orientation else None,
                    'orientation_ab_score': ab_score,
                    'orientation_ba_score': ba_score,
                    'gene_a_distal_patches': gene_a_distal,
                    'gene_a_proximal_patches': gene_a_proximal,
                    'gene_b_distal_patches': gene_b_distal,
                    'gene_b_proximal_patches': gene_b_proximal,
                }

                logger.info(f"  Pair score: {result.final_pair_score:.3f}")
                logger.info(f"  Both valid: {result.both_valid}")
                logger.info(f"  Best orientation: {result.best_orientation.label if result.best_orientation else 'None'}")
            else:
                cyno_results[threshold] = {
                    'pair_score': 0.0,
                    'both_valid': False,
                    'best_orientation': None,
                    'orientation_ab_score': 0.0,
                    'orientation_ba_score': 0.0,
                    'gene_a_distal_patches': 0,
                    'gene_a_proximal_patches': 0,
                    'gene_b_distal_patches': 0,
                    'gene_b_proximal_patches': 0,
                }
                logger.info("  No valid pair results")

        finally:
            # Restore original config
            config.MAX_CYNO_MISMATCH_PERCENT = original_cyno

    # Sweep 2: Vary human specificity (cyno fixed at 15%)
    logger.info("\n" + "=" * 80)
    logger.info("SWEEP 2: Human Specificity Thresholds")
    logger.info("(Cyno Conservation fixed at 15%)")
    logger.info("=" * 80)

    for threshold in thresholds:
        logger.info("\n" + "-" * 80)
        logger.info(f"Testing specificity threshold: {threshold}% (cyno=15%)")
        logger.info("-" * 80 + "\n")

        # Temporarily modify config
        original_spec = config.MAX_NONSPECIFIC_PERCENT
        config.MAX_NONSPECIFIC_PERCENT = threshold

        try:
            # Run bispecific analysis
            # Parse pair name (e.g., "CEACAM5:ERBB2") into tuple
            gene_a, gene_b = pair_name.split(":")

            # Create run name to nest inside parameter sweep folder
            custom_run_name = f"{run_dir.name}/specificity_sweep/spec_{int(threshold)}pct"
            output = run_bispecific([(gene_a, gene_b)], run_name=custom_run_name)
            pair_results = output.get("pair_results", [])

            if pair_results and len(pair_results) > 0:
                result = pair_results[0]

                # Extract orientation scores
                ab_score = result.orientation_ab.orientation_score if result.orientation_ab else 0.0
                ba_score = result.orientation_ba.orientation_score if result.orientation_ba else 0.0

                # Extract patch counts (gene_a distal = ab orientation, gene_a proximal = ba orientation)
                gene_a_distal = len(result.orientation_ab.distal_zone.specificity_result.specific_patches) if result.orientation_ab and result.orientation_ab.distal_zone.specificity_result else 0
                gene_b_proximal = len(result.orientation_ab.proximal_zone.specificity_result.specific_patches) if result.orientation_ab and result.orientation_ab.proximal_zone.specificity_result else 0
                gene_a_proximal = len(result.orientation_ba.proximal_zone.specificity_result.specific_patches) if result.orientation_ba and result.orientation_ba.proximal_zone.specificity_result else 0
                gene_b_distal = len(result.orientation_ba.distal_zone.specificity_result.specific_patches) if result.orientation_ba and result.orientation_ba.distal_zone.specificity_result else 0

                spec_results[threshold] = {
                    'pair_score': result.final_pair_score,
                    'both_valid': result.both_valid,
                    'best_orientation': result.best_orientation.label if result.best_orientation else None,
                    'orientation_ab_score': ab_score,
                    'orientation_ba_score': ba_score,
                    'gene_a_distal_patches': gene_a_distal,
                    'gene_a_proximal_patches': gene_a_proximal,
                    'gene_b_distal_patches': gene_b_distal,
                    'gene_b_proximal_patches': gene_b_proximal,
                }

                logger.info(f"  Pair score: {result.final_pair_score:.3f}")
                logger.info(f"  Both valid: {result.both_valid}")
                logger.info(f"  Best orientation: {result.best_orientation.label if result.best_orientation else 'None'}")
            else:
                spec_results[threshold] = {
                    'pair_score': 0.0,
                    'both_valid': False,
                    'best_orientation': None,
                    'orientation_ab_score': 0.0,
                    'orientation_ba_score': 0.0,
                    'gene_a_distal_patches': 0,
                    'gene_a_proximal_patches': 0,
                    'gene_b_distal_patches': 0,
                    'gene_b_proximal_patches': 0,
                }
                logger.info("  No valid pair results")

        finally:
            # Restore original config
            config.MAX_NONSPECIFIC_PERCENT = original_spec

    # Generate summary figures
    logger.info("\n" + "=" * 80)
    logger.info("Generating summary figures...")
    logger.info("=" * 80 + "\n")

    _plot_cyno_sweep_summary(cyno_results, thresholds, pair_name, summary_dir)
    _plot_spec_sweep_summary(spec_results, thresholds, pair_name, summary_dir)
    _write_summary_report(cyno_results, spec_results, thresholds, pair_name, summary_dir)

    logger.info(f"\nParameter sweep complete!")
    logger.info(f"Results in: {run_dir}")


def _plot_cyno_sweep_summary(results, thresholds, pair_name, output_dir):
    """Generate summary plots for cyno conservation sweep."""
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    fig.suptitle(f'Cyno Conservation Threshold Sweep\n{pair_name} Bispecific Pair\n(Human Specificity fixed at 15%)',
                 fontsize=16, fontweight='bold')

    # Plot 1: Pair score
    ax1 = axes[0]
    pair_scores = [results[t]['pair_score'] for t in thresholds]
    ax1.plot(thresholds, pair_scores, 'o-', linewidth=2, markersize=8,
            color=CARTO_PURPLE)
    ax1.set_ylabel('Bispecific Pair Score', fontsize=12, fontweight='bold')
    ax1.set_ylim(bottom=-0.05, top=1.2)
    ax1.grid(True, alpha=0.3)
    ax1.axvline(x=15.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)
    ax1.axhline(y=1.0, color='gray', linestyle=':', linewidth=1, alpha=0.3)

    # Plot 2: Orientation scores
    ax2 = axes[1]
    orientation_ab_scores = [results[t]['orientation_ab_score'] for t in thresholds]
    orientation_ba_scores = [results[t]['orientation_ba_score'] for t in thresholds]
    ax2.plot(thresholds, orientation_ab_scores, 'o-', linewidth=2, markersize=8,
            color=CARTO_TEAL, label='Orientation AB')
    ax2.plot(thresholds, orientation_ba_scores, 's-', linewidth=2, markersize=8,
            color=CARTO_BLUE, label='Orientation BA')
    ax2.set_ylabel('Orientation Scores', fontsize=12, fontweight='bold')
    ax2.set_ylim(bottom=-0.05, top=1.05)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=11)
    ax2.axvline(x=15.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)

    # Plot 3: Patch counts
    ax3 = axes[2]
    gene_a, gene_b = pair_name.split(':')
    a_distal = [results[t]['gene_a_distal_patches'] for t in thresholds]
    a_proximal = [results[t]['gene_a_proximal_patches'] for t in thresholds]
    b_distal = [results[t]['gene_b_distal_patches'] for t in thresholds]
    b_proximal = [results[t]['gene_b_proximal_patches'] for t in thresholds]

    width = 1.0
    x = np.array(thresholds)
    ax3.bar(x - width*1.5, a_distal, width, label=f'{gene_a} distal', color=CARTO_TEAL, alpha=0.7)
    ax3.bar(x - width*0.5, a_proximal, width, label=f'{gene_a} proximal', color=CARTO_TEAL, alpha=0.4)
    ax3.bar(x + width*0.5, b_distal, width, label=f'{gene_b} distal', color=CARTO_PURPLE, alpha=0.7)
    ax3.bar(x + width*1.5, b_proximal, width, label=f'{gene_b} proximal', color=CARTO_PURPLE, alpha=0.4)

    ax3.set_xlabel('Max Cyno Mismatch Percent (%)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Number of Patches', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')
    ax3.legend(fontsize=10, ncol=2)
    ax3.axvline(x=15.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5,
                label='Current default')

    plt.tight_layout()
    output_path = output_dir / 'cyno_threshold_sweep.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path.name}")


def _plot_spec_sweep_summary(results, thresholds, pair_name, output_dir):
    """Generate summary plots for human specificity sweep."""
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    fig.suptitle(f'Human Specificity Threshold Sweep\n{pair_name} Bispecific Pair\n(Cyno Conservation fixed at 15%)',
                 fontsize=16, fontweight='bold')

    # Plot 1: Pair score
    ax1 = axes[0]
    pair_scores = [results[t]['pair_score'] for t in thresholds]
    ax1.plot(thresholds, pair_scores, 'o-', linewidth=2, markersize=8,
            color=CARTO_PURPLE)
    ax1.set_ylabel('Bispecific Pair Score', fontsize=12, fontweight='bold')
    ax1.set_ylim(bottom=-0.05, top=1.2)
    ax1.grid(True, alpha=0.3)
    ax1.axvline(x=15.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)
    ax1.axhline(y=1.0, color='gray', linestyle=':', linewidth=1, alpha=0.3)

    # Plot 2: Orientation scores
    ax2 = axes[1]
    orientation_ab_scores = [results[t]['orientation_ab_score'] for t in thresholds]
    orientation_ba_scores = [results[t]['orientation_ba_score'] for t in thresholds]
    ax2.plot(thresholds, orientation_ab_scores, 'o-', linewidth=2, markersize=8,
            color=CARTO_TEAL, label='Orientation AB')
    ax2.plot(thresholds, orientation_ba_scores, 's-', linewidth=2, markersize=8,
            color=CARTO_BLUE, label='Orientation BA')
    ax2.set_ylabel('Orientation Scores', fontsize=12, fontweight='bold')
    ax2.set_ylim(bottom=-0.05, top=1.05)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=11)
    ax2.axvline(x=15.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)

    # Plot 3: Patch counts
    ax3 = axes[2]
    gene_a, gene_b = pair_name.split(':')
    a_distal = [results[t]['gene_a_distal_patches'] for t in thresholds]
    a_proximal = [results[t]['gene_a_proximal_patches'] for t in thresholds]
    b_distal = [results[t]['gene_b_distal_patches'] for t in thresholds]
    b_proximal = [results[t]['gene_b_proximal_patches'] for t in thresholds]

    width = 1.0
    x = np.array(thresholds)
    ax3.bar(x - width*1.5, a_distal, width, label=f'{gene_a} distal', color=CARTO_TEAL, alpha=0.7)
    ax3.bar(x - width*0.5, a_proximal, width, label=f'{gene_a} proximal', color=CARTO_TEAL, alpha=0.4)
    ax3.bar(x + width*0.5, b_distal, width, label=f'{gene_b} distal', color=CARTO_PURPLE, alpha=0.7)
    ax3.bar(x + width*1.5, b_proximal, width, label=f'{gene_b} proximal', color=CARTO_PURPLE, alpha=0.4)

    ax3.set_xlabel('Max Non-Specific Percent (%)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Number of Patches', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')
    ax3.legend(fontsize=10, ncol=2)
    ax3.axvline(x=15.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)

    plt.tight_layout()
    output_path = output_dir / 'specificity_threshold_sweep.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path.name}")


def _write_summary_report(cyno_results, spec_results, thresholds, pair_name, output_dir):
    """Write text summary of parameter sweep results."""
    output_path = output_dir / 'parameter_sweep_summary.txt'

    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(f"BISPECIFIC PARAMETER SWEEP: {pair_name}\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Test Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"Thresholds tested: {', '.join(f'{t}%' for t in thresholds)}\n\n")

        # Cyno sweep summary
        f.write("=" * 80 + "\n")
        f.write("CYNO CONSERVATION THRESHOLD SWEEP\n")
        f.write("(Human Specificity fixed at 15%)\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"  Threshold | Pair Score | Both Valid | Best Orientation\n")
        f.write(f"  ----------|------------|------------|-----------------\n")
        for threshold in thresholds:
            data = cyno_results[threshold]
            marker = " <-- current default" if threshold == 15.0 else ""
            validity_str = "Yes" if data['both_valid'] else "No"
            f.write(f"  {threshold:>6.0f}%   | {data['pair_score']:>10.3f} | {validity_str:<10} | "
                   f"{data['best_orientation']}{marker}\n")
        f.write("\n")

        # Specificity sweep summary
        f.write("=" * 80 + "\n")
        f.write("HUMAN SPECIFICITY THRESHOLD SWEEP\n")
        f.write("(Cyno Conservation fixed at 15%)\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"  Threshold | Pair Score | Both Valid | Best Orientation\n")
        f.write(f"  ----------|------------|------------|-----------------\n")
        for threshold in thresholds:
            data = spec_results[threshold]
            marker = " <-- current default" if threshold == 15.0 else ""
            validity_str = "Yes" if data['both_valid'] else "No"
            f.write(f"  {threshold:>6.0f}%   | {data['pair_score']:>10.3f} | {validity_str:<10} | "
                   f"{data['best_orientation']}{marker}\n")
        f.write("\n")

        # Recommendations
        f.write("=" * 80 + "\n")
        f.write("RECOMMENDATIONS\n")
        f.write("=" * 80 + "\n\n")

        f.write("Based on this bispecific parameter sweep:\n\n")

        # Analyze results
        cyno_15 = cyno_results[15.0]['pair_score']
        spec_15 = spec_results[15.0]['pair_score']

        f.write(f"Cyno Conservation (15% default):\n")
        f.write(f"  - Pair score at 15%: {cyno_15:.3f}\n")
        f.write(f"  - Both valid: {'Yes' if cyno_results[15.0]['both_valid'] else 'No'}\n")
        if cyno_15 == 0.0:
            # Find first working threshold
            for t in thresholds:
                if cyno_results[t]['pair_score'] > 0.0:
                    f.write(f"  - First working threshold: {t}% (score={cyno_results[t]['pair_score']:.3f})\n")
                    break

        f.write(f"\nHuman Specificity (15% default):\n")
        f.write(f"  - Pair score at 15%: {spec_15:.3f}\n")
        f.write(f"  - Both valid: {'Yes' if spec_results[15.0]['both_valid'] else 'No'}\n")
        if spec_15 == 0.0:
            # Find first working threshold
            for t in thresholds:
                if spec_results[t]['pair_score'] > 0.0:
                    f.write(f"  - First working threshold: {t}% (score={spec_results[t]['pair_score']:.3f})\n")
                    break

    print(f"  Saved: {output_path.name}")


def main():
    """Run bispecific parameter sweep on CEACAM5:ERBB2 pair."""
    # Setup
    timestamp = datetime.now().strftime("%y%m%d_%H%M")
    run_name = f"{timestamp}_param_sweep_ceacam5_erbb2"
    run_dir = Path(__file__).parent.parent / "runs" / run_name
    run_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    log_file = run_dir / "parameter_sweep.log"
    setup_logging(verbose=True)

    # Add file handler
    logger = logging.getLogger()
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    ))
    logger.addHandler(file_handler)

    # Run sweep
    pair_name = 'CEACAM5:ERBB2'
    run_bispecific_threshold_sweep(pair_name, run_dir)


if __name__ == '__main__':
    main()
