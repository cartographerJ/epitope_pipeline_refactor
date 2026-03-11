"""
Parameter sweep for cyno conservation and human specificity thresholds.

Tests CEACAM5 and ERBB2 (HER2) at multiple threshold values to determine
optimal defaults. Generates diagnostic figures and summary plots.

This is a temporary analysis script - results inform config.py defaults but
the script itself doesn't modify any production code.
"""

import logging
import shutil
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from epitope_pipeline import config
from epitope_pipeline.run import run_pipeline
from epitope_pipeline.utils import setup_logging

# Cartography color palette
CARTO_GRAY = '#D3D3D3'
CARTO_MINT = '#E0F3DB'
CARTO_GREEN = '#6BC291'
CARTO_TEAL = '#18B5CB'
CARTO_BLUE = '#2E95D2'
CARTO_PURPLE = '#28154C'


def run_threshold_sweep(targets, run_dir):
    """
    Run parameter sweep on conservation and specificity thresholds.

    Two separate sweeps:
    1. Cyno conservation: 10%, 15%, 20%, 25%, 30% (specificity fixed at 15%)
    2. Human specificity: 10%, 15%, 20%, 25%, 30% (conservation fixed at 15%)

    Args:
        targets: List of gene names to test
        run_dir: Output directory for all results
    """
    logger = logging.getLogger(__name__)

    # Create summary directory (individual runs will go in standard runs/ folder)
    summary_dir = run_dir / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    thresholds = [10.0, 15.0, 20.0, 25.0, 30.0]

    # Results storage: {target: {threshold: {patches, score, identity}}}
    cyno_results = {t: {} for t in targets}
    spec_results = {t: {} for t in targets}

    logger.info("=" * 80)
    logger.info("PARAMETER SWEEP: Cyno Conservation Thresholds")
    logger.info("=" * 80)

    # Sweep 1: Vary cyno conservation (specificity fixed at 15%)
    for threshold in thresholds:
        logger.info("\n" + "=" * 80)
        logger.info(f"Testing cyno threshold: {threshold}% (specificity fixed at 15%)")
        logger.info("=" * 80 + "\n")

        # Temporarily modify config
        original_cyno = config.MAX_CYNO_MISMATCH_PERCENT
        config.MAX_CYNO_MISMATCH_PERCENT = threshold

        try:
            # Run pipeline for each target
            for target in targets:
                logger.info(f"\n--- Running {target} at cyno={threshold}% ---\n")

                # Run pipeline (creates its own dated folder in runs/)
                results = run_pipeline([target])

                # Extract results
                if target in results and results[target]:
                    result = results[target]
                    n_patches = len(result.get('scores', []))
                    best_score = max([s['composite_score'] for s in result.get('scores', [])], default=0.0)

                    # Get overall cyno identity
                    conservation = result.get('conservation')
                    overall_identity = conservation.overall_identity if conservation else 0.0

                    cyno_results[target][threshold] = {
                        'n_patches': n_patches,
                        'best_score': best_score,
                        'overall_identity': overall_identity,
                    }

                    logger.info(f"  {target}: {n_patches} patches, best score {best_score:.3f}, "
                               f"overall cyno identity {overall_identity*100:.1f}%")
                else:
                    cyno_results[target][threshold] = {
                        'n_patches': 0,
                        'best_score': 0.0,
                        'overall_identity': 0.0,
                    }
                    logger.info(f"  {target}: 0 patches")

        finally:
            # Restore original config
            config.MAX_CYNO_MISMATCH_PERCENT = original_cyno

    logger.info("\n" + "=" * 80)
    logger.info("PARAMETER SWEEP: Human Specificity Thresholds")
    logger.info("=" * 80)

    # Sweep 2: Vary human specificity (cyno fixed at 15%)
    for threshold in thresholds:
        logger.info("\n" + "=" * 80)
        logger.info(f"Testing specificity threshold: {threshold}% (cyno fixed at 15%)")
        logger.info("=" * 80 + "\n")

        # Temporarily modify config
        original_spec = config.MAX_NONSPECIFIC_PERCENT
        config.MAX_NONSPECIFIC_PERCENT = threshold

        try:
            # Run pipeline for each target
            for target in targets:
                logger.info(f"\n--- Running {target} at specificity={threshold}% ---\n")

                # Run pipeline (creates its own dated folder in runs/)
                results = run_pipeline([target])

                # Extract results
                if target in results and results[target]:
                    result = results[target]
                    n_patches = len(result.get('scores', []))
                    best_score = max([s['composite_score'] for s in result.get('scores', [])], default=0.0)

                    # Get specificity metrics
                    specificity = result.get('specificity')
                    n_specific = len(specificity.specific_patches) if specificity else 0
                    n_rejected = len(specificity.rejected_patches) if specificity else 0

                    spec_results[target][threshold] = {
                        'n_patches': n_patches,
                        'best_score': best_score,
                        'n_specific': n_specific,
                        'n_rejected': n_rejected,
                    }

                    logger.info(f"  {target}: {n_patches} patches, best score {best_score:.3f}, "
                               f"{n_specific} specific, {n_rejected} rejected")
                else:
                    spec_results[target][threshold] = {
                        'n_patches': 0,
                        'best_score': 0.0,
                        'n_specific': 0,
                        'n_rejected': 0,
                    }
                    logger.info(f"  {target}: 0 patches")

        finally:
            # Restore original config
            config.MAX_NONSPECIFIC_PERCENT = original_spec

    # Generate summary figures
    logger.info("\n" + "=" * 80)
    logger.info("Generating summary figures...")
    logger.info("=" * 80 + "\n")

    _plot_cyno_sweep_summary(cyno_results, thresholds, summary_dir)
    _plot_spec_sweep_summary(spec_results, thresholds, summary_dir)

    # Write summary text report
    _write_summary_report(cyno_results, spec_results, thresholds, summary_dir)

    logger.info(f"\nParameter sweep complete. Results in: {run_dir}")


def _plot_cyno_sweep_summary(results, thresholds, output_dir):
    """Generate summary plots for cyno conservation sweep."""
    targets = list(results.keys())

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle('Cyno Conservation Threshold Sweep\n(Human Specificity fixed at 15%)',
                 fontsize=16, fontweight='bold')

    # Plot 1: Number of patches
    ax1 = axes[0]
    for i, target in enumerate(targets):
        patch_counts = [results[target][t]['n_patches'] for t in thresholds]
        color = CARTO_PURPLE if i == 0 else CARTO_TEAL
        ax1.plot(thresholds, patch_counts, 'o-', linewidth=2, markersize=8,
                color=color, label=target)

    ax1.set_ylabel('Number of Qualifying Patches', fontsize=12, fontweight='bold')
    ax1.set_ylim(bottom=-0.5)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=11, loc='upper left')
    ax1.axvline(x=15.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5,
                label='Current default (15%)')

    # Plot 2: Best composite score
    ax2 = axes[1]
    for i, target in enumerate(targets):
        scores = [results[target][t]['best_score'] for t in thresholds]
        color = CARTO_PURPLE if i == 0 else CARTO_TEAL
        ax2.plot(thresholds, scores, 'o-', linewidth=2, markersize=8,
                color=color, label=target)

    ax2.set_xlabel('Max Cyno Mismatch Percent (%)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Best Composite Score', fontsize=12, fontweight='bold')
    ax2.set_ylim(bottom=-0.05, top=1.05)
    ax2.grid(True, alpha=0.3)
    ax2.axvline(x=15.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)

    plt.tight_layout()
    output_path = output_dir / 'cyno_threshold_sweep.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path.name}")


def _plot_spec_sweep_summary(results, thresholds, output_dir):
    """Generate summary plots for human specificity sweep."""
    targets = list(results.keys())

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle('Human Specificity Threshold Sweep\n(Cyno Conservation fixed at 15%)',
                 fontsize=16, fontweight='bold')

    # Plot 1: Number of patches
    ax1 = axes[0]
    for i, target in enumerate(targets):
        patch_counts = [results[target][t]['n_patches'] for t in thresholds]
        color = CARTO_PURPLE if i == 0 else CARTO_TEAL
        ax1.plot(thresholds, patch_counts, 'o-', linewidth=2, markersize=8,
                color=color, label=target)

    ax1.set_ylabel('Number of Qualifying Patches', fontsize=12, fontweight='bold')
    ax1.set_ylim(bottom=-0.5)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=11, loc='upper left')
    ax1.axvline(x=15.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5,
                label='Current default (15%)')

    # Plot 2: Best composite score
    ax2 = axes[1]
    for i, target in enumerate(targets):
        scores = [results[target][t]['best_score'] for t in thresholds]
        color = CARTO_PURPLE if i == 0 else CARTO_TEAL
        ax2.plot(thresholds, scores, 'o-', linewidth=2, markersize=8,
                color=color, label=target)

    ax2.set_xlabel('Max Non-Specific Percent (%)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Best Composite Score', fontsize=12, fontweight='bold')
    ax2.set_ylim(bottom=-0.05, top=1.05)
    ax2.grid(True, alpha=0.3)
    ax2.axvline(x=15.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)

    plt.tight_layout()
    output_path = output_dir / 'specificity_threshold_sweep.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path.name}")


def _write_summary_report(cyno_results, spec_results, thresholds, output_dir):
    """Write text summary of parameter sweep results."""
    output_path = output_dir / 'parameter_sweep_summary.txt'

    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("PARAMETER SWEEP SUMMARY\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Test Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"Targets: {', '.join(cyno_results.keys())}\n")
        f.write(f"Thresholds tested: {', '.join(f'{t}%' for t in thresholds)}\n\n")

        # Cyno sweep summary
        f.write("=" * 80 + "\n")
        f.write("CYNO CONSERVATION THRESHOLD SWEEP\n")
        f.write("(Human Specificity fixed at 15%)\n")
        f.write("=" * 80 + "\n\n")

        for target in cyno_results:
            f.write(f"{target}:\n")
            f.write(f"  Overall cyno identity: {cyno_results[target][15.0]['overall_identity']*100:.1f}%\n")
            f.write(f"\n  Threshold | Patches | Best Score\n")
            f.write(f"  ----------|---------|------------\n")
            for threshold in thresholds:
                data = cyno_results[target][threshold]
                marker = " <-- current default" if threshold == 15.0 else ""
                f.write(f"  {threshold:>6.0f}%   | {data['n_patches']:>7} | {data['best_score']:>10.3f}{marker}\n")
            f.write("\n")

        # Specificity sweep summary
        f.write("=" * 80 + "\n")
        f.write("HUMAN SPECIFICITY THRESHOLD SWEEP\n")
        f.write("(Cyno Conservation fixed at 15%)\n")
        f.write("=" * 80 + "\n\n")

        for target in spec_results:
            f.write(f"{target}:\n")
            f.write(f"\n  Threshold | Patches | Best Score | Specific | Rejected\n")
            f.write(f"  ----------|---------|------------|----------|----------\n")
            for threshold in thresholds:
                data = spec_results[target][threshold]
                marker = " <-- current default" if threshold == 15.0 else ""
                f.write(f"  {threshold:>6.0f}%   | {data['n_patches']:>7} | {data['best_score']:>10.3f} | "
                       f"{data['n_specific']:>8} | {data['n_rejected']:>8}{marker}\n")
            f.write("\n")

        # Recommendations
        f.write("=" * 80 + "\n")
        f.write("RECOMMENDATIONS\n")
        f.write("=" * 80 + "\n\n")

        f.write("Based on this parameter sweep:\n\n")

        # Analyze cyno results
        f.write("Cyno Conservation:\n")
        for target in cyno_results:
            patch_at_15 = cyno_results[target][15.0]['n_patches']
            if patch_at_15 == 0:
                # Find first threshold that produces patches
                first_working = None
                for t in thresholds:
                    if cyno_results[target][t]['n_patches'] > 0:
                        first_working = t
                        break
                if first_working:
                    f.write(f"  - {target}: 0 patches at 15%, first patches at {first_working}%\n")
                else:
                    f.write(f"  - {target}: 0 patches at all tested thresholds\n")
            else:
                f.write(f"  - {target}: {patch_at_15} patches at 15% (current default working)\n")

        f.write("\nHuman Specificity:\n")
        for target in spec_results:
            patch_at_15 = spec_results[target][15.0]['n_patches']
            if patch_at_15 == 0:
                # Find first threshold that produces patches
                first_working = None
                for t in thresholds:
                    if spec_results[target][t]['n_patches'] > 0:
                        first_working = t
                        break
                if first_working:
                    f.write(f"  - {target}: 0 patches at 15%, first patches at {first_working}%\n")
                else:
                    f.write(f"  - {target}: 0 patches at all tested thresholds\n")
            else:
                f.write(f"  - {target}: {patch_at_15} patches at 15% (current default working)\n")

    print(f"  Saved: {output_path.name}")


def main():
    """Run parameter sweep on CEACAM5 and ERBB2."""
    # Setup
    timestamp = datetime.now().strftime("%y%m%d_%H%M")
    run_name = f"{timestamp}_parameter_sweep_ceacam5_erbb2"
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
    targets = ['CEACAM5', 'ERBB2']
    run_threshold_sweep(targets, run_dir)


if __name__ == '__main__':
    main()
