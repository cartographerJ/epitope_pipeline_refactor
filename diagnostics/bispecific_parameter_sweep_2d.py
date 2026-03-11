"""
2D Parameter sweep for bispecific pairs - test ALL combinations of cyno and specificity thresholds.

Tests CEACAM5:ERBB2 at all combinations:
- Cyno conservation: 0%, 20%, 40%, 60%, 80%, 100%
- Human specificity: 0%, 20%, 40%, 60%, 80%, 100%
- Total: 6×6 = 36 threshold combinations

All results contained in single timestamped run folder.
"""

import logging
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Add parent directory to path
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


def run_bispecific_2d_sweep(pair_name, run_dir):
    """
    Run full 2D parameter sweep testing all combinations of cyno and specificity thresholds.

    Args:
        pair_name: Pair in format "GENE_A:GENE_B"
        run_dir: Output directory for all results
    """
    logger = logging.getLogger(__name__)

    # Parse pair
    gene_a, gene_b = pair_name.split(':')

    # Create subdirectories
    sweep_dir = run_dir / "threshold_combinations"
    summary_dir = run_dir / "summary"
    sweep_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    thresholds = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0]

    # Results storage: {(cyno, spec): {pair_score, both_valid, patch_counts, ...}}
    results = {}

    logger.info("=" * 80)
    logger.info(f"2D PARAMETER SWEEP: {pair_name}")
    logger.info("=" * 80)
    logger.info(f"Testing {len(thresholds)}×{len(thresholds)} = {len(thresholds)**2} combinations")
    logger.info(f"Output directory: {run_dir}")
    logger.info("=" * 80)

    total = len(thresholds) ** 2
    count = 0

    for cyno_thresh in thresholds:
        for spec_thresh in thresholds:
            count += 1
            logger.info("\n" + "-" * 80)
            logger.info(f"[{count}/{total}] Testing cyno={cyno_thresh}%, spec={spec_thresh}%")
            logger.info("-" * 80 + "\n")

            try:
                # Create run name
                custom_run_name = f"{run_dir.name}/threshold_combinations/cyno_{int(cyno_thresh)}pct_spec_{int(spec_thresh)}pct"

                # Pass thresholds as parameters to run_bispecific
                output = run_bispecific(
                    [(gene_a, gene_b)],
                    run_name=custom_run_name,
                    cyno_mismatch_percent=cyno_thresh,
                    nonspecific_percent=spec_thresh
                )
                pair_results = output.get("pair_results", [])

                if pair_results and len(pair_results) > 0:
                    result = pair_results[0]

                    # Extract metrics
                    ab_score = result.orientation_ab.orientation_score if result.orientation_ab else 0.0
                    ba_score = result.orientation_ba.orientation_score if result.orientation_ba else 0.0

                    gene_a_distal = len(result.orientation_ab.distal_zone.specificity_result.specific_patches) if result.orientation_ab and result.orientation_ab.distal_zone.specificity_result else 0
                    gene_b_proximal = len(result.orientation_ab.proximal_zone.specificity_result.specific_patches) if result.orientation_ab and result.orientation_ab.proximal_zone.specificity_result else 0
                    gene_a_proximal = len(result.orientation_ba.proximal_zone.specificity_result.specific_patches) if result.orientation_ba and result.orientation_ba.proximal_zone.specificity_result else 0
                    gene_b_distal = len(result.orientation_ba.distal_zone.specificity_result.specific_patches) if result.orientation_ba and result.orientation_ba.distal_zone.specificity_result else 0

                    results[(cyno_thresh, spec_thresh)] = {
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
                else:
                    results[(cyno_thresh, spec_thresh)] = {
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

            except Exception as e:
                logger.error(f"  Error running bispecific analysis: {e}")
                results[(cyno_thresh, spec_thresh)] = {
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

    # Generate summary outputs
    logger.info("\n" + "=" * 80)
    logger.info("Generating summary outputs...")
    logger.info("=" * 80 + "\n")

    _plot_2d_heatmaps(results, thresholds, pair_name, summary_dir)
    _write_2d_summary_report(results, thresholds, pair_name, summary_dir)

    logger.info(f"\n2D parameter sweep complete. Results in: {run_dir}")


def _plot_2d_heatmaps(results, thresholds, pair_name, output_dir):
    """Generate 2D heatmaps showing pair scores and patch counts across parameter space."""

    # Build matrices
    n = len(thresholds)
    pair_score_matrix = np.zeros((n, n))
    ceacam5_distal_matrix = np.zeros((n, n))
    ceacam5_proximal_matrix = np.zeros((n, n))
    erbb2_distal_matrix = np.zeros((n, n))
    erbb2_proximal_matrix = np.zeros((n, n))

    for i, cyno in enumerate(thresholds):
        for j, spec in enumerate(thresholds):
            data = results.get((cyno, spec), {})
            pair_score_matrix[i, j] = data.get('pair_score', 0.0)
            ceacam5_distal_matrix[i, j] = data.get('gene_a_distal_patches', 0)
            ceacam5_proximal_matrix[i, j] = data.get('gene_a_proximal_patches', 0)
            erbb2_distal_matrix[i, j] = data.get('gene_b_distal_patches', 0)
            erbb2_proximal_matrix[i, j] = data.get('gene_b_proximal_patches', 0)

    # Plot 1: Pair score heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(pair_score_matrix, cmap='RdYlGn', vmin=0, vmax=1, origin='lower')
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels([f'{int(t)}%' for t in thresholds])
    ax.set_yticklabels([f'{int(t)}%' for t in thresholds])
    ax.set_xlabel('Human Specificity Threshold', fontsize=12, fontweight='bold')
    ax.set_ylabel('Cyno Conservation Threshold', fontsize=12, fontweight='bold')
    ax.set_title(f'{pair_name} Bispecific Pair Score\n2D Parameter Sweep', fontsize=14, fontweight='bold')

    # Add text annotations
    for i in range(n):
        for j in range(n):
            text = ax.text(j, i, f'{pair_score_matrix[i, j]:.2f}',
                          ha="center", va="center", color="black", fontsize=9)

    # Mark default (15%, 15%) if in range
    if 15.0 in thresholds:
        default_idx = thresholds.index(15.0)
        ax.add_patch(plt.Rectangle((default_idx-0.5, default_idx-0.5), 1, 1,
                                   fill=False, edgecolor='blue', linewidth=3))

    plt.colorbar(im, ax=ax, label='Pair Score')
    plt.tight_layout()
    output_path = output_dir / 'pair_score_heatmap.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path.name}")

    # Plot 2: Patch count heatmaps (2x2 grid)
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(f'{pair_name} Patch Counts by Zone\n2D Parameter Sweep', fontsize=14, fontweight='bold')

    matrices = [
        (ceacam5_distal_matrix, 'CEACAM5 Distal', axes[0, 0]),
        (ceacam5_proximal_matrix, 'CEACAM5 Proximal', axes[0, 1]),
        (erbb2_distal_matrix, 'ERBB2 Distal', axes[1, 0]),
        (erbb2_proximal_matrix, 'ERBB2 Proximal', axes[1, 1]),
    ]

    for matrix, title, ax in matrices:
        im = ax.imshow(matrix, cmap='Blues', vmin=0, vmax=max(3, matrix.max()), origin='lower')
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels([f'{int(t)}%' for t in thresholds])
        ax.set_yticklabels([f'{int(t)}%' for t in thresholds])
        ax.set_xlabel('Specificity Threshold', fontsize=10)
        ax.set_ylabel('Cyno Threshold', fontsize=10)
        ax.set_title(title, fontsize=11, fontweight='bold')

        # Add text annotations
        for i in range(n):
            for j in range(n):
                ax.text(j, i, f'{int(matrix[i, j])}',
                       ha="center", va="center",
                       color="white" if matrix[i, j] > 1.5 else "black",
                       fontsize=8)

        # Mark default if in range
        if 15.0 in thresholds:
            default_idx = thresholds.index(15.0)
            ax.add_patch(plt.Rectangle((default_idx-0.5, default_idx-0.5), 1, 1,
                                       fill=False, edgecolor='red', linewidth=2))
        plt.colorbar(im, ax=ax, label='# Patches')

    plt.tight_layout()
    output_path = output_dir / 'patch_counts_heatmaps.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path.name}")


def _write_2d_summary_report(results, thresholds, pair_name, output_dir):
    """Write text summary of 2D parameter sweep."""
    output_path = output_dir / 'parameter_sweep_2d_summary.txt'

    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(f"2D PARAMETER SWEEP: {pair_name}\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Test Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"Thresholds tested: {', '.join(f'{int(t)}%' for t in thresholds)}\n")
        f.write(f"Total combinations: {len(thresholds)}×{len(thresholds)} = {len(results)}\n\n")

        # Matrix format
        f.write("=" * 80 + "\n")
        f.write("PAIR SCORE MATRIX\n")
        f.write("=" * 80 + "\n\n")

        f.write("Cyno↓ / Spec→ |")
        for spec in thresholds:
            f.write(f"  {int(spec):>3}% |")
        f.write("\n")
        f.write("-" * 80 + "\n")

        for cyno in thresholds:
            f.write(f"    {int(cyno):>3}%      |")
            for spec in thresholds:
                score = results.get((cyno, spec), {}).get('pair_score', 0.0)
                marker = " *" if (cyno == 15.0 and spec == 15.0 and 15.0 in thresholds) else ""
                f.write(f" {score:>5.3f}{marker} |")
            f.write("\n")

        if 15.0 in thresholds:
            f.write("\n* = default (15%, 15%)\n\n")
        else:
            f.write("\n")

        # Patch count matrices
        gene_a, gene_b = pair_name.split(':')

        for zone, key in [
            (f"{gene_a} Distal", 'gene_a_distal_patches'),
            (f"{gene_a} Proximal", 'gene_a_proximal_patches'),
            (f"{gene_b} Distal", 'gene_b_distal_patches'),
            (f"{gene_b} Proximal", 'gene_b_proximal_patches'),
        ]:
            f.write("=" * 80 + "\n")
            f.write(f"{zone.upper()} PATCH COUNTS\n")
            f.write("=" * 80 + "\n\n")

            f.write("Cyno↓ / Spec→ |")
            for spec in thresholds:
                f.write(f" {int(spec):>3}% |")
            f.write("\n")
            f.write("-" * 80 + "\n")

            for cyno in thresholds:
                f.write(f"    {int(cyno):>3}%      |")
                for spec in thresholds:
                    count = results.get((cyno, spec), {}).get(key, 0)
                    marker = " *" if (cyno == 15.0 and spec == 15.0 and 15.0 in thresholds) else ""
                    f.write(f" {count:>4}{marker} |")
                f.write("\n")
            f.write("\n")

    print(f"  Saved: {output_path.name}")


def main():
    """Run 2D parameter sweep on CEACAM5:ERBB2."""
    timestamp = datetime.now().strftime("%y%m%d_%H%M")
    run_name = f"{timestamp}_param_sweep_2d_ceacam5_erbb2"
    run_dir = Path(__file__).parent.parent / "runs" / run_name
    run_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    log_file = run_dir / "parameter_sweep_2d.log"
    setup_logging(verbose=True)

    # Add file handler
    logger = logging.getLogger()
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    ))
    logger.addHandler(file_handler)

    # Run 2D sweep
    pair_name = 'CEACAM5:ERBB2'
    run_bispecific_2d_sweep(pair_name, run_dir)


if __name__ == '__main__':
    main()
