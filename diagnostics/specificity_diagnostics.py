#!/usr/bin/env python3
"""
Diagnostic plotting for specificity sliding window analysis.

Shows how the per-residue identity scores and sliding window scans
determine which patches pass/fail specificity screening.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Cartography colors
CARTO_MINT = '#E0F3DB'
CARTO_GREEN = '#6BC291'
CARTO_TEAL = '#18B5CB'
CARTO_BLUE = '#2E95D2'
CARTO_PURPLE = '#28154C'
CARTO_RED = '#FF6B6B'
CARTO_GRAY = '#D3D3D3'
CARTO_AMBER = '#FFB74D'


def plot_specificity_diagnostics(
    gene_name,
    residue_identity,
    patches_before_trim,
    patches_after_trim,
    window_counts,
    threshold,
    max_nonspecific_per_window,
    output_dir,
):
    """
    Generate diagnostic plots showing specificity sliding window analysis.

    Args:
        gene_name: Target gene name.
        residue_identity: Dict {resnum: float} — BLAST identity scores (0.0-1.0)
        patches_before_trim: List of SurfacePatch objects before trimming.
        patches_after_trim: List of SurfacePatch objects after trimming.
        window_counts: Dict {resnum: int} — count of non-specific residues in window.
        threshold: Per-residue HSP identity threshold (legacy; the main
            pipeline now uses per-paralog patch match fractions).
        max_nonspecific_per_window: MAX_NONSPECIFIC_PER_600A2 (e.g., 2).
        output_dir: Path to save diagnostic plots.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Convert to arrays for plotting
    all_residues = sorted(residue_identity.keys())
    identity_scores = [residue_identity.get(r, 0.0) for r in all_residues]
    window_count_values = [window_counts.get(r, 0) for r in all_residues]

    # Create 3-panel diagnostic plot
    fig, axes = plt.subplots(3, 1, figsize=(16, 10))
    fig.suptitle(f'{gene_name}: Specificity Sliding Window Diagnostics',
                 fontsize=16, fontweight='bold', y=0.995)

    # Track 1: Per-residue identity scores
    ax1 = axes[0]

    # Color by threshold
    colors = [CARTO_RED if s >= threshold else CARTO_MINT for s in identity_scores]
    ax1.bar(all_residues, identity_scores, color=colors, width=1.0, alpha=0.7)
    ax1.axhline(threshold, color='black', linestyle='--', linewidth=1.5,
                label=f'Threshold ({threshold:.0%})')

    ax1.set_ylabel('BLAST Identity', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, 1.0)
    ax1.set_title(f'Per-Residue Identity to Human Proteome (≥{threshold:.0%} = non-specific)',
                  fontsize=11, loc='left')
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(axis='y', alpha=0.3)
    ax1.set_xlim(0, max(all_residues) + 10)

    # Track 2: Sliding window counts
    ax2 = axes[1]

    # Color by max threshold
    colors = [CARTO_RED if c > max_nonspecific_per_window else CARTO_MINT
              for c in window_count_values]
    ax2.bar(all_residues, window_count_values, color=colors, width=1.0, alpha=0.7)
    ax2.axhline(max_nonspecific_per_window, color='black', linestyle='--',
                linewidth=1.5, label=f'Max allowed ({max_nonspecific_per_window})')

    ax2.set_ylabel('Non-Specific\nResidues in\n600A² Window',
                   fontsize=12, fontweight='bold')
    ax2.set_title(f'Sliding Window Scan (>{max_nonspecific_per_window} non-specific → trim residue)',
                  fontsize=11, loc='left')
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(axis='y', alpha=0.3)
    ax2.set_xlim(0, max(all_residues) + 10)

    # Track 3: Patch boundaries (before/after trimming)
    ax3 = axes[2]

    # Draw patches before trimming (gray background)
    for i, patch in enumerate(patches_before_trim):
        resnums = patch.residue_numbers
        if not resnums:
            continue
        min_res = min(resnums)
        max_res = max(resnums)
        width = max_res - min_res + 1
        ax3.add_patch(mpatches.Rectangle(
            (min_res - 0.5, 0.4), width, 0.2,
            facecolor=CARTO_GRAY, edgecolor='black', linewidth=1.5, alpha=0.4,
            label='Before trim' if i == 0 else None
        ))

    # Draw patches after trimming (green)
    for i, patch in enumerate(patches_after_trim):
        resnums = patch.residue_numbers
        if not resnums:
            continue
        min_res = min(resnums)
        max_res = max(resnums)
        width = max_res - min_res + 1
        ax3.add_patch(mpatches.Rectangle(
            (min_res - 0.5, 0.6), width, 0.2,
            facecolor=CARTO_GREEN, edgecolor='black', linewidth=1.5, alpha=0.8,
            label='After trim' if i == 0 else None
        ))

    ax3.set_ylim(0, 1.2)
    ax3.set_yticks([])
    ax3.set_ylabel('Patches', fontsize=12, fontweight='bold')
    ax3.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
    ax3.set_title(f'Patch Boundaries: {len(patches_before_trim)} before → {len(patches_after_trim)} after trimming',
                  fontsize=11, loc='left')
    ax3.legend(loc='upper right', fontsize=10)
    ax3.grid(axis='x', alpha=0.3)
    ax3.set_xlim(0, max(all_residues) + 10)

    plt.tight_layout()

    # Save figure
    output_path = output_dir / f"{gene_name.lower()}_specificity_diagnostics.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved specificity diagnostics: {output_path.name}")

    # Also create a summary stats file
    stats_path = output_dir / f"{gene_name.lower()}_specificity_stats.txt"
    with open(stats_path, 'w') as f:
        f.write(f"Specificity Diagnostics: {gene_name}\n")
        f.write(f"{'='*60}\n\n")

        f.write(f"Threshold Settings:\n")
        f.write(f"  Identity threshold: {threshold:.1%}\n")
        f.write(f"  Max non-specific per 600A²: {max_nonspecific_per_window}\n\n")

        f.write(f"Residue-Level Stats:\n")
        n_assessed = len([s for s in identity_scores if s > 0])
        n_specific = len([s for s in identity_scores if 0 < s < threshold])
        n_nonspecific = len([s for s in identity_scores if s >= threshold])
        f.write(f"  Total assessed: {n_assessed}\n")
        f.write(f"  Specific (<{threshold:.0%}): {n_specific} ({n_specific/n_assessed*100 if n_assessed else 0:.1f}%)\n")
        f.write(f"  Non-specific (≥{threshold:.0%}): {n_nonspecific} ({n_nonspecific/n_assessed*100 if n_assessed else 0:.1f}%)\n\n")

        f.write(f"Sliding Window Stats:\n")
        n_failing = len([c for c in window_count_values if c > max_nonspecific_per_window])
        f.write(f"  Residues exceeding threshold: {n_failing}\n\n")

        f.write(f"Patch-Level Stats:\n")
        f.write(f"  Patches before trimming: {len(patches_before_trim)}\n")
        if patches_before_trim:
            total_before = sum(len(p.residue_numbers) for p in patches_before_trim)
            f.write(f"  Total residues before: {total_before}\n")
        f.write(f"  Patches after trimming: {len(patches_after_trim)}\n")
        if patches_after_trim:
            total_after = sum(len(p.residue_numbers) for p in patches_after_trim)
            f.write(f"  Total residues after: {total_after}\n")
            f.write(f"  Residues trimmed: {total_before - total_after if patches_before_trim else 0}\n")

    print(f"  Saved specificity stats: {stats_path.name}")


def plot_specificity_comparison(
    gene_name,
    residue_identity,
    patches_old_method,
    patches_new_method,
    threshold,
    output_dir,
):
    """
    Compare old (HSP range) vs new (per-residue) specificity filtering.

    Args:
        gene_name: Target gene name.
        residue_identity: Dict {resnum: float} — BLAST identity scores.
        patches_old_method: Patches from old HSP-range method.
        patches_new_method: Patches from new per-residue method.
        threshold: Per-residue HSP identity threshold (legacy).
        output_dir: Path to save comparison plot.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(16, 8))
    fig.suptitle(f'{gene_name}: Old vs New Specificity Filtering',
                 fontsize=16, fontweight='bold')

    all_residues = sorted(residue_identity.keys())
    identity_scores = [residue_identity.get(r, 0.0) for r in all_residues]

    # Track 1: Old method (HSP range - BROKEN)
    ax1 = axes[0]
    colors = [CARTO_RED if s >= threshold else CARTO_MINT for s in identity_scores]
    ax1.bar(all_residues, identity_scores, color=colors, width=1.0, alpha=0.5)

    # Draw old patches
    for i, patch in enumerate(patches_old_method):
        resnums = patch.residue_numbers
        if not resnums:
            continue
        min_res = min(resnums)
        max_res = max(resnums)
        ax1.axvspan(min_res, max_res, alpha=0.2, color=CARTO_GRAY,
                   label='Old method patch' if i == 0 else None)

    ax1.axhline(threshold, color='black', linestyle='--', linewidth=1.5)
    ax1.set_ylabel('BLAST Identity', fontsize=12, fontweight='bold')
    ax1.set_title('OLD: HSP Range Overlap (patch-based binary)',
                  fontsize=11, loc='left')
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(axis='y', alpha=0.3)
    ax1.set_xlim(0, max(all_residues) + 10)

    # Track 2: New method (per-residue + sliding window)
    ax2 = axes[1]
    colors = [CARTO_RED if s >= threshold else CARTO_MINT for s in identity_scores]
    ax2.bar(all_residues, identity_scores, color=colors, width=1.0, alpha=0.5)

    # Draw new patches
    for i, patch in enumerate(patches_new_method):
        resnums = patch.residue_numbers
        if not resnums:
            continue
        min_res = min(resnums)
        max_res = max(resnums)
        ax2.axvspan(min_res, max_res, alpha=0.3, color=CARTO_GREEN,
                   label='New method patch' if i == 0 else None)

    ax2.axhline(threshold, color='black', linestyle='--', linewidth=1.5)
    ax2.set_ylabel('BLAST Identity', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
    ax2.set_title('NEW: Per-Residue Identity + Sliding Window Trim',
                  fontsize=11, loc='left')
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(axis='y', alpha=0.3)
    ax2.set_xlim(0, max(all_residues) + 10)

    plt.tight_layout()

    output_path = output_dir / f"{gene_name.lower()}_specificity_comparison.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved specificity comparison: {output_path.name}")
