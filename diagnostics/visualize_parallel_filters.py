#!/usr/bin/env python3
"""
Diagnostic visualization showing how cyno conservation and human specificity
filters operate independently vs sequentially.

This script:
1. Loads pre-filtered residues (spatial, surface, distance)
2. Runs cyno conservation filter independently
3. Runs specificity filter independently (on ALL qualifying residues)
4. Shows intersection and comparison with sequential results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.patches as mpatches

# Cartography colors
CARTO_MINT = '#E0F3DB'
CARTO_RED = '#FF6B6B'
CARTO_TEAL = '#18B5CB'
CARTO_GRAY = '#D3D3D3'
CARTO_PURPLE = '#28154C'

def load_ceacam5_residues():
    """Load CEACAM5 residue table from previous run."""
    csv_path = Path("epitope_pipeline/runs/260226_1737_ceacam5/annotated_sequences/ceacam5_residue_table.csv")
    if not csv_path.exists():
        raise FileNotFoundError(f"CEACAM5 residue table not found: {csv_path}")

    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} CEACAM5 residues")
    return df

def apply_spatial_filters(df):
    """Apply spatial/surface/distance filters (Steps 1-5)."""
    # Extracellular only (case-insensitive)
    df_filtered = df[df['topology'].str.lower() == 'extracellular'].copy()
    print(f"After extracellular filter: {len(df_filtered)} residues")

    # Surface-exposed (SASA > 0)
    df_filtered = df_filtered[df_filtered['sasa_A2'] > 0].copy()
    print(f"After surface filter: {len(df_filtered)} residues")

    # Distance >= 80A from membrane
    df_filtered = df_filtered[df_filtered['distance_from_membrane_A'] >= 80].copy()
    print(f"After distance filter: {len(df_filtered)} residues")

    return df_filtered

def analyze_cyno_conservation(df):
    """Analyze cyno conservation status."""
    # cyno_conserved column: 'yes', 'no', or NaN/blank
    cyno_yes = df['cyno_conserved'] == 'yes'
    cyno_no = df['cyno_conserved'] == 'no'
    cyno_unknown = ~(cyno_yes | cyno_no)

    print(f"\nCyno Conservation:")
    print(f"  Conserved (yes): {cyno_yes.sum()}")
    print(f"  Mismatch (no): {cyno_no.sum()}")
    print(f"  Unknown/NaN: {cyno_unknown.sum()}")

    return {
        'conserved': df[cyno_yes].copy(),
        'mismatch': df[cyno_no].copy(),
        'unknown': df[cyno_unknown].copy()
    }

def analyze_specificity(df):
    """Analyze human specificity status."""
    # human_specific column: 'yes', 'no', or NaN/blank
    spec_yes = df['human_specific'] == 'yes'
    spec_no = df['human_specific'] == 'no'
    spec_unknown = ~(spec_yes | spec_no)

    print(f"\nHuman Specificity:")
    print(f"  Human-specific (yes): {spec_yes.sum()}")
    print(f"  Non-specific (no): {spec_no.sum()}")
    print(f"  Not assessed (NaN): {spec_unknown.sum()}")

    return {
        'specific': df[spec_yes].copy(),
        'nonspecific': df[spec_no].copy(),
        'not_assessed': df[spec_unknown].copy()
    }

def plot_filter_comparison(df_spatial, cyno_results, spec_results):
    """Create visualization comparing sequential vs parallel filtering."""

    fig, axes = plt.subplots(4, 1, figsize=(14, 12))
    fig.suptitle('CEACAM5: Sequential vs Parallel Filtering',
                 fontsize=16, fontweight='bold', y=0.995)

    residue_nums = df_spatial['residue_num'].values

    # Track 1: Spatial filtering results (baseline)
    ax1 = axes[0]
    ax1.scatter(residue_nums, [1]*len(residue_nums),
                c=CARTO_TEAL, s=20, alpha=0.6, marker='s')
    ax1.set_ylim(0.5, 1.5)
    ax1.set_yticks([])
    ax1.set_ylabel('Spatial\nFiltering', fontsize=11, fontweight='bold')
    ax1.set_title(f'Step 1-5: Spatial/Surface/Distance ≥80A → {len(df_spatial)} qualifying residues',
                  fontsize=11, loc='left')
    ax1.grid(axis='x', alpha=0.3)
    ax1.set_xlim(0, 702)

    # Track 2: Cyno conservation (independent)
    ax2 = axes[1]
    cyno_conserved = cyno_results['conserved']['residue_num'].values
    cyno_mismatch = cyno_results['mismatch']['residue_num'].values

    ax2.scatter(cyno_conserved, [1]*len(cyno_conserved),
                c=CARTO_MINT, s=20, alpha=0.8, marker='s', label='Cyno conserved')
    ax2.scatter(cyno_mismatch, [1]*len(cyno_mismatch),
                c=CARTO_RED, s=20, alpha=0.8, marker='s', label='Cyno mismatch')

    ax2.set_ylim(0.5, 1.5)
    ax2.set_yticks([])
    ax2.set_ylabel('Cyno\nConservation', fontsize=11, fontweight='bold')
    ax2.set_title(f'Step 6A: Cyno Conservation → {len(cyno_conserved)} conserved, {len(cyno_mismatch)} mismatch ({len(cyno_mismatch)/len(df_spatial)*100:.1f}%)',
                  fontsize=11, loc='left')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(axis='x', alpha=0.3)
    ax2.set_xlim(0, 702)

    # Track 3: Human specificity (independent)
    ax3 = axes[2]
    spec_specific = spec_results['specific']['residue_num'].values
    spec_nonspecific = spec_results['nonspecific']['residue_num'].values
    spec_not_assessed = spec_results['not_assessed']['residue_num'].values

    ax3.scatter(spec_specific, [1]*len(spec_specific),
                c=CARTO_MINT, s=20, alpha=0.8, marker='s', label='Human-specific')
    if len(spec_nonspecific) > 0:
        ax3.scatter(spec_nonspecific, [1]*len(spec_nonspecific),
                    c=CARTO_RED, s=20, alpha=0.8, marker='s', label='Non-specific')
    ax3.scatter(spec_not_assessed, [1]*len(spec_not_assessed),
                c=CARTO_GRAY, s=20, alpha=0.5, marker='s', label='Not assessed')

    ax3.set_ylim(0.5, 1.5)
    ax3.set_yticks([])
    ax3.set_ylabel('Human\nSpecificity', fontsize=11, fontweight='bold')
    ax3.set_title(f'Step 6B: Human Specificity (BLAST) → {len(spec_specific)} specific, {len(spec_nonspecific)} non-specific, {len(spec_not_assessed)} not assessed',
                  fontsize=11, loc='left')
    ax3.legend(loc='upper right', fontsize=9)
    ax3.grid(axis='x', alpha=0.3)
    ax3.set_xlim(0, 702)

    # Track 4: Intersection (both filters pass)
    ax4 = axes[3]

    # Find residues passing BOTH filters
    cyno_pass_set = set(cyno_conserved)
    spec_pass_set = set(spec_specific)
    intersection = sorted(cyno_pass_set & spec_pass_set)

    # Find residues failing each filter
    cyno_fail_set = set(cyno_mismatch)
    spec_fail_set = set(spec_nonspecific)

    # Residues failing cyno only
    cyno_only_fail = sorted(cyno_fail_set - spec_fail_set)
    # Residues failing specificity only
    spec_only_fail = sorted(spec_fail_set - cyno_fail_set)
    # Residues failing both
    both_fail = sorted(cyno_fail_set & spec_fail_set)

    ax4.scatter(intersection, [1]*len(intersection),
                c=CARTO_MINT, s=20, alpha=0.8, marker='s', label='Pass both')

    if len(cyno_only_fail) > 0:
        ax4.scatter(cyno_only_fail, [0.9]*len(cyno_only_fail),
                    c='#FFA07A', s=15, alpha=0.6, marker='v', label='Fail cyno only')
    if len(spec_only_fail) > 0:
        ax4.scatter(spec_only_fail, [0.85]*len(spec_only_fail),
                    c='#FF8C42', s=15, alpha=0.6, marker='^', label='Fail specificity only')
    if len(both_fail) > 0:
        ax4.scatter(both_fail, [0.8]*len(both_fail),
                    c=CARTO_RED, s=20, alpha=0.8, marker='X', label='Fail both')

    ax4.set_ylim(0.75, 1.5)
    ax4.set_yticks([])
    ax4.set_ylabel('Parallel\nIntersection', fontsize=11, fontweight='bold')
    ax4.set_xlabel('Residue Number', fontsize=11, fontweight='bold')
    ax4.set_title(f'Step 7: Parallel Filter Intersection → {len(intersection)} residues pass BOTH filters',
                  fontsize=11, loc='left')
    ax4.legend(loc='upper right', fontsize=9, ncol=2)
    ax4.grid(axis='x', alpha=0.3)
    ax4.set_xlim(0, 702)

    plt.tight_layout()

    # Save figure
    output_dir = Path("epitope_pipeline/runs/260310_1148_specificity_sweep/figures")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "ceacam5_parallel_filter_analysis.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved filter comparison: {output_path}")

    return {
        'pass_both': len(intersection),
        'fail_cyno_only': len(cyno_only_fail),
        'fail_spec_only': len(spec_only_fail),
        'fail_both': len(both_fail)
    }

def analyze_not_assessed_overlap(df_spatial, cyno_results, spec_results):
    """Analyze the 'not assessed' specificity residues - are they cyno mismatches or orphans?"""

    not_assessed_nums = set(spec_results['not_assessed']['residue_num'].values)
    cyno_mismatch_nums = set(cyno_results['mismatch']['residue_num'].values)
    cyno_conserved_nums = set(cyno_results['conserved']['residue_num'].values)

    # Split not_assessed by cyno status
    not_assessed_but_cyno_mismatch = sorted(not_assessed_nums & cyno_mismatch_nums)
    not_assessed_but_cyno_conserved = sorted(not_assessed_nums & cyno_conserved_nums)

    print(f"\n{'='*60}")
    print("CRITICAL FINDING: 'Not Assessed' Specificity Breakdown")
    print(f"{'='*60}")
    print(f"Total 'not assessed' for specificity: {len(not_assessed_nums)}")
    print(f"  └─ Cyno MISMATCH (explains absence): {len(not_assessed_but_cyno_mismatch)}")
    print(f"  └─ Cyno CONSERVED (orphaned!): {len(not_assessed_but_cyno_conserved)}")
    print(f"\n{len(not_assessed_but_cyno_conserved)} cyno-conserved residues never got BLAST-screened!")
    print("These were likely orphaned in sub-threshold clusters after cyno reclustering.")

    return {
        'not_assessed_cyno_mismatch': len(not_assessed_but_cyno_mismatch),
        'not_assessed_cyno_conserved': len(not_assessed_but_cyno_conserved)
    }

def main():
    print("="*60)
    print("CEACAM5 Parallel vs Sequential Filter Analysis")
    print("="*60)

    # Load residue data
    df = load_ceacam5_residues()

    # Apply spatial filters (Steps 1-5)
    df_spatial = apply_spatial_filters(df)

    # Analyze cyno conservation (Step 6A - independent)
    cyno_results = analyze_cyno_conservation(df_spatial)

    # Analyze human specificity (Step 6B - independent)
    spec_results = analyze_specificity(df_spatial)

    # Analyze the overlap of "not assessed" with cyno status
    orphan_stats = analyze_not_assessed_overlap(df_spatial, cyno_results, spec_results)

    # Create visualization
    intersection_stats = plot_filter_comparison(df_spatial, cyno_results, spec_results)

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY: Parallel Filtering Results")
    print(f"{'='*60}")
    print(f"Starting residues (spatial filter): {len(df_spatial)}")
    print(f"  Cyno conserved: {len(cyno_results['conserved'])} ({len(cyno_results['conserved'])/len(df_spatial)*100:.1f}%)")
    print(f"  Cyno mismatch: {len(cyno_results['mismatch'])} ({len(cyno_results['mismatch'])/len(df_spatial)*100:.1f}%)")
    print(f"  Human-specific: {len(spec_results['specific'])} ({len(spec_results['specific'])/len(df_spatial)*100:.1f}%)")
    print(f"  Human non-specific: {len(spec_results['nonspecific'])} ({len(spec_results['nonspecific'])/len(df_spatial)*100:.1f}%)")
    print(f"\nIntersection (pass BOTH filters): {intersection_stats['pass_both']} residues")
    print(f"  Fail cyno only: {intersection_stats['fail_cyno_only']}")
    print(f"  Fail specificity only: {intersection_stats['fail_spec_only']}")
    print(f"  Fail both: {intersection_stats['fail_both']}")

    print(f"\n{'='*60}")
    print("KEY INSIGHT: Sequential vs Parallel")
    print(f"{'='*60}")
    print("Current SEQUENTIAL pipeline (Cyno → Specificity):")
    print(f"  • {len(cyno_results['conserved'])} residues pass cyno → sent to BLAST")
    print(f"  • {len(cyno_results['mismatch'])} residues fail cyno → NEVER BLAST-screened")
    print(f"  • Result: incomplete data, can't distinguish filters")
    print("")
    print("Proposed PARALLEL pipeline (Cyno + Specificity independently):")
    print(f"  • ALL {len(df_spatial)} residues get cyno assessment")
    print(f"  • ALL {len(df_spatial)} residues get BLAST assessment")
    print(f"  • Intersection: {intersection_stats['pass_both']} residues pass BOTH")
    print(f"  • Result: complete data, can definitively attribute failures")

if __name__ == '__main__':
    main()
