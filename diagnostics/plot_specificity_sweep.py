#!/usr/bin/env python3
"""
Bespoke plot: Available patches vs. specificity parameter conditions.
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys

# Cartography palette
PALETTE = {
    "light_gray":  "#D3D3D3",
    "mint":        "#E0F3DB",
    "green":       "#6BC291",
    "teal":        "#18B5CB",
    "blue":        "#2E95D2",
    "dark_purple": "#28154C",
}

# Find most recent results file
results_dir = Path("epitope_pipeline/specificity_sweep_results")
if not results_dir.exists():
    print(f"ERROR: Results directory not found: {results_dir}")
    sys.exit(1)

csv_files = list(results_dir.glob("param_sweep_*.csv"))
if not csv_files:
    print(f"ERROR: No results files found in {results_dir}")
    sys.exit(1)

latest_csv = max(csv_files, key=lambda p: p.stat().st_mtime)
print(f"Loading results from: {latest_csv}")

# Load data
df = pd.read_csv(latest_csv)

# Convert scores > 0 to "has patches" (1) vs "no patches" (0)
df['has_patches'] = (df['pair_score'] > 0).astype(int)

# Create figure with subplots
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle('Epitope Patch Availability vs. Specificity Parameters',
             fontsize=16, fontweight='bold', y=0.98)

pairs = df['pair'].unique()
max_nonspecific_values = sorted(df['max_nonspecific'].unique())
specificity_thresholds = sorted(df['specificity_threshold'].unique())

for idx, pair in enumerate(pairs):
    ax = axes[idx]
    pair_df = df[df['pair'] == pair]

    # Create heatmap data
    heatmap_data = np.zeros((len(max_nonspecific_values), len(specificity_thresholds)))

    for i, max_nonspec in enumerate(max_nonspecific_values):
        for j, spec_thresh in enumerate(specificity_thresholds):
            row = pair_df[
                (pair_df['max_nonspecific'] == max_nonspec) &
                (pair_df['specificity_threshold'] == spec_thresh)
            ]
            if not row.empty:
                heatmap_data[i, j] = row.iloc[0]['pair_score']

    # Create heatmap
    im = ax.imshow(heatmap_data, cmap='RdYlGn', aspect='auto',
                   vmin=0, vmax=1.2, interpolation='nearest')

    # Set ticks and labels
    ax.set_xticks(np.arange(len(specificity_thresholds)))
    ax.set_yticks(np.arange(len(max_nonspecific_values)))
    ax.set_xticklabels([f'{t:.0%}' for t in specificity_thresholds], fontsize=10)
    ax.set_yticklabels([f'{m}' for m in max_nonspecific_values], fontsize=10)

    # Labels
    ax.set_xlabel('BLAST Identity Threshold', fontsize=11, fontweight='bold')
    ax.set_ylabel('Max Non-Specific Residues per 600A²', fontsize=11, fontweight='bold')
    ax.set_title(pair.replace(':', ' × '), fontsize=12, fontweight='bold', pad=10)

    # Add text annotations with scores
    for i in range(len(max_nonspecific_values)):
        for j in range(len(specificity_thresholds)):
            score = heatmap_data[i, j]
            text_color = 'white' if score > 0.4 else 'black'
            text = f'{score:.3f}' if score > 0 else '0.000'
            ax.text(j, i, text, ha="center", va="center",
                   color=text_color, fontsize=9, fontweight='bold')

    # Grid
    ax.set_xticks(np.arange(len(specificity_thresholds))-.5, minor=True)
    ax.set_yticks(np.arange(len(max_nonspecific_values))-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=2)
    ax.tick_params(which="minor", size=0)

# Add colorbar
cbar = fig.colorbar(im, ax=axes, orientation='horizontal',
                    pad=0.08, shrink=0.8, aspect=30)
cbar.set_label('Bispecific Pair Score', fontsize=12, fontweight='bold')

# Adjust layout
plt.tight_layout(rect=[0, 0.05, 1, 0.96])

# Save figure
output_file = results_dir / f"specificity_heatmap_{latest_csv.stem.split('_')[-1]}.png"
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"✓ Heatmap saved to: {output_file}")

# Create summary bar chart
fig2, ax2 = plt.subplots(figsize=(14, 6))

# Count viable pairs (score > 0) for each parameter combination
summary = df.groupby(['max_nonspecific', 'specificity_threshold'])['has_patches'].sum().reset_index()
summary['param_label'] = (summary['max_nonspecific'].astype(str) + ' res, ' +
                          (summary['specificity_threshold'] * 100).astype(int).astype(str) + '%')

# Bar chart
x = np.arange(len(summary))
colors = [PALETTE['green'] if val == 3 else PALETTE['teal'] if val == 2 else
          PALETTE['blue'] if val == 1 else PALETTE['light_gray']
          for val in summary['has_patches']]

bars = ax2.bar(x, summary['has_patches'], color=colors, edgecolor='black', linewidth=1.5)

# Add value labels on bars
for i, (bar, val) in enumerate(zip(bars, summary['has_patches'])):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
            f'{int(val)}/3',
            ha='center', va='bottom', fontsize=9, fontweight='bold')

ax2.set_xlabel('Parameter Combination (Max Non-Specific, BLAST Threshold)',
               fontsize=12, fontweight='bold')
ax2.set_ylabel('Number of Viable Pairs (out of 3)', fontsize=12, fontweight='bold')
ax2.set_title('Patch Availability Across Parameter Space',
              fontsize=14, fontweight='bold', pad=15)
ax2.set_xticks(x)
ax2.set_xticklabels(summary['param_label'], rotation=45, ha='right', fontsize=9)
ax2.set_ylim(0, 3.5)
ax2.set_yticks([0, 1, 2, 3])
ax2.grid(axis='y', alpha=0.3, linestyle='--')

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=PALETTE['green'], edgecolor='black', label='3/3 pairs viable'),
    Patch(facecolor=PALETTE['teal'], edgecolor='black', label='2/3 pairs viable'),
    Patch(facecolor=PALETTE['blue'], edgecolor='black', label='1/3 pairs viable'),
    Patch(facecolor=PALETTE['light_gray'], edgecolor='black', label='0/3 pairs viable')
]
ax2.legend(handles=legend_elements, loc='upper right', frameon=True, fontsize=10)

plt.tight_layout()

# Save figure
output_file2 = results_dir / f"specificity_summary_{latest_csv.stem.split('_')[-1]}.png"
plt.savefig(output_file2, dpi=150, bbox_inches='tight')
print(f"✓ Summary chart saved to: {output_file2}")

print("\nVisualization complete!")
print(f"  - Heatmap: {output_file}")
print(f"  - Summary: {output_file2}")
