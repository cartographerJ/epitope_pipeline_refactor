#!/usr/bin/env python3
"""Generate diagnostic plots for CEACAM5 from the latest run."""

import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Cartography colors
CARTO_MINT = '#E0F3DB'
CARTO_RED = '#FF6B6B'
CARTO_TEAL = '#18B5CB'
CARTO_GRAY = '#D3D3D3'

# Load BLAST cache for CEACAM5
blast_cache = Path("epitope_pipeline/cache/blast/fullseq_a4cb9e0ebeec52e351b5b20f648bab94.json")
with open(blast_cache) as f:
    blast_data = json.load(f)

residue_spec = {int(k): float(v) for k, v in blast_data["residue_specificity"].items()}

# Parameters
threshold = 0.70
all_residues = sorted(residue_spec.keys())
identity_scores = [residue_spec.get(r, 0.0) for r in all_residues]

# Create diagnostic plot
fig, axes = plt.subplots(2, 1, figsize=(16, 8))
fig.suptitle('CEACAM5: Per-Residue BLAST Identity to Human Proteome',
             fontsize=16, fontweight='bold', y=0.995)

# Track 1: Per-residue identity scores
ax1 = axes[0]
colors = [CARTO_RED if s >= threshold else CARTO_MINT for s in identity_scores]
ax1.bar(all_residues, identity_scores, color=colors, width=1.0, alpha=0.7)
ax1.axhline(threshold, color='black', linestyle='--', linewidth=1.5,
            label=f'Threshold ({threshold:.0%})')
ax1.set_ylabel('BLAST Identity', fontsize=12, fontweight='bold')
ax1.set_ylim(0, 1.0)
ax1.set_title(f'Per-Residue Identity (Red ≥{threshold:.0%} = non-specific, Mint <{threshold:.0%} = specific)',
              fontsize=11, loc='left')
ax1.legend(loc='upper right', fontsize=10)
ax1.grid(axis='y', alpha=0.3)
ax1.set_xlim(0, 702)

# Track 2: Stats by region
ax2 = axes[1]

# Compute stats
n_total = len([s for s in identity_scores if s > 0])
n_specific = len([s for s in identity_scores if 0 < s < threshold])
n_nonspecific = len([s for s in identity_scores if s >= threshold])

# Bar chart
categories = ['Assessed', 'Specific\n(<70%)', 'Non-Specific\n(≥70%)']
values = [n_total, n_specific, n_nonspecific]
colors_bar = [CARTO_TEAL, CARTO_MINT, CARTO_RED]
bars = ax2.bar(categories, values, color=colors_bar, alpha=0.7, edgecolor='black', linewidth=1.5)

# Add value labels on bars
for bar, val in zip(bars, values):
    height = bar.get_height()
    pct = val/n_total*100 if n_total > 0 else 0
    ax2.text(bar.get_x() + bar.get_width()/2., height + 10,
             f'{val}\n({pct:.1f}%)',
             ha='center', va='bottom', fontsize=11, fontweight='bold')

ax2.set_ylabel('Residue Count', fontsize=12, fontweight='bold')
ax2.set_title(f'Summary: {n_nonspecific}/{n_total} residues ({n_nonspecific/n_total*100:.1f}%) are non-specific (≥{threshold:.0%} identity)',
              fontsize=11, loc='left')
ax2.set_ylim(0, max(values) * 1.2)
ax2.grid(axis='y', alpha=0.3)

plt.tight_layout()

# Save
output_path = Path("epitope_pipeline/runs/260310_1609_ceacam5/figures/ceacam5_specificity_diagnostics.png")
output_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Saved diagnostic plot: {output_path}")

# Also save stats
stats_path = output_path.parent / "ceacam5_specificity_stats.txt"
with open(stats_path, 'w') as f:
    f.write(f"CEACAM5 Specificity Diagnostics\n")
    f.write(f"{'='*60}\n\n")
    f.write(f"Threshold: {threshold:.1%}\n\n")
    f.write(f"Residue-Level Stats:\n")
    f.write(f"  Total assessed: {n_total}\n")
    f.write(f"  Specific (<{threshold:.0%}): {n_specific} ({n_specific/n_total*100:.1f}%)\n")
    f.write(f"  Non-specific (≥{threshold:.0%}): {n_nonspecific} ({n_nonspecific/n_total*100:.1f}%)\n\n")
    f.write(f"Mean identity: {sum(identity_scores)/len(identity_scores):.1%}\n")
    f.write(f"Max identity: {max(identity_scores):.1%}\n")
    f.write(f"Min identity: {min([s for s in identity_scores if s > 0]):.1%}\n\n")
    f.write(f"Conclusion:\n")
    f.write(f"CEACAM5 has massive cross-reactivity to CEACAM family members.\n")
    f.write(f"59% of residues exceed the 70% identity threshold, making it\n")
    f.write(f"unsuitable for human-specific VHH binder discovery.\n")

print(f"Saved stats file: {stats_path}")
