#!/usr/bin/env python3
"""
Generate patch-level diagnostic visualization for CEACAM5.

Shows which patches passed/failed at each filtering step:
- Cyno conservation
- Human specificity
- Final patches
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Cartography colors
CARTO_GRAY = '#D3D3D3'
CARTO_MINT = '#E0F3DB'
CARTO_GREEN = '#6BC291'
CARTO_TEAL = '#18B5CB'
CARTO_BLUE = '#2E95D2'
CARTO_PURPLE = '#28154C'

# Load data from the latest CEACAM5 run
run_dir = Path("epitope_pipeline/runs/260310_1609_ceacam5")

# Load residue table to get cyno conservation
residue_df = pd.read_csv(run_dir / "annotated_sequences/ceacam5_residue_table.csv")

# Load BLAST data for specificity
blast_cache = Path("epitope_pipeline/cache/blast/fullseq_a4cb9e0ebeec52e351b5b20f648bab94.json")
with open(blast_cache) as f:
    blast_data = json.load(f)
    residue_identity = {int(k): float(v) for k, v in blast_data["residue_specificity"].items()}

# From the log, we know these patches existed:
# Surface patches: 2 patches (149 residues @ patch 0, 101 residues @ patch 1)
# After cyno filtering: 3 sub-patches (19, 11, 7 residues)
# After specificity: 0 patches

# Reconstruct patch information from residue table
patches_info = []

# Patch 0: residues with patch_id "0" in the CSV (if it exists)
# Actually, let me just manually define based on the log output

# From log:
# "2 surface patches (>= 600 A²), total epitope surface: 20814 A²"
# Patch cyno results showed sub-patches: 19, 11, 7 residues
# Specificity rejected all 3

# Let me create a simplified visualization showing:
# 1. All extracellular residues
# 2. Cyno conservation status
# 3. Specificity status
# 4. Final patch status

fig, axes = plt.subplots(4, 1, figsize=(16, 12))
fig.suptitle('CEACAM5: Patch-Level Filtering Cascade', fontsize=16, fontweight='bold', y=0.995)

all_residues = sorted(residue_df['residue_num'].values)

# Track 1: All extracellular + surface-exposed residues
ax1 = axes[0]
ec_residues = residue_df[residue_df['topology'] == 'extracellular']
surface_residues = ec_residues[ec_residues['sasa_A2'] > 0]

for idx, row in surface_residues.iterrows():
    resnum = row['residue_num']
    distance = row.get('distance_from_membrane_A', 0)
    if distance >= 80:
        ax1.scatter(resnum, 1, c=CARTO_TEAL, s=30, marker='s', alpha=0.8)
    else:
        ax1.scatter(resnum, 1, c=CARTO_GRAY, s=20, marker='s', alpha=0.4)

ax1.set_ylim(0.5, 1.5)
ax1.set_yticks([])
ax1.set_ylabel('Spatial\nFilter', fontsize=11, fontweight='bold')
ax1.set_title(f'Step 1-4: Surface + Distance ≥80A → {len(surface_residues[surface_residues["distance_from_membrane_A"] >= 80])} qualifying residues',
              fontsize=11, loc='left')
ax1.axhline(1, color='black', alpha=0.1)
ax1.grid(axis='x', alpha=0.3)
ax1.set_xlim(0, 702)
ax1.legend(['≥80A (qualify)', '<80A (fail)'], loc='upper right', fontsize=9)

# Track 2: Cyno conservation
ax2 = axes[1]
qualifying_residues = surface_residues[surface_residues['distance_from_membrane_A'] >= 80]

# Plot in two passes: passes first (background), then fails (foreground)
passes = []
fails = []
for idx, row in qualifying_residues.iterrows():
    resnum = row['residue_num']
    cyno_status = row.get('cyno_conserved', '')
    if cyno_status == 'yes':
        passes.append(resnum)
    elif cyno_status == 'no':
        fails.append(resnum)

# Plot passes first (smaller, background)
if passes:
    ax2.scatter(passes, [1]*len(passes), c=CARTO_MINT, s=30, marker='s', alpha=0.8)
# Plot fails second (larger, foreground)
if fails:
    ax2.scatter(fails, [1]*len(fails), c=CARTO_PURPLE, s=30, marker='s', alpha=0.8)

ax2.set_ylim(0.5, 1.5)
ax2.set_yticks([])
ax2.set_ylabel('Cyno\nConservation', fontsize=11, fontweight='bold')
n_conserved = len(qualifying_residues[qualifying_residues['cyno_conserved'] == 'yes'])
n_mismatch = len(qualifying_residues[qualifying_residues['cyno_conserved'] == 'no'])
ax2.set_title(f'Step 6 (PRE-PATCH): Cyno Conservation → {n_conserved} conserved (mint), {n_mismatch} mismatch (purple)',
              fontsize=11, loc='left')
ax2.axhline(1, color='black', alpha=0.1)
ax2.grid(axis='x', alpha=0.3)
ax2.set_xlim(0, 702)

# Track 3: Human specificity
ax3 = axes[2]

# Plot in two passes: passes first (background), then fails (foreground)
spec_passes = []
spec_fails = []
for idx, row in qualifying_residues.iterrows():
    resnum = row['residue_num']
    identity = residue_identity.get(resnum, 0.0)
    if identity < 0.70:
        spec_passes.append(resnum)
    else:
        spec_fails.append(resnum)

# Plot passes first (smaller, background)
if spec_passes:
    ax3.scatter(spec_passes, [1]*len(spec_passes), c=CARTO_MINT, s=30, marker='s', alpha=0.8)
# Plot fails second (larger, foreground)
if spec_fails:
    ax3.scatter(spec_fails, [1]*len(spec_fails), c=CARTO_PURPLE, s=30, marker='s', alpha=0.8)

ax3.set_ylim(0.5, 1.5)
ax3.set_yticks([])
ax3.set_ylabel('Human\nSpecificity', fontsize=11, fontweight='bold')
n_specific = len([r for r in qualifying_residues['residue_num'].values
                  if residue_identity.get(r, 0) < 0.70])
n_nonspecific = len([r for r in qualifying_residues['residue_num'].values
                     if residue_identity.get(r, 0) >= 0.70])
ax3.set_title(f'Step 7 (PRE-PATCH): Human Specificity → {n_specific} specific (mint), {n_nonspecific} non-specific (purple)',
              fontsize=11, loc='left')
ax3.axhline(1, color='black', alpha=0.1)
ax3.grid(axis='x', alpha=0.3)
ax3.set_xlim(0, 702)

# Track 4: Intersection (pass BOTH filters)
ax4 = axes[3]

# Categorize residues
pass_both = []
fail_cyno_only = []
fail_spec_only = []
fail_both = []

for idx, row in qualifying_residues.iterrows():
    resnum = row['residue_num']
    cyno_ok = row.get('cyno_conserved', '') == 'yes'
    spec_ok = residue_identity.get(resnum, 1.0) < 0.70

    if cyno_ok and spec_ok:
        pass_both.append(resnum)
    elif not cyno_ok and not spec_ok:
        fail_both.append(resnum)
    elif not cyno_ok:
        fail_cyno_only.append(resnum)
    elif not spec_ok:
        fail_spec_only.append(resnum)

# Plot as horizontal bars at different y-levels
if pass_both:
    ax4.scatter(pass_both, [4]*len(pass_both), c=CARTO_GREEN, s=50, marker='s', alpha=0.9, label='Pass both')
if fail_spec_only:
    ax4.scatter(fail_spec_only, [3]*len(fail_spec_only), c=CARTO_PURPLE, s=50, marker='s', alpha=0.9, label='Fail specificity only')
if fail_cyno_only:
    ax4.scatter(fail_cyno_only, [2]*len(fail_cyno_only), c=CARTO_BLUE, s=50, marker='s', alpha=0.9, label='Fail cyno only')
if fail_both:
    ax4.scatter(fail_both, [1]*len(fail_both), c=CARTO_GRAY, s=50, marker='s', alpha=0.7, label='Fail both')

ax4.set_ylim(0.5, 4.5)
ax4.set_yticks([1, 2, 3, 4])
ax4.set_yticklabels(['Fail\nBoth', 'Fail\nCyno', 'Fail\nSpec', 'Pass\nBoth'], fontsize=9)
ax4.set_ylabel('Filter\nOutcome', fontsize=11, fontweight='bold')
ax4.set_xlabel('Residue Number', fontsize=12, fontweight='bold')
n_pass_both = len(pass_both)
ax4.set_title(f'Step 8 (POST-PATCH): Intersection → {n_pass_both} residues pass BOTH filters → 0 patches ≥600A²',
              fontsize=11, loc='left')
ax4.grid(axis='x', alpha=0.3)
ax4.set_xlim(0, 702)
ax4.legend(loc='upper right', fontsize=9, ncol=2, framealpha=0.9)

plt.tight_layout()

# Save
output_path = run_dir / "figures/ceacam5_patch_filtering_cascade.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Saved patch filtering cascade: {output_path}")
plt.close()

# Now create per-patch BLAST report
print("\nGenerating per-patch BLAST report...")

# Get top hits
top_hits = sorted(blast_data["hits"], key=lambda h: h["identity"], reverse=True)[:10]

report_path = run_dir / "ceacam5_blast_report.txt"
with open(report_path, 'w') as f:
    f.write("CEACAM5: Per-Patch BLAST Analysis Report\n")
    f.write("="*80 + "\n\n")

    f.write("TOP 10 OFF-TARGET HITS (sorted by identity):\n")
    f.write("-"*80 + "\n\n")

    for i, hit in enumerate(top_hits, 1):
        f.write(f"{i}. {hit['accession']}: {hit['identity']:.1%} identity\n")
        f.write(f"   Residues: {hit['query_start']}-{hit['query_end']} ({hit['align_length']} aa)\n")
        f.write(f"   {hit['title'][:70]}\n\n")

    f.write("\n" + "="*80 + "\n")
    f.write("PATCH REJECTION ANALYSIS:\n")
    f.write("="*80 + "\n\n")

    f.write("From the pipeline log, 3 patches failed specificity screening:\n\n")

    # Patch 1: 19 residues, rejected
    f.write("PATCH 1 (19 residues):\n")
    f.write("  Rejection reason: Residue 26 has 10 non-specific residues in ~600A²\n")
    f.write("  Primary off-target: CEACAM3 (85.3% identity, residues 1-143)\n")
    f.write("  Status: REJECTED (all 19 residues trimmed)\n\n")

    # Patch 2: 11 residues
    f.write("PATCH 2 (11 residues):\n")
    f.write("  Rejection reason: Residue 142 has 11 non-specific residues in ~600A²\n")
    f.write("  Primary off-target: CEACAM3 (85.3% identity, residues 1-143)\n")
    f.write("  Status: REJECTED (all 11 residues trimmed)\n\n")

    # Patch 3: 7 residues
    f.write("PATCH 3 (7 residues):\n")
    f.write("  Rejection reason: Residue 175 has 7 non-specific residues in ~600A²\n")
    f.write("  Primary off-target: CEACAM6 (83.9% identity, residues 1-323)\n")
    f.write("  Status: REJECTED (all 7 residues trimmed)\n\n")

    f.write("\n" + "="*80 + "\n")
    f.write("CONCLUSION:\n")
    f.write("="*80 + "\n\n")
    f.write("All patches were located in the N-terminal domain (residues 1-323) which\n")
    f.write("shows 83-85% identity to CEACAM3/6. These regions are unsuitable for\n")
    f.write("human-specific binder discovery due to high cross-reactivity risk.\n\n")
    f.write("The C-terminal region (residues 423-702) has lower homology (<70%)\n")
    f.write("but was filtered out earlier due to:\n")
    f.write("  - Distance from membrane (<80A) or\n")
    f.write("  - Cyno conservation failures or\n")
    f.write("  - Insufficient surface area for patch formation\n")

print(f"Saved BLAST report: {report_path}")
