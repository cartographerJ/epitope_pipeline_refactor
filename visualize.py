"""
Step 10: Visualization — matplotlib figures for epitope analysis results.

Generates:
  - Per-target linear epitope map with topology, SASA, conservation,
    distance, and patch annotations
  - Multi-target scoring summary bar chart

Uses the Cartography palette and figure style consistent with ablang2.py.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from epitope_pipeline.config import (
    COLOR_EPITOPE_PATCH,
    COLOR_MISMATCH,
    COLOR_TRANSMEMBRANE,
    PALETTE,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plot_epitope_map(target, membrane, spatial_filter, surface_analysis,
                     conservation_result, specificity_result, scores,
                     output_path, title_suffix=None):
    """
    Generate a linear sequence map annotating the epitope analysis results.

    Tracks (top to bottom):
      1. Topology + Domains — topology coloring with labeled domain blocks
      2. Distance from membrane (line plot with threshold)
      3. SASA bar chart
      4. Cyno conservation — per-residue identity band
      5. Human specificity — per-patch BLAST result band (same style as cyno)
      6. Target epitope — final qualifying patches

    Args:
        target: TargetInfo object (with features for domain labels).
        membrane: MembraneAnnotation.
        spatial_filter: SpatialFilter.
        surface_analysis: SurfaceAnalysis.
        conservation_result: ConservationResult.
        specificity_result: SpecificityResult (or None).
        scores: List of EpitopeScore objects.
        output_path: Path to save PNG.
    """
    seq_len = target.sequence_length
    positions = np.arange(1, seq_len + 1)

    # ---- Pre-compute per-residue data ----

    # Use NaN for residues not in the structure so unresolved regions
    # appear as gaps rather than misleading zeros
    distances = np.array([
        spatial_filter.residue_distances.get(i, np.nan)
        for i in range(1, seq_len + 1)
    ])

    sasa_dict = surface_analysis.residue_sasa if surface_analysis else {}
    sasa_vals = np.array([
        sasa_dict.get(i, np.nan)
        for i in range(1, seq_len + 1)
    ])

    # Conservation per-residue
    conservation_colors = []
    if conservation_result:
        for i in range(1, seq_len + 1):
            cons = conservation_result.residue_conservation.get(i)
            if cons is True:
                conservation_colors.append(PALETTE["mint"])
            elif cons is False:
                conservation_colors.append(COLOR_MISMATCH)
            else:
                conservation_colors.append("#FFFFFF")
    else:
        conservation_colors = ["#FFFFFF"] * seq_len

    # Specificity per-residue (binary from ectodomain 3D patch screening)
    # True/None = specific (mint), False = non-specific (red)
    specificity_colors = []
    if specificity_result and specificity_result.residue_specificity:
        for i in range(1, seq_len + 1):
            val = specificity_result.residue_specificity.get(i)
            if val is False:
                specificity_colors.append(COLOR_MISMATCH)
            else:
                specificity_colors.append(PALETTE["mint"])
    else:
        specificity_colors = [PALETTE["mint"]] * seq_len

    # Qualifying patch residue set (for target epitope track)
    patch_residues = set()
    for s in scores:
        for r in s.patch.residue_numbers:
            patch_residues.add(r)

    # ---- Create 6-track figure ----
    fig, axes = plt.subplots(
        6, 1,
        figsize=(max(18, seq_len * 0.025), 11),
        gridspec_kw={"height_ratios": [0.7, 1.4, 1.4, 0.4, 0.4, 0.5]},
        sharex=True,
    )
    title = "{} ({}) \u2014 Epitope Analysis".format(target.gene_name, target.uniprot_id)
    if title_suffix:
        title += " \u2014 {}".format(title_suffix)
    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.98)

    # ==================================================================
    # Track 1: Domains (compact blocks with topology brackets above)
    # ==================================================================
    ax_dom = axes[0]
    ax_dom.set_xlim(0.5, seq_len + 0.5)
    ax_dom.set_ylim(0, 1)
    ax_dom.set_yticks([])
    ax_dom.set_ylabel("Domains", fontsize=11, rotation=0, ha="right", va="center")

    # --- Topology brackets above the domain blocks ---
    _draw_topology_brackets(ax_dom, target, membrane, seq_len)

    # --- Collect and render domain blocks ---
    domain_blocks = _collect_domain_blocks(target)

    # Draw largest first so smaller features appear on top
    domain_blocks.sort(key=lambda x: -(x[1] - x[0]))

    block_y = 0.05
    block_h = 0.55
    for start, end, label, fcolor in domain_blocks:
        span = end - start + 1
        rect = plt.Rectangle(
            (start - 0.5, block_y), span, block_h,
            facecolor=fcolor, edgecolor="#333333", linewidth=0.7,
            alpha=0.85, zorder=3,
        )
        ax_dom.add_patch(rect)
        center = (start + end) / 2.0
        text_color = "white" if fcolor in (COLOR_TRANSMEMBRANE, PALETTE["blue"]) else "black"
        fontsize = 7 if span < 80 else 9
        txt = ax_dom.text(
            center, block_y + block_h / 2, label,
            ha="center", va="center", fontsize=fontsize, fontweight="bold",
            color=text_color, zorder=4, clip_on=True,
        )
        # Right-only clip: let text overflow left but clip at domain right edge
        right_clip = plt.Rectangle(
            (center - 5000, block_y), 5000 + (end - center + 0.5), block_h,
            transform=ax_dom.transData, visible=False,
        )
        ax_dom.add_patch(right_clip)
        txt.set_clip_path(right_clip)

    # ==================================================================
    # Track 2: Distance from membrane (abs value, ECD teal, non-ECD gray)
    # ==================================================================
    ax_dist = axes[1]

    # Split into ECD vs non-ECD for distinct styling
    ecd_dists = np.full(seq_len, np.nan)
    nonecd_dists = np.full(seq_len, np.nan)
    for i in range(seq_len):
        resnum = i + 1
        d = spatial_filter.residue_distances.get(resnum)
        if d is None:
            continue
        topo = membrane.residue_topology.get(resnum)
        if topo == "extracellular":
            ecd_dists[i] = d
        else:
            nonecd_dists[i] = d

    # ECD fill — teal, prominent
    ax_dist.fill_between(positions, ecd_dists, alpha=0.3, color=PALETTE["teal"])
    # Non-ECD fill — light gray, clearly distinct from ECD
    ax_dist.fill_between(positions, nonecd_dists, alpha=0.08, color="#999999")
    # One continuous line through all residues
    ax_dist.plot(positions, distances, color=PALETTE["teal"], linewidth=0.5, alpha=0.5)

    # Threshold line — proximal uses max_distance_threshold, distal uses min
    if getattr(spatial_filter, "max_distance_threshold", None) is not None:
        thresh = spatial_filter.max_distance_threshold
    else:
        thresh = spatial_filter.min_distance_threshold
    ax_dist.axhline(
        y=thresh, color=PALETTE["dark_purple"],
        linestyle="--", linewidth=1, alpha=0.7,
    )
    ax_dist.set_ylabel("Dist from\nmembrane (A)", fontsize=11)
    ax_dist.tick_params(labelsize=10)

    # ==================================================================
    # Track 3: SASA
    # ==================================================================
    ax_sasa = axes[2]
    valid_sasa = ~np.isnan(sasa_vals)
    if np.any(valid_sasa):
        ax_sasa.bar(
            positions[valid_sasa], sasa_vals[valid_sasa], width=1.0,
            color=PALETTE["teal"], alpha=0.6, linewidth=0,
        )

    ax_sasa.set_ylabel("SASA (A\u00B2)", fontsize=11)
    ax_sasa.set_ylim(bottom=0)
    ax_sasa.tick_params(labelsize=10)

    # ==================================================================
    # Track 4: Cyno conservation (per-residue band)
    # ==================================================================
    ax_cons = axes[3]
    for i, color in enumerate(conservation_colors):
        ax_cons.axvspan(i + 0.5, i + 1.5, color=color, alpha=0.8)
    ax_cons.set_yticks([])
    ax_cons.set_ylabel("Cyno\nConserv.", fontsize=11, rotation=0, ha="right", va="center")

    # Label when no cyno ortholog exists
    if not conservation_result and not target.cyno_sequence:
        ax_cons.text(
            seq_len / 2.0, 0.5, "No cyno ortholog found",
            ha="center", va="center", fontsize=11, fontstyle="italic",
            color="#888888", transform=ax_cons.get_xaxis_transform(),
        )

    # ==================================================================
    # Track 5: Human specificity (per-residue band, same style as cyno)
    # ==================================================================
    ax_spec = axes[4]
    for i, color in enumerate(specificity_colors):
        ax_spec.axvspan(i + 0.5, i + 1.5, color=color, alpha=0.8)
    ax_spec.set_yticks([])
    ax_spec.set_ylabel("Human\nSpecific.", fontsize=11, rotation=0, ha="right", va="center")

    # ==================================================================
    # Track 6: Target epitope (final qualifying patches)
    # ==================================================================
    ax_patch = axes[5]

    for i in range(seq_len):
        resnum = i + 1
        if resnum in patch_residues:
            ax_patch.axvspan(i + 0.5, i + 1.5, color=COLOR_EPITOPE_PATCH, alpha=0.9)

    ax_patch.set_yticks([])
    ax_patch.set_ylabel("Target\nEpitope", fontsize=11, rotation=0, ha="right", va="center")
    ax_patch.set_xlabel("Residue Position", fontsize=12)

    # Patch labels (score + rank)
    if scores:
        for s in scores:
            center = np.mean(s.patch.residue_numbers)
            label = "Patch {} \u2014 {:.0f} A\u00B2 | score {:.2f}".format(
                s.rank, s.patch_area_a2, s.composite_score)
            ax_patch.text(
                center, 0.5, label,
                ha="center", va="center", fontsize=8.5, fontweight="bold",
                transform=ax_patch.get_xaxis_transform(),
            )
    else:
        ax_patch.text(
            seq_len / 2.0, 0.5, "No qualifying epitope patches",
            ha="center", va="center", fontsize=11, fontstyle="italic",
            color="#888888", transform=ax_patch.get_xaxis_transform(),
        )

    # ==================================================================
    # Legend
    # ==================================================================
    legend_patches = [
        mpatches.Patch(color=PALETTE["teal"], label="InterPro Domain"),
        mpatches.Patch(color=PALETTE["blue"], label="UniProt Domain"),
        mpatches.Patch(color=COLOR_TRANSMEMBRANE, label="Transmembrane"),
        mpatches.Patch(color=PALETTE["mint"], label="Conserved / Specific"),
        mpatches.Patch(color=COLOR_MISMATCH, label="Mismatch / Non-specific"),
        mpatches.Patch(color=COLOR_EPITOPE_PATCH, label="Target Epitope"),
    ]
    fig.legend(
        handles=legend_patches, loc="lower center", ncol=3, fontsize=10,
        frameon=True, fancybox=True,
    )

    plt.tight_layout(rect=[0, 0.06, 1, 0.96])

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Wrote epitope map: %s", output_path.name)


# ---------------------------------------------------------------------------
# Domain track helpers
# ---------------------------------------------------------------------------

def _draw_topology_brackets(ax, target, membrane, seq_len):
    """
    Draw thin bracket annotations for Extracellular / Cytoplasmic spans
    at the top of the domain track.

    For GPI-anchored proteins that lack UniProt "Topological domain" features,
    infers a single "Extracellular" bracket spanning the mature protein.
    """
    bracket_y = 0.78
    tick_h = 0.08
    label_y = 0.82
    color = "#555555"

    # Collect topology spans from UniProt features
    topo_spans = []
    for feat in target.features:
        if feat.get("type") != "Topological domain":
            continue
        desc = feat.get("description", "")
        start = feat.get("start")
        end = feat.get("end")
        if start is None or end is None:
            continue

        if "extracellular" in desc.lower():
            topo_spans.append(("Extracellular", start, end))
        elif "cytoplasmic" in desc.lower():
            topo_spans.append(("Cytoplasmic", start, end))

    # Fallback for GPI-anchored proteins (no topological domain annotations)
    if not topo_spans and membrane and membrane.topology_type == "gpi_anchored":
        # Signal peptide end (mature protein start)
        sp_end = 1
        for feat in target.features:
            if feat.get("type") in ("Signal peptide", "Signal") and feat.get("end"):
                sp_end = feat["end"] + 1
                break
        topo_spans.append(("Extracellular", sp_end, seq_len))

    for label, start, end in topo_spans:
        # Horizontal line
        ax.plot([start, end], [bracket_y, bracket_y],
                color=color, linewidth=1.0, clip_on=False, zorder=5)
        # Downward ticks at ends
        ax.plot([start, start], [bracket_y, bracket_y - tick_h],
                color=color, linewidth=1.0, clip_on=False, zorder=5)
        ax.plot([end, end], [bracket_y, bracket_y - tick_h],
                color=color, linewidth=1.0, clip_on=False, zorder=5)
        # Centered label above
        center = (start + end) / 2.0
        ax.text(center, label_y, label,
                ha="center", va="bottom", fontsize=9, color=color,
                fontweight="bold", zorder=5)


def _collect_domain_blocks(target):
    """
    Collect all domain features to render as blocks.

    Returns list of (start, end, label, color) tuples.
    Filters out overlapping sub-domains and deduplicates InterPro vs UniProt.
    """
    _FEAT_COLORS = {
        "Signal peptide": "#777777",
        "Signal": "#777777",
        "Transmembrane": COLOR_TRANSMEMBRANE,
        "Lipidation": COLOR_TRANSMEMBRANE,
        "Domain": PALETTE["blue"],
        "Region": "#BBBBBB",
        "InterPro": PALETTE["teal"],
        "Chain": PALETTE["blue"],
    }
    _ABBREV = {
        "Signal peptide": "SP", "Signal": "SP", "Transmembrane": "TM",
    }

    # Include Chain features only for cleaved precursors (2+ chains)
    chain_features = [f for f in target.features
                      if f.get("type") == "Chain"
                      and f.get("start") is not None
                      and f.get("end") is not None]
    include_chains = len(chain_features) >= 2

    raw_blocks = []
    for feat in target.features:
        ftype = feat.get("type", "")
        start = feat.get("start")
        end = feat.get("end")
        desc = feat.get("description", "")
        if start is None or end is None:
            continue
        span = end - start + 1

        # Skip topology (handled by brackets) and transit peptide
        if ftype in ("Topological domain", "Intramembrane",
                      "Transit peptide"):
            continue
        # Include Chain only for cleaved precursors
        if ftype == "Chain" and not include_chains:
            continue

        if ftype in _ABBREV:
            label = _ABBREV[ftype]
        elif ftype == "Lipidation" and "gpi" in desc.lower():
            label = "GPI"
        elif ftype in ("Domain", "InterPro", "Chain") and desc:
            label = desc
        elif ftype == "Region" and desc and span > 50:
            label = desc
        else:
            continue

        # Skip very small features (except TM, SP, and GPI)
        if span < 15 and ftype not in ("Transmembrane", "Signal peptide",
                                        "Signal", "Lipidation"):
            continue

        if len(label) > 22:
            label = label[:20] + ".."

        fcolor = _FEAT_COLORS.get(ftype, "#CCCCCC")
        raw_blocks.append((start, end, label, fcolor, ftype))

    # Filter overlapping entries:
    # - UniProt features (TM, SP, Domain, Region) always kept
    # - InterPro entries removed if >70% overlaps with any UniProt feature
    #   or a larger InterPro entry
    uniprot_blocks = [(s, e, l, c, t) for s, e, l, c, t in raw_blocks
                      if t != "InterPro"]
    interpro_blocks = [(s, e, l, c, t) for s, e, l, c, t in raw_blocks
                       if t == "InterPro"]

    # Keep all UniProt blocks
    filtered = [(s, e, l, c) for s, e, l, c, t in uniprot_blocks]

    # Filter InterPro: skip if overlaps significantly with any UniProt feature
    # or a larger InterPro entry already kept.
    # Use 50% threshold vs UniProt (catch TM-juxtamembrane etc.)
    # and 70% threshold vs other InterPro entries.
    all_keepers = list(uniprot_blocks)
    interpro_blocks.sort(key=lambda x: -(x[1] - x[0]))

    for s1, e1, l1, c1, t1 in interpro_blocks:
        span1 = e1 - s1 + 1
        contained = False
        for s2, e2, l2, c2, t2 in all_keepers:
            overlap = max(0, min(e1, e2) - max(s1, s2) + 1)
            threshold = 0.50 if t2 != "InterPro" else 0.70
            if overlap / span1 > threshold:
                contained = True
                break
        if not contained:
            filtered.append((s1, e1, l1, c1))
            all_keepers.append((s1, e1, l1, c1, t1))

    return filtered


def plot_scoring_summary(all_scores, target_metrics, targets, output_path):
    """
    Multi-target comparison bar chart.

    Shows for each target:
      - Total qualifying epitope surface area
      - Number of patches
      - Best composite score

    Args:
        all_scores: Dict {uniprot_id: [EpitopeScore, ...]}.
        target_metrics: Dict {uniprot_id: metric_dict}.
        targets: List of TargetInfo objects.
        output_path: Path to save PNG.
    """
    gene_names = [t.gene_name for t in targets]
    n_targets = len(targets)

    total_areas = []
    n_patches_list = []
    best_scores = []

    for t in targets:
        uid = t.uniprot_id
        metric = target_metrics.get(uid, {})
        total_areas.append(metric.get("total_epitope_area_a2", 0.0))
        n_patches_list.append(metric.get("n_patches", 0))
        best_scores.append(metric.get("best_score", 0.0))

    fig, axes = plt.subplots(1, 3, figsize=(4 * n_targets, 5))
    fig.suptitle("Epitope Space Summary", fontsize=16, fontweight="bold")

    x = np.arange(n_targets)
    bar_width = 0.6

    # Panel 1: Total epitope area
    ax1 = axes[0]
    bars1 = ax1.bar(x, total_areas, bar_width, color=PALETTE["teal"], alpha=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(gene_names, fontsize=12)
    ax1.set_ylabel("Total Epitope Area (A²)", fontsize=12)
    ax1.set_title("Epitope Surface", fontsize=13)
    for bar, val in zip(bars1, total_areas):
        if val > 0:
            ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                     "{:.0f}".format(val), ha="center", va="bottom", fontsize=11)

    # Panel 2: Number of patches
    ax2 = axes[1]
    bars2 = ax2.bar(x, n_patches_list, bar_width, color=PALETTE["green"], alpha=0.8)
    ax2.set_xticks(x)
    ax2.set_xticklabels(gene_names, fontsize=12)
    ax2.set_ylabel("Number of Patches", fontsize=12)
    ax2.set_title("Qualifying Patches", fontsize=13)
    ax2.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    for bar, val in zip(bars2, n_patches_list):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                 str(val), ha="center", va="bottom", fontsize=11)

    # Panel 3: Best composite score
    ax3 = axes[2]
    bars3 = ax3.bar(x, best_scores, bar_width, color=PALETTE["dark_purple"], alpha=0.8)
    ax3.set_xticks(x)
    ax3.set_xticklabels(gene_names, fontsize=12)
    ax3.set_ylabel("Composite Score", fontsize=12)
    ax3.set_title("Best Patch Score", fontsize=13)
    ax3.set_ylim(0, 1.0)
    for bar, val in zip(bars3, best_scores):
        if val > 0:
            ax3.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                     "{:.3f}".format(val), ha="center", va="bottom", fontsize=11)

    plt.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Wrote scoring summary: %s", output_path.name)
