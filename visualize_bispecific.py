"""
Bispecific visualization — dual-target epitope maps and pair summary charts.
"""

import logging
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

from epitope_pipeline.config import (
    COLOR_EPITOPE_PATCH,
    COLOR_MISMATCH,
    COLOR_TRANSMEMBRANE,
    PALETTE,
)
from epitope_pipeline.visualize import (
    _draw_topology_brackets,
    _collect_domain_blocks,
)

logger = logging.getLogger("epitope_pipeline")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plot_bispecific_epitope_map(orientation, output_path):
    """
    Two-panel stacked figure: distal target on top, proximal on bottom.

    Each panel has 6 tracks (matching the single-target epitope map):
      1. Domains — topology brackets + domain blocks
      2. Distance — line plot with zone-specific threshold
      3. SASA — per-residue bar chart
      4. Cyno conservation — per-residue band
      5. Human specificity — per-residue band
      6. Qualifying patches — green patches with labels

    Args:
        orientation: OrientationResult with distal_zone and proximal_zone.
        output_path: Path to save PNG.
    """
    dz = orientation.distal_zone
    pz = orientation.proximal_zone

    fig = plt.figure(figsize=(max(18, max(
        dz.target.sequence_length, pz.target.sequence_length) * 0.025), 20))

    fig.suptitle(
        "BISPECIFIC: {} (distal) \u00d7 {} (proximal)".format(
            dz.target.gene_name, pz.target.gene_name),
        fontsize=16, fontweight="bold", y=0.98,
    )

    # Add pair score annotation
    fig.text(
        0.98, 0.98,
        "Pair score: {:.3f}".format(orientation.orientation_score),
        ha="right", va="top", fontsize=13, fontweight="bold",
        color=PALETTE["dark_purple"],
    )

    # Two panels with GridSpec (top = distal, bottom = proximal)
    gs = gridspec.GridSpec(2, 1, hspace=0.30, top=0.95, bottom=0.05)

    # Distal panel (top)
    gs_distal = gridspec.GridSpecFromSubplotSpec(
        6, 1, subplot_spec=gs[0],
        height_ratios=[0.7, 1.4, 1.4, 0.4, 0.4, 0.5], hspace=0.10,
    )
    _draw_zone_panel(fig, gs_distal, dz, zone="distal")

    # Proximal panel (bottom)
    gs_proximal = gridspec.GridSpecFromSubplotSpec(
        6, 1, subplot_spec=gs[1],
        height_ratios=[0.7, 1.4, 1.4, 0.4, 0.4, 0.5], hspace=0.10,
    )
    _draw_zone_panel(fig, gs_proximal, pz, zone="proximal")

    # Legend
    legend_patches = [
        mpatches.Patch(color=PALETTE["teal"], label="InterPro Domain"),
        mpatches.Patch(color=PALETTE["blue"], label="UniProt Domain"),
        mpatches.Patch(color=COLOR_TRANSMEMBRANE, label="Transmembrane"),
        mpatches.Patch(color=PALETTE["mint"], label="Conserved / Specific"),
        mpatches.Patch(color=COLOR_MISMATCH, label="Mismatch / Non-specific"),
        mpatches.Patch(color=COLOR_EPITOPE_PATCH, label="Epitope Patch"),
    ]
    fig.legend(
        handles=legend_patches, loc="lower center", ncol=3, fontsize=10,
        frameon=True, fancybox=True,
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Wrote bispecific map: %s", output_path.name)


def plot_bispecific_summary(pair_results, output_path):
    """
    Multi-pair comparison bar chart.

    Three panels:
      1. Final pair scores (with flexibility bonus marker)
      2. Best distal score per pair
      3. Best proximal score per pair

    Args:
        pair_results: List of BispecificPairResult sorted by final_pair_score.
        output_path: Path to save PNG.
    """
    if not pair_results:
        return

    n_pairs = len(pair_results)
    pair_labels = [
        "{} \u00d7 {}".format(pr.target_a.gene_name, pr.target_b.gene_name)
        for pr in pair_results
    ]

    final_scores = [pr.final_pair_score for pr in pair_results]
    distal_scores = [pr.best_orientation.distal_zone.best_score
                     if pr.best_orientation.distal_zone else 0.0
                     for pr in pair_results]
    proximal_scores = [pr.best_orientation.proximal_zone.best_score
                       if pr.best_orientation.proximal_zone else 0.0
                       for pr in pair_results]
    dual_valid = [pr.both_valid for pr in pair_results]

    # Vertical layout: 3 rows, 1 column — narrow horizontal bars for slides
    row_h = max(0.22 * n_pairs, 1.5)
    fig, axes = plt.subplots(3, 1, figsize=(2.5, row_h * 3 + 1.0))
    fig.suptitle("Bispecific Pair Summary", fontsize=7, fontweight="bold")

    y = np.arange(n_pairs)
    bar_height = 0.55

    # Panel 1: Final pair score
    ax1 = axes[0]
    bars1 = ax1.barh(y, final_scores, bar_height, color=PALETTE["teal"], alpha=0.8)
    ax1.set_yticks(y)
    ax1.set_yticklabels(pair_labels, fontsize=5)
    ax1.set_xlabel("Final Pair Score", fontsize=6)
    ax1.set_title("Pair Score", fontsize=6.5)
    ax1.set_xlim(0, 1.15)
    ax1.tick_params(axis="x", labelsize=5)
    ax1.invert_yaxis()
    for i, (bar, val) in enumerate(zip(bars1, final_scores)):
        label = "{:.3f}".format(val)
        if dual_valid[i]:
            label += " *"
        ax1.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                 label, ha="left", va="center", fontsize=4.5)

    # Panel 2: Best distal score
    ax2 = axes[1]
    bars2 = ax2.barh(y, distal_scores, bar_height, color=PALETTE["green"], alpha=0.8)
    ax2.set_yticks(y)
    ax2.set_yticklabels(pair_labels, fontsize=5)
    ax2.set_xlabel("Composite Score", fontsize=6)
    ax2.set_title("Best Distal Score", fontsize=6.5)
    ax2.set_xlim(0, 1.15)
    ax2.tick_params(axis="x", labelsize=5)
    ax2.invert_yaxis()
    for bar, val in zip(bars2, distal_scores):
        if val > 0:
            ax2.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                     "{:.3f}".format(val), ha="left", va="center", fontsize=4.5)

    # Panel 3: Best proximal score
    ax3 = axes[2]
    bars3 = ax3.barh(y, proximal_scores, bar_height, color=PALETTE["blue"], alpha=0.8)
    ax3.set_yticks(y)
    ax3.set_yticklabels(pair_labels, fontsize=5)
    ax3.set_xlabel("Composite Score", fontsize=6)
    ax3.set_title("Best Proximal Score", fontsize=6.5)
    ax3.set_xlim(0, 1.15)
    ax3.tick_params(axis="x", labelsize=5)
    ax3.invert_yaxis()
    for bar, val in zip(bars3, proximal_scores):
        if val > 0:
            ax3.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                     "{:.3f}".format(val), ha="left", va="center", fontsize=4.5)

    # Footnote for dual-valid marker
    if any(dual_valid):
        fig.text(0.5, 0.003, "* = dual-valid (both orientations, 20% bonus)",
                 ha="center", fontsize=4.5, fontstyle="italic", color="#666666")

    plt.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Wrote bispecific summary: %s", output_path.name)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _draw_zone_panel(fig, gs, zone_result, zone):
    """
    Draw a 6-track panel for one target in one zone.

    Tracks (matching single-target epitope map):
      1. Domains
      2. Distance (with zone-specific threshold)
      3. SASA bar chart
      4. Cyno conservation band
      5. Human specificity band
      6. Qualifying patches
    """
    target = zone_result.target
    membrane = zone_result.membrane
    spatial = zone_result.spatial_filter
    surface = zone_result.surface_analysis
    conservation = zone_result.conservation_result
    specificity = zone_result.specificity_result
    scores = zone_result.scores

    seq_len = target.sequence_length
    positions = np.arange(1, seq_len + 1)

    # Pre-compute data
    distances = np.array([
        spatial.residue_distances.get(i, np.nan)
        for i in range(1, seq_len + 1)
    ])

    sasa_dict = surface.residue_sasa if surface else {}
    sasa_vals = np.array([
        sasa_dict.get(i, np.nan)
        for i in range(1, seq_len + 1)
    ])

    conservation_colors = []
    if conservation:
        for i in range(1, seq_len + 1):
            cons = conservation.residue_conservation.get(i)
            if cons is True:
                conservation_colors.append(PALETTE["mint"])
            elif cons is False:
                conservation_colors.append(COLOR_MISMATCH)
            else:
                conservation_colors.append("#FFFFFF")
    else:
        conservation_colors = ["#FFFFFF"] * seq_len

    specificity_colors = []
    if specificity and specificity.residue_specificity:
        for i in range(1, seq_len + 1):
            val = specificity.residue_specificity.get(i)
            if val is False:
                specificity_colors.append(COLOR_MISMATCH)
            else:
                specificity_colors.append(PALETTE["mint"])
    else:
        specificity_colors = [PALETTE["mint"]] * seq_len

    patch_residues = set()
    for s in scores:
        for r in s.patch.residue_numbers:
            patch_residues.add(r)

    # Create axes (6 tracks)
    ax_dom = fig.add_subplot(gs[0])
    ax_dist = fig.add_subplot(gs[1], sharex=ax_dom)
    ax_sasa = fig.add_subplot(gs[2], sharex=ax_dom)
    ax_cons = fig.add_subplot(gs[3], sharex=ax_dom)
    ax_spec = fig.add_subplot(gs[4], sharex=ax_dom)
    ax_patch = fig.add_subplot(gs[5], sharex=ax_dom)

    # Zone label
    if zone == "distal":
        zone_label = "DISTAL \u2265{:.0f}\u00c5".format(
            spatial.min_distance_threshold)
        zone_color = PALETTE["teal"]
    else:
        zone_label = "PROXIMAL \u2264{:.0f}\u00c5".format(
            spatial.max_distance_threshold)
        zone_color = PALETTE["blue"]

    # Target name as panel title
    ax_dom.set_title(
        "{} ({}) \u2014 {}".format(
            target.gene_name, target.uniprot_id, zone_label),
        fontsize=13, fontweight="bold", color=zone_color,
        loc="left",
    )

    # ==================================================================
    # Track 1: Domains
    # ==================================================================
    ax_dom.set_xlim(0.5, seq_len + 0.5)
    ax_dom.set_ylim(0, 1)
    ax_dom.set_yticks([])
    ax_dom.set_ylabel("Domains", fontsize=10, rotation=0, ha="right", va="center")

    _draw_topology_brackets(ax_dom, target, membrane, seq_len)

    domain_blocks = _collect_domain_blocks(target)
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
        right_clip = plt.Rectangle(
            (center - 5000, block_y), 5000 + (end - center + 0.5), block_h,
            transform=ax_dom.transData, visible=False,
        )
        ax_dom.add_patch(right_clip)
        txt.set_clip_path(right_clip)

    # ==================================================================
    # Track 2: Distance
    # ==================================================================
    ecd_dists = np.full(seq_len, np.nan)
    nonecd_dists = np.full(seq_len, np.nan)
    for i in range(seq_len):
        resnum = i + 1
        d = spatial.residue_distances.get(resnum)
        if d is None:
            continue
        topo = membrane.residue_topology.get(resnum)
        if topo == "extracellular":
            ecd_dists[i] = d
        else:
            nonecd_dists[i] = d

    ax_dist.fill_between(positions, ecd_dists, alpha=0.3, color=PALETTE["teal"])
    ax_dist.fill_between(positions, nonecd_dists, alpha=0.08, color="#999999")
    ax_dist.plot(positions, distances, color=PALETTE["teal"], linewidth=0.5, alpha=0.5)

    # Zone-specific threshold line
    if zone == "distal":
        thresh = spatial.min_distance_threshold
    else:
        thresh = spatial.max_distance_threshold
    ax_dist.axhline(y=thresh, color=PALETTE["dark_purple"],
                    linestyle="--", linewidth=1, alpha=0.7)

    ax_dist.set_ylabel("Dist from\nmembrane (A)", fontsize=10)
    ax_dist.tick_params(labelsize=9)

    # ==================================================================
    # Track 3: SASA
    # ==================================================================
    valid_sasa = ~np.isnan(sasa_vals)
    if np.any(valid_sasa):
        ax_sasa.bar(
            positions[valid_sasa], sasa_vals[valid_sasa], width=1.0,
            color=PALETTE["teal"], alpha=0.6, linewidth=0,
        )
    ax_sasa.set_ylabel("SASA (A\u00B2)", fontsize=10)
    ax_sasa.set_ylim(bottom=0)
    ax_sasa.tick_params(labelsize=9)

    # ==================================================================
    # Track 4: Cyno conservation band
    # ==================================================================
    for i, color in enumerate(conservation_colors):
        ax_cons.axvspan(i + 0.5, i + 1.5, color=color, alpha=0.8)
    ax_cons.set_yticks([])
    ax_cons.set_ylabel("Cyno\nConserv.", fontsize=10, rotation=0, ha="right", va="center")

    if not conservation and not target.cyno_sequence:
        ax_cons.text(
            seq_len / 2.0, 0.5, "No cyno ortholog",
            ha="center", va="center", fontsize=10, fontstyle="italic",
            color="#888888", transform=ax_cons.get_xaxis_transform(),
        )

    # ==================================================================
    # Track 5: Human specificity band
    # ==================================================================
    for i, color in enumerate(specificity_colors):
        ax_spec.axvspan(i + 0.5, i + 1.5, color=color, alpha=0.8)
    ax_spec.set_yticks([])
    ax_spec.set_ylabel("Human\nSpecific.", fontsize=10, rotation=0, ha="right", va="center")

    # ==================================================================
    # Track 6: Qualifying patches
    # ==================================================================
    for i in range(seq_len):
        resnum = i + 1
        if resnum in patch_residues:
            ax_patch.axvspan(i + 0.5, i + 1.5, color=COLOR_EPITOPE_PATCH, alpha=0.9)

    ax_patch.set_yticks([])
    ax_patch.set_ylabel("Target\nEpitope", fontsize=10, rotation=0, ha="right", va="center")
    ax_patch.set_xlabel("Residue Position", fontsize=11)

    if scores:
        # Limit labels to top 5 patches to avoid clutter
        labeled = scores[:5]
        n = len(labeled)
        if n == 1:
            y_levels = [0.5]
        elif n == 2:
            y_levels = [0.75, 0.25]
        else:
            y_levels = [0.85 - i * 0.7 / (n - 1) for i in range(n)]
        for i, s in enumerate(labeled):
            y_pos = y_levels[i]
            center = np.mean(s.patch.residue_numbers)
            label = "P{} {:.0f}A\u00B2 | {:.2f}".format(
                s.rank, s.patch_area_a2, s.composite_score)
            ax_patch.text(
                center, y_pos, label,
                ha="center", va="center", fontsize=8.5, fontweight="bold",
                transform=ax_patch.get_xaxis_transform(),
            )
    else:
        ax_patch.text(
            seq_len / 2.0, 0.5, "No qualifying patches",
            ha="center", va="center", fontsize=10, fontstyle="italic",
            color="#888888", transform=ax_patch.get_xaxis_transform(),
        )

    # Hide x ticks for all but bottom track
    plt.setp(ax_dom.get_xticklabels(), visible=False)
    plt.setp(ax_dist.get_xticklabels(), visible=False)
    plt.setp(ax_sasa.get_xticklabels(), visible=False)
    plt.setp(ax_cons.get_xticklabels(), visible=False)
    plt.setp(ax_spec.get_xticklabels(), visible=False)


# ---------------------------------------------------------------------------
# Combined orientation figure (horizontal, both orientations side-by-side)
# ---------------------------------------------------------------------------

def plot_bispecific_combined(pair_result, output_path):
    """
    Combined figure showing both orientations side-by-side (horizontal).
    Only domain + epitope tracks are shown per target.

    Layout:
        Left half  = Orientation A (Target A distal, Target B proximal)
        Right half = Orientation B (Target B distal, Target A proximal)

    Each orientation shows two targets stacked:
        [Distal target:  domains row + epitope row]
        [Proximal target: domains row + epitope row]

    Args:
        pair_result: BispecificPairResult with orientation_ab and orientation_ba.
        output_path: Path to save PNG.
    """
    ori_ab = pair_result.orientation_ab
    ori_ba = pair_result.orientation_ba

    gene_a = pair_result.target_a.gene_name
    gene_b = pair_result.target_b.gene_name

    orientations = [
        (ori_ab, "{} (distal) \u00d7 {} (proximal)".format(gene_a, gene_b)),
        (ori_ba, "{} (distal) \u00d7 {} (proximal)".format(gene_b, gene_a)),
    ]

    # Figure width scales with longest sequence
    all_seq_lens = [
        ori_ab.distal_zone.target.sequence_length,
        ori_ab.proximal_zone.target.sequence_length,
        ori_ba.distal_zone.target.sequence_length,
        ori_ba.proximal_zone.target.sequence_length,
    ]
    max_seq = max(all_seq_lens)
    fig_width = max(18, max_seq * 0.025)
    fig_height = 13

    fig = plt.figure(figsize=(fig_width, fig_height))

    fig.suptitle(
        "BISPECIFIC: {} \u00d7 {}".format(gene_a, gene_b),
        fontsize=16, fontweight="bold", y=0.99,
    )
    fig.text(
        0.98, 0.99,
        "Pair score: {:.3f}{}".format(
            pair_result.final_pair_score,
            " [DUAL-VALID]" if pair_result.both_valid else "",
        ),
        ha="right", va="top", fontsize=13, fontweight="bold",
        color=PALETTE["dark_purple"],
    )

    # Outer grid: 1 row x 2 cols (one per orientation), with gap
    gs_outer = gridspec.GridSpec(
        1, 2, wspace=0.15,
        top=0.93, bottom=0.10, left=0.09, right=0.98,
    )

    for col_idx, (ori, col_title) in enumerate(orientations):
        dz = ori.distal_zone
        pz = ori.proximal_zone

        # Inner grid: 3 tracks per target (domains, distance, epitope) + gap
        gs_inner = gridspec.GridSpecFromSubplotSpec(
            7, 1, subplot_spec=gs_outer[col_idx],
            height_ratios=[0.7, 1.2, 0.5, 0.12, 0.7, 1.2, 0.5], hspace=0.08,
        )

        # --- Distal target ---
        ax_dom_d = fig.add_subplot(gs_inner[0])
        ax_dist_d = fig.add_subplot(gs_inner[1], sharex=ax_dom_d)
        ax_epi_d = fig.add_subplot(gs_inner[2], sharex=ax_dom_d)
        _draw_horizontal_domains(ax_dom_d, dz)
        _draw_horizontal_distance(ax_dist_d, dz)
        _draw_horizontal_epitope(ax_epi_d, dz)

        # Title above distal domains
        ax_dom_d.set_title(col_title, fontsize=12, fontweight="bold",
                           color=PALETTE["teal"], loc="left")

        # Zone + gene label on left
        zone_label = "DISTAL \u2265{:.0f}\u00c5".format(
            dz.spatial_filter.min_distance_threshold)
        ax_dom_d.set_ylabel(
            "{}\n{}".format(dz.target.gene_name, zone_label),
            fontsize=10, fontweight="bold", rotation=0, ha="right", va="center",
            color=PALETTE["teal"],
        )

        plt.setp(ax_dom_d.get_xticklabels(), visible=False)
        plt.setp(ax_dist_d.get_xticklabels(), visible=False)

        # --- Gap separator ---
        ax_gap = fig.add_subplot(gs_inner[3])
        ax_gap.set_visible(False)

        # --- Proximal target ---
        ax_dom_p = fig.add_subplot(gs_inner[4])
        ax_dist_p = fig.add_subplot(gs_inner[5], sharex=ax_dom_p)
        ax_epi_p = fig.add_subplot(gs_inner[6], sharex=ax_dom_p)
        _draw_horizontal_domains(ax_dom_p, pz)
        _draw_horizontal_distance(ax_dist_p, pz)
        _draw_horizontal_epitope(ax_epi_p, pz)

        zone_label_p = "PROXIMAL \u2264{:.0f}\u00c5".format(
            pz.spatial_filter.max_distance_threshold)
        ax_dom_p.set_ylabel(
            "{}\n{}".format(pz.target.gene_name, zone_label_p),
            fontsize=10, fontweight="bold", rotation=0, ha="right", va="center",
            color=PALETTE["blue"],
        )

        plt.setp(ax_dom_p.get_xticklabels(), visible=False)
        plt.setp(ax_dist_p.get_xticklabels(), visible=False)
        ax_epi_p.set_xlabel("Residue Position", fontsize=11)

        # Orientation score below
        score_label = "score: {:.3f}".format(ori.orientation_score)
        if not ori.is_valid:
            score_label += "  (no patches)"
        fig.text(
            gs_outer[col_idx].get_position(fig).x0
            + gs_outer[col_idx].get_position(fig).width / 2,
            0.04,
            score_label, ha="center", fontsize=12, fontweight="bold",
            color=PALETTE["dark_purple"] if ori.is_valid else "#999999",
        )

    # Legend
    legend_patches = [
        mpatches.Patch(color=PALETTE["teal"], label="InterPro Domain"),
        mpatches.Patch(color=PALETTE["blue"], label="UniProt Domain"),
        mpatches.Patch(color=COLOR_TRANSMEMBRANE, label="Transmembrane"),
        mpatches.Patch(color=COLOR_EPITOPE_PATCH, label="Epitope Patch"),
    ]
    fig.legend(
        handles=legend_patches, loc="lower center", ncol=4, fontsize=10,
        frameon=True, fancybox=True, bbox_to_anchor=(0.5, -0.01),
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Wrote combined bispecific: %s", output_path.name)


def _draw_horizontal_domains(ax, zone_result):
    """Draw the domain track for one target (horizontal, residues on x-axis)."""
    target = zone_result.target
    seq_len = target.sequence_length
    membrane = zone_result.membrane

    ax.set_xlim(0.5, seq_len + 0.5)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_ylabel("Domains", fontsize=9, rotation=0, ha="right", va="center")

    # Topology brackets
    _draw_topology_brackets(ax, target, membrane, seq_len)

    # Domain blocks
    domain_blocks = _collect_domain_blocks(target)
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
        ax.add_patch(rect)
        center = (start + end) / 2.0
        text_color = ("white" if fcolor in (COLOR_TRANSMEMBRANE, PALETTE["blue"])
                      else "black")
        fontsize = 7 if span < 80 else 9
        if span >= 20:
            txt = ax.text(
                center, block_y + block_h / 2, label,
                ha="center", va="center", fontsize=fontsize, fontweight="bold",
                color=text_color, zorder=4, clip_on=True,
            )
            right_clip = plt.Rectangle(
                (center - 5000, block_y), 5000 + (end - center + 0.5), block_h,
                transform=ax.transData, visible=False,
            )
            ax.add_patch(right_clip)
            txt.set_clip_path(right_clip)

    ax.tick_params(labelsize=9)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("#CCCCCC")


def _draw_horizontal_distance(ax, zone_result):
    """Draw the membrane-distance line plot for one target (horizontal)."""
    target = zone_result.target
    spatial = zone_result.spatial_filter
    membrane = zone_result.membrane
    seq_len = target.sequence_length
    zone = zone_result.zone

    positions = np.arange(1, seq_len + 1)
    distances = np.array([
        spatial.residue_distances.get(i, np.nan)
        for i in range(1, seq_len + 1)
    ])

    # Split into ECD vs non-ECD for separate fill colors
    ecd_dists = np.full(seq_len, np.nan)
    nonecd_dists = np.full(seq_len, np.nan)
    for i in range(seq_len):
        resnum = i + 1
        d = spatial.residue_distances.get(resnum)
        if d is None:
            continue
        topo = membrane.residue_topology.get(resnum)
        if topo == "extracellular":
            ecd_dists[i] = d
        else:
            nonecd_dists[i] = d

    ax.set_xlim(0.5, seq_len + 0.5)
    ax.fill_between(positions, ecd_dists, alpha=0.3, color=PALETTE["teal"])
    ax.fill_between(positions, nonecd_dists, alpha=0.08, color="#999999")
    ax.plot(positions, distances, color=PALETTE["teal"], linewidth=0.5, alpha=0.5)

    # Zone-specific threshold line
    if zone == "distal":
        thresh = spatial.min_distance_threshold
    else:
        thresh = spatial.max_distance_threshold
    ax.axhline(y=thresh, color=PALETTE["dark_purple"],
               linestyle="--", linewidth=1, alpha=0.7)

    ax.set_ylabel("Dist (\u00c5)", fontsize=9, rotation=0, ha="right", va="center")
    ax.tick_params(labelsize=8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("#CCCCCC")


def _draw_horizontal_epitope(ax, zone_result):
    """Draw the epitope-patch track for one target (horizontal)."""
    target = zone_result.target
    scores = zone_result.scores
    seq_len = target.sequence_length

    ax.set_xlim(0.5, seq_len + 0.5)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_ylabel("Epitope", fontsize=9, rotation=0, ha="right", va="center")

    # Collect patch residues
    patch_residues = set()
    for s in scores:
        for r in s.patch.residue_numbers:
            patch_residues.add(r)

    if patch_residues:
        for i in range(seq_len):
            resnum = i + 1
            if resnum in patch_residues:
                ax.axvspan(i + 0.5, i + 1.5, color=COLOR_EPITOPE_PATCH, alpha=0.9)

        # Limit labels to top 5 patches to avoid clutter
        labeled = scores[:5]
        n = len(labeled)
        if n == 1:
            y_levels = [0.5]
        elif n == 2:
            y_levels = [0.75, 0.25]
        else:
            y_levels = [0.85 - i * 0.7 / (n - 1) for i in range(n)]
        for i, s in enumerate(labeled):
            y_pos = y_levels[i]
            center = np.mean(s.patch.residue_numbers)
            label = "P{} {:.0f}A\u00B2 | {:.2f}".format(
                s.rank, s.patch_area_a2, s.composite_score)
            ax.text(
                center, y_pos, label,
                ha="center", va="center", fontsize=8.5, fontweight="bold",
                transform=ax.get_xaxis_transform(),
            )
    else:
        ax.text(
            seq_len / 2.0, 0.5, "No qualifying patches",
            ha="center", va="center", fontsize=10, fontstyle="italic",
            color="#888888", transform=ax.get_xaxis_transform(),
        )

    ax.tick_params(labelsize=9)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("#CCCCCC")
