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
    # Font sizes scaled for slide legibility at 30-40% of slide width
    fig, axes = plt.subplots(
        6, 1,
        figsize=(28.6, 16),
        gridspec_kw={"height_ratios": [0.7, 1.4, 1.4, 0.4, 0.4, 0.5]},
        sharex=True,
    )
    title = "{} ({}) \u2014 Epitope Analysis".format(target.gene_name, target.uniprot_id)
    if title_suffix:
        title += " \u2014 {}".format(title_suffix)
    fig.suptitle(title, fontsize=36, fontweight="bold", y=0.94)

    # Y-axis labels: place at a fixed x in figure coords so they left-align
    _YLABEL_X = 0.06  # fixed x position (figure fraction)
    _YLABEL_KW = dict(fontsize=24, ha="left", va="center", fontweight="bold")

    def _ylabel(ax, text):
        """Place y-axis label at fixed x so all labels left-align."""
        mid_y = (ax.get_position().y0 + ax.get_position().y1) / 2.0
        fig.text(_YLABEL_X, mid_y, text, transform=fig.transFigure, **_YLABEL_KW)

    # ==================================================================
    # Track 1: Domains (compact blocks with topology brackets above)
    # ==================================================================
    ax_dom = axes[0]
    ax_dom.set_xlim(0.5, seq_len + 0.5)
    ax_dom.set_ylim(0, 1)
    ax_dom.set_yticks([])

    # --- Topology brackets above the domain blocks ---
    _draw_topology_brackets(ax_dom, target, membrane, seq_len)

    # --- Collect and render domain blocks ---
    domain_blocks = _collect_domain_blocks(target)

    # Draw largest first so smaller features appear on top
    domain_blocks.sort(key=lambda x: -(x[1] - x[0]))

    block_y = 0.05
    block_h = 0.55

    # Collect TM ranges for visible-gap calculation
    _tm_ranges = [(s, e) for s, e, l, c in domain_blocks
                  if l == "TM" or l == "GPI"]

    for start, end, label, fcolor in domain_blocks:
        span = end - start + 1
        is_tm = label in ("TM", "SP", "GPI")
        rect = plt.Rectangle(
            (start - 0.5, block_y), span, block_h,
            facecolor=fcolor, edgecolor="#333333", linewidth=0.7,
            alpha=0.85, zorder=4 if is_tm else 3,
        )
        ax_dom.add_patch(rect)
        text_color = "white" if fcolor == COLOR_TRANSMEMBRANE else "black"

        # For blocks that overlap TM segments (e.g. InterPro rhodopsin),
        # use the largest visible gap between TMs for text sizing
        effective_span = span
        if not is_tm:
            overlapping_tms = [(max(s, start), min(e, end))
                               for s, e in _tm_ranges
                               if s <= end and e >= start]
            if overlapping_tms:
                boundaries = sorted(set([start] + [e + 1 for _, e in overlapping_tms]
                                        + [s for s, _ in overlapping_tms] + [end + 1]))
                gaps = []
                for i in range(len(boundaries) - 1):
                    gap_s, gap_e = boundaries[i], boundaries[i + 1]
                    in_tm = any(ts <= gap_s and gap_e <= te + 1
                                for ts, te in overlapping_tms)
                    if not in_tm:
                        gaps.append(gap_e - gap_s)
                if gaps:
                    effective_span = max(gaps)

        center = (start + end) / 2.0
        fig_width_pts = 28.6 * 72  # 2059 points
        block_pts = (effective_span / seq_len) * fig_width_pts
        char_width = 30.0  # conservative for slide-legible bold fontsize
        max_chars = max(0, int(block_pts / char_width) - 1)

        if max_chars < 2:
            continue  # too small for any label
        fontsize = 16 if span < 80 else 20
        display_label = label[:max_chars] + ".." if len(label) > max_chars else label

        # Place text in the center of the largest visible gap
        if fcolor not in (COLOR_TRANSMEMBRANE,) and label not in ("TM", "SP", "GPI") and _tm_ranges:
            overlapping_tms = [(max(s, start), min(e, end))
                               for s, e in _tm_ranges
                               if s <= end and e >= start]
            if overlapping_tms:
                boundaries = sorted(set([start] + [e + 1 for _, e in overlapping_tms]
                                        + [s for s, _ in overlapping_tms] + [end + 1]))
                best_gap_center = center
                best_gap_size = 0
                for i in range(len(boundaries) - 1):
                    gap_s, gap_e = boundaries[i], boundaries[i + 1]
                    in_tm = any(ts <= gap_s and gap_e <= te + 1
                                for ts, te in overlapping_tms)
                    if not in_tm and (gap_e - gap_s) > best_gap_size:
                        best_gap_size = gap_e - gap_s
                        best_gap_center = (gap_s + gap_e) / 2.0
                center = best_gap_center

        txt = ax_dom.text(
            center, block_y + block_h / 2, display_label,
            ha="center", va="center", fontsize=fontsize, fontweight="bold",
            color=text_color, zorder=4, clip_on=True,
        )
        txt.set_clip_path(rect)

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
        linestyle="--", linewidth=2, alpha=0.7,
    )
    # y-label set after tight_layout
    ax_dist.tick_params(labelsize=20)

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

    # y-label set after tight_layout
    ax_sasa.set_ylim(bottom=0)
    ax_sasa.tick_params(labelsize=20)

    # ==================================================================
    # Track 4: Cyno conservation (per-residue band)
    # ==================================================================
    ax_cons = axes[3]
    for i, color in enumerate(conservation_colors):
        ax_cons.axvspan(i + 0.5, i + 1.5, color=color, alpha=0.8)
    ax_cons.set_yticks([])
    # y-label set after tight_layout

    # Label when no cyno ortholog exists
    if not conservation_result and not target.cyno_sequence:
        ax_cons.text(
            seq_len / 2.0, 0.5, "No cyno ortholog found",
            ha="center", va="center", fontsize=22, fontstyle="italic",
            color="#888888", transform=ax_cons.get_xaxis_transform(),
        )

    # ==================================================================
    # Track 5: Human specificity (per-residue band, same style as cyno)
    # ==================================================================
    ax_spec = axes[4]
    for i, color in enumerate(specificity_colors):
        ax_spec.axvspan(i + 0.5, i + 1.5, color=color, alpha=0.8)
    ax_spec.set_yticks([])
    # y-label set after tight_layout

    # ==================================================================
    # Track 6: Target epitope (final qualifying patches)
    # ==================================================================
    ax_patch = axes[5]

    for i in range(seq_len):
        resnum = i + 1
        if resnum in patch_residues:
            ax_patch.axvspan(i + 0.5, i + 1.5, color=COLOR_EPITOPE_PATCH, alpha=0.9)

    ax_patch.set_yticks([])
    # y-label set after tight_layout
    ax_patch.set_xlabel("Residue Position", fontsize=26)
    ax_patch.tick_params(labelsize=20)

    # Patch labels (score + rank) — stagger y positions to avoid overlap
    if scores:
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
            label = "P{} {:.0f}\u00c5\u00B2 | {:.2f}".format(
                s.rank, s.patch_area_a2, s.composite_score)
            ax_patch.text(
                center, y_pos, label,
                ha="center", va="center", fontsize=20, fontweight="bold",
                transform=ax_patch.get_xaxis_transform(),
            )
    else:
        ax_patch.text(
            seq_len / 2.0, 0.5, "No qualifying epitope patches",
            ha="center", va="center", fontsize=22, fontstyle="italic",
            color="#888888", transform=ax_patch.get_xaxis_transform(),
        )

    # ==================================================================
    # Legend
    # ==================================================================
    legend_patches = [
        mpatches.Patch(color=PALETTE["teal"], label="Domain"),
        mpatches.Patch(color=COLOR_TRANSMEMBRANE, label="Transmembrane"),
        mpatches.Patch(color=PALETTE["mint"], label="Conserved / Specific"),
        mpatches.Patch(color=COLOR_MISMATCH, label="Mismatch / Non-specific"),
        mpatches.Patch(color=COLOR_EPITOPE_PATCH, label="Target Epitope"),
    ]
    fig.legend(
        handles=legend_patches, loc="lower center", ncol=3, fontsize=18,
        frameon=True, fancybox=True,
    )

    # Filter parameters subtitle (under title)
    from epitope_pipeline import config
    if getattr(spatial_filter, "max_distance_threshold", None) is not None:
        dist_str = "\u2264{:.0f}\u00c5 from membrane".format(spatial_filter.max_distance_threshold)
    else:
        dist_str = "\u2265{:.0f}\u00c5 from membrane".format(spatial_filter.min_distance_threshold)
    cyno_pct = 100.0 - config.MAX_CYNO_MISMATCH_PERCENT
    spec_pct = 100.0 - config.MAX_NONSPECIFIC_PERCENT
    subtitle = "Filters: {} | SASA >{:.0f}% | \u2265{:.0f}% cyno conserved | \u2265{:.0f}% specific (scaled)".format(
        dist_str, config.SURFACE_EXPOSURE_THRESHOLD * 100, cyno_pct, spec_pct)
    fig.text(0.5, 0.895, subtitle, ha="center", fontsize=16, color="#888888", fontstyle="italic")

    plt.tight_layout(rect=[0.14, 0.06, 1, 0.93])

    # Place y-axis labels AFTER tight_layout so axis positions are finalized
    _ylabel(ax_dom, "Domains")
    _ylabel(ax_dist, "Membrane\ndist. (\u00c5)")
    _ylabel(ax_sasa, "SASA (\u00c5\u00b2)")
    _ylabel(ax_cons, "Cyno\nConserv.")
    _ylabel(ax_spec, "Human\nSpecific.")
    _ylabel(ax_patch, "Target\nEpitope")

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
        span = end - start + 1
        # Horizontal line
        ax.plot([start, end], [bracket_y, bracket_y],
                color=color, linewidth=1.5, clip_on=False, zorder=5)
        # Downward ticks at ends
        ax.plot([start, start], [bracket_y, bracket_y - tick_h],
                color=color, linewidth=1.5, clip_on=False, zorder=5)
        ax.plot([end, end], [bracket_y, bracket_y - tick_h],
                color=color, linewidth=1.5, clip_on=False, zorder=5)
        # Centered label above — skip if span too narrow for text
        fig_width_pts = 20 * 72
        block_pts = (span / seq_len) * fig_width_pts
        max_chars = max(0, int(block_pts / 14.0) - 1)
        if max_chars >= 3:
            display_label = label[:max_chars] + ".." if len(label) > max_chars else label
            center = (start + end) / 2.0
            ax.text(center, label_y, display_label,
                    ha="center", va="bottom", fontsize=18, color=color,
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
        "Domain": PALETTE["teal"],
        "Region": "#BBBBBB",
        "InterPro": PALETTE["teal"],
        "Chain": PALETTE["teal"],
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


def plot_scoring_summary(all_scores, target_metrics, targets, output_path,
                         distance_label=None, distance_value=None, distance_mode=None):
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
    # Collect data, sort by target score descending
    entries = []
    for t in targets:
        uid = t.uniprot_id
        metric = target_metrics.get(uid, {})
        entries.append({
            "gene": t.gene_name,
            "target_score": metric.get("target_score", 0.0),
            "area_component": metric.get("area_component", 0.0),
            "quality_component": metric.get("quality_component", 0.0),
            "total_area": metric.get("total_epitope_area_a2", 0.0),
            "n_patches": metric.get("n_patches", 0),
        })
    entries.sort(key=lambda e: e["target_score"], reverse=True)

    gene_names = [e["gene"] for e in entries]
    target_scores = [e["target_score"] for e in entries]
    area_components = [e["area_component"] for e in entries]
    quality_components = [e["quality_component"] for e in entries]
    total_areas = [e["total_area"] for e in entries]
    n_patches_list = [e["n_patches"] for e in entries]
    n_targets = len(entries)

    fig, axes = plt.subplots(1, 2, figsize=(max(10, 2.5 * n_targets), 5))
    fig.suptitle("Epitope Druggability Summary", fontsize=16, fontweight="bold")

    x = np.arange(n_targets)
    bar_width = 0.6

    # Panel 1: Target druggability score
    ax1 = axes[0]
    bars1 = ax1.bar(x, target_scores, bar_width, color=PALETTE["dark_purple"], alpha=0.85)
    ax1.set_xticks(x)
    ax1.set_xticklabels(gene_names, fontsize=12)
    ax1.set_ylabel("Target Score", fontsize=12)
    ax1.set_title("Druggability Score", fontsize=13)
    ax1.set_ylim(0, 1.05)
    for bar, val in zip(bars1, target_scores):
        if val > 0:
            ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                     "{:.2f}".format(val), ha="center", va="bottom", fontsize=11)

    # Panel 2: Total epitope area with patch count annotations
    ax2 = axes[1]
    bars2 = ax2.bar(x, total_areas, bar_width, color=PALETTE["teal"], alpha=0.8)
    ax2.set_xticks(x)
    ax2.set_xticklabels(gene_names, fontsize=12)
    ax2.set_ylabel("Total Epitope Area (\u00c5\u00b2)", fontsize=12)
    ax2.set_title("Epitope Surface", fontsize=13)
    for bar, val, n in zip(bars2, total_areas, n_patches_list):
        if val > 0:
            label = "{:.0f} ({})".format(val, n)
            ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                     label, ha="center", va="bottom", fontsize=10)

    # Footer with filter parameters
    from epitope_pipeline import config
    if distance_mode == "proximal":
        dist_str = "\u2264{:.0f}\u00c5 from membrane".format(distance_value or 40)
    else:
        dist_str = "\u2265{:.0f}\u00c5 from membrane".format(distance_value or 80)
    cyno_conserved = 100.0 - config.MAX_CYNO_MISMATCH_PERCENT
    specific_pct = 100.0 - config.MAX_NONSPECIFIC_PERCENT
    footer = "Filtered: {}, SASA >{:.0f}%, \u2265{:.0f}% cyno conserved (scaled), \u2265{:.0f}% specific (scaled)".format(
        dist_str, config.SURFACE_EXPOSURE_THRESHOLD * 100,
        cyno_conserved, specific_pct,
    )
    fig.text(0.5, -0.02, footer, ha="center", fontsize=9, color="#888888", fontstyle="italic")

    plt.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Wrote scoring summary: %s", output_path.name)


# ======================================================================
# BLAST off-target dot plot
# ======================================================================

def plot_blast_offtargets(target, specificity_result, output_path,
                          max_proteins=25, identity_threshold=0.40):
    """
    Dot plot of BLAST off-target HSPs along the target sequence.

    Y-axis = off-target protein names (sorted by max identity, descending).
    X-axis = position along target sequence.
    Each dot = one HSP, colored by identity (teal gradient).
    Top hits (>=70% identity) highlighted with a red box.

    Args:
        target: TargetInfo.
        specificity_result: SpecificityResult with full_blast_hits.
        output_path: Path to save the figure.
        max_proteins: Max number of off-target proteins to show.
        identity_threshold: Min HSP identity to include (default 0.40).
    """
    if not specificity_result or not specificity_result.full_blast_hits:
        return

    hits = specificity_result.full_blast_hits
    seq_len = target.sequence_length

    # Filter by identity threshold and group by protein
    protein_hsps = {}  # {gene_name: [hsp_dicts]}
    for h in hits:
        if h["identity"] < identity_threshold:
            continue
        # Extract conventional gene name from BLAST title
        # Title format: "sp|Q9BYE2|TMPSD_HUMAN ... GN=TMPRSS13 ..."
        # Prefer GN= field (conventional gene name), fall back to UniProt short name
        title = h.get("title", "")
        import re as _re
        gn_match = _re.search(r'GN=(\S+)', title)
        if gn_match:
            gene = gn_match.group(1)
        elif "|" in title:
            pipe_parts = title.split("|")
            gene = pipe_parts[2].split("_")[0] if len(pipe_parts) >= 3 else pipe_parts[-1].split("_")[0]
        else:
            parts = title.split()
            gene = parts[0].split("_")[0] if parts and "_" in parts[0] else h.get("accession", "?")
        if gene not in protein_hsps:
            protein_hsps[gene] = []
        protein_hsps[gene].append(h)

    if not protein_hsps:
        return

    # Sort proteins by max identity (descending), take top N
    protein_max_id = {g: max(h["identity"] for h in hsps)
                      for g, hsps in protein_hsps.items()}
    sorted_proteins = sorted(protein_max_id, key=protein_max_id.get, reverse=True)
    sorted_proteins = sorted_proteins[:max_proteins]

    # Build plot data
    y_labels = sorted_proteins
    n_proteins = len(y_labels)
    gene_to_y = {g: i for i, g in enumerate(reversed(y_labels))}

    fig_height = max(4, 0.35 * n_proteins + 1.5)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    # Identify non-specific threshold (>=70% identity)
    nonspecific_genes = [g for g in sorted_proteins if protein_max_id[g] >= 0.70]

    # Plot dots
    xs, ys, colors = [], [], []
    for gene in sorted_proteins:
        y = gene_to_y[gene]
        for h in protein_hsps[gene]:
            # Plot dot at midpoint of HSP range
            mid = (h["query_start"] + h["query_end"]) / 2.0
            xs.append(mid)
            ys.append(y)
            colors.append(PALETTE["teal"])

    ax.scatter(xs, ys, c=colors, s=30, alpha=0.7, edgecolors="none", zorder=3)

    # Red box around non-specific hits (>=70% identity)
    if nonspecific_genes:
        n_nonspec = len(nonspecific_genes)
        box_y_min = gene_to_y[nonspecific_genes[-1]] - 0.5
        box_y_max = gene_to_y[nonspecific_genes[0]] + 0.5
        from matplotlib.patches import Rectangle
        rect = Rectangle(
            (0, box_y_min), seq_len, box_y_max - box_y_min,
            linewidth=2, edgecolor="red", facecolor="none", zorder=2,
        )
        ax.add_patch(rect)
        ax.text(
            seq_len, box_y_max - 0.1,
            " \u226570% identity ({})".format(n_nonspec),
            fontsize=8, color="red", va="top", ha="left",
        )

    # Styling
    ax.set_xlim(0, seq_len * 1.15)  # extra space for label
    ax.set_ylim(-0.5, n_proteins - 0.5)
    ax.set_yticks(range(n_proteins))
    ax.set_yticklabels(reversed(y_labels), fontsize=9)
    ax.set_xlabel("Position", fontsize=12)
    ax.set_title("{} \u2014 BLAST Off-Targets".format(
        target.gene_name), fontsize=13, fontweight="bold")
    ax.grid(axis="y", alpha=0.3, linewidth=0.5)
    ax.tick_params(axis="x", labelsize=10)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    plt.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("  Wrote BLAST off-target plot: %s", output_path.name)
