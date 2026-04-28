"""
Epitope Pipeline — Local Streamlit App

Run: streamlit run epitope_pipeline/app.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")

import streamlit as st
import re
import logging
import io
from pathlib import Path

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="Epitope Pipeline",
    page_icon="\U0001f9ec",
    layout="wide",
)

COLORS = {
    "gray": "#D3D3D3",
    "mint": "#E0F3DB",
    "green": "#6BC291",
    "teal": "#18B5CB",
    "blue": "#2E95D2",
    "purple": "#28154C",
}

# ---------------------------------------------------------------------------
# Session state
# ---------------------------------------------------------------------------
if "results" not in st.session_state:
    st.session_state.results = None
    st.session_state.mode = "Single-target"

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
st.sidebar.markdown(
    "<h2 style='color: {};'>Epitope Pipeline</h2>".format(COLORS["teal"]),
    unsafe_allow_html=True,
)

mode = st.sidebar.radio("Mode", ["Single-target", "Bispecific"], horizontal=True)
st.session_state.mode = mode

st.sidebar.markdown("---")

if mode == "Single-target":
    target_text = st.sidebar.text_area(
        "Target genes",
        placeholder="ERBB2\nEGFR\nNECTIN4",
        help="Gene names or UniProt IDs, one per line or comma-separated",
        height=120,
    )
    no_dist_filter = st.sidebar.checkbox("No distance filter", value=False,
        help="Include all extracellular residues regardless of distance from membrane")
    if not no_dist_filter:
        max_dist = st.sidebar.slider(
            "Max distance from membrane (\u00c5)",
            min_value=20, max_value=120, value=50, step=5,
            help="Find epitopes within this distance of the membrane (proximal mode)",
        )
    else:
        max_dist = None
else:
    target_text = st.sidebar.text_area(
        "Target pairs",
        placeholder="ERBB2:NECTIN4\nROR1:EGFR",
        help="Colon-separated pairs, one per line",
        height=120,
    )
    distal_min = st.sidebar.slider(
        "Distal min distance (\u00c5)",
        min_value=40, max_value=120, value=60, step=5,
    )
    proximal_max = st.sidebar.slider(
        "Proximal max distance (\u00c5)",
        min_value=20, max_value=80, value=50, step=5,
    )

st.sidebar.markdown("---")

cyno_mismatch_pct = st.sidebar.slider(
    "Max cyno mismatch %",
    min_value=5, max_value=50, value=15, step=5,
    help="Max % cyno mismatches to accept a patch (e.g., 15% = 3 mismatches in a 20-residue patch)",
)

nonspecific_pct = st.sidebar.slider(
    "Max % shared with off-targets",
    min_value=50, max_value=95, value=85, step=5,
    help="Max % of patch residues that can be shared with off-target human proteins. Higher = more permissive. E.g., 85% = patch can share up to 85% of residues with a paralog (15% must be unique).",
)

with st.sidebar.expander("Advanced"):
    force_experimental = st.checkbox("Force experimental PDB", value=False)

st.sidebar.markdown("---")

run_clicked = st.sidebar.button(
    "Run Pipeline",
    type="primary",
    use_container_width=True,
)

# ---------------------------------------------------------------------------
# Run pipeline (synchronous — blocks until done, shows live log)
# ---------------------------------------------------------------------------
if run_clicked:
    # Parse targets
    raw = target_text.strip()
    if not raw:
        st.sidebar.error("Enter at least one target")
        st.stop()

    if mode == "Single-target":
        identifiers = [t.strip().upper() for t in re.split(r"[,\n\s]+", raw) if t.strip()]
        if not identifiers:
            st.sidebar.error("No valid targets found")
            st.stop()
    else:
        pairs = []
        for line in raw.splitlines():
            line = line.strip()
            if not line:
                continue
            if ":" not in line:
                st.sidebar.error("Pairs must be colon-separated (e.g. ERBB2:NECTIN4)")
                st.stop()
            a, b = line.split(":", 1)
            pairs.append((a.strip().upper(), b.strip().upper()))
        if not pairs:
            st.sidebar.error("No valid pairs found")
            st.stop()

    # Capture log output
    log_stream = io.StringIO()
    logger = logging.getLogger("epitope_pipeline")
    handler = logging.StreamHandler(log_stream)
    handler.setFormatter(logging.Formatter("%(message)s"))
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)

    try:
        with st.status("Running pipeline...", expanded=True) as status:
            if mode == "Single-target":
                from epitope_pipeline.run import run_pipeline
                results = run_pipeline(
                    identifiers,
                    max_distance_a=float(max_dist) if max_dist is not None else None,
                    no_distance_filter=no_dist_filter,
                    cyno_mismatch_percent=float(cyno_mismatch_pct),
                    nonspecific_percent=float(nonspecific_pct),
                    force_experimental=force_experimental,
                )
            else:
                from epitope_pipeline.bispecific import run_bispecific
                results = run_bispecific(
                    pairs,
                    distal_min_distance_a=float(distal_min),
                    proximal_max_distance_a=float(proximal_max),
                    cyno_mismatch_percent=float(cyno_mismatch_pct),
                    nonspecific_percent=float(nonspecific_pct),
                    force_experimental=force_experimental,
                )

            st.session_state.results = results
            status.update(label="Pipeline complete", state="complete")

    except Exception as e:
        st.error("Pipeline failed: {}".format(e))
        import traceback
        st.code(traceback.format_exc())
        st.session_state.results = None

    finally:
        logger.removeHandler(handler)

    # Show log
    log_text = log_stream.getvalue()
    if log_text:
        with st.expander("Pipeline log", expanded=False):
            st.code(log_text, language=None)

# ---------------------------------------------------------------------------
# Results display
# ---------------------------------------------------------------------------
results = st.session_state.results
if results is not None:
    run_dir = Path(str(results["run_dir"]))

    st.markdown("### Results \u2014 `{}`".format(run_dir.name))

    # Open in Finder button
    col_path, col_btn = st.columns([4, 1])
    col_path.code(str(run_dir))
    if col_btn.button("Open in Finder"):
        import subprocess
        subprocess.run(["open", str(run_dir)])

    st.markdown("---")

    if st.session_state.mode == "Single-target":
        scores = results.get("scores", {})
        metrics = results.get("metrics", {})
        rows = []
        for target in results.get("targets", []):
            uid = target.uniprot_id
            n_patches = len(scores.get(uid, []))
            best = max((s.composite_score for s in scores.get(uid, [])), default=0.0)
            rows.append({
                "Gene": target.gene_name,
                "UniProt": uid,
                "Patches": n_patches,
                "Best Score": round(best, 3),
                "Total Area (\u00c5\u00b2)": round(
                    sum(s.patch_area_a2 for s in scores.get(uid, [])), 0
                ),
            })
        if rows:
            st.dataframe(rows, use_container_width=True, hide_index=True)

    else:
        rows = []
        for pr in results.get("pair_results", []):
            rows.append({
                "Target A": pr.target_a.gene_name,
                "Target B": pr.target_b.gene_name,
                "Pair Score": round(pr.final_pair_score, 3),
                "Dual-Valid": "Yes" if pr.both_valid else "No",
                "Best Orientation": pr.best_orientation.label,
            })
        if rows:
            st.dataframe(rows, use_container_width=True, hide_index=True)

    # Figures
    fig_dir = run_dir / "Figures"
    if fig_dir.exists():
        pngs = sorted(fig_dir.glob("*.png"))
        zone_dir = fig_dir / "zone_details"
        if zone_dir.exists():
            pngs.extend(sorted(zone_dir.glob("*.png")))

        if pngs:
            st.markdown("### Figures")
            for png in pngs:
                st.image(str(png), caption=png.stem, use_container_width=True)
