# tests/test_viz_footer.py
"""Unit test for the distance-footer string helper used by plot_scoring_summary."""

from epitope_pipeline.viz.visualize import _distance_footer_str


def test_footer_whole_ecd():
    assert _distance_footer_str("whole_ecd", None) == (
        "whole ectodomain (no distance filter)"
    )
    # distance_value is ignored in whole-ECD mode
    assert _distance_footer_str("whole_ecd", 99.0) == (
        "whole ectodomain (no distance filter)"
    )


def test_footer_proximal():
    assert _distance_footer_str("proximal", 50.0) == "\u226450\u00c5 from membrane"


def test_footer_distal():
    assert _distance_footer_str("distal", 80.0) == "\u226580\u00c5 from membrane"


def test_footer_proximal_default_value():
    assert _distance_footer_str("proximal", None) == "\u226440\u00c5 from membrane"


def test_footer_distal_default_value():
    assert _distance_footer_str("distal", None) == "\u226580\u00c5 from membrane"
