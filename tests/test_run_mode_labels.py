# tests/test_run_mode_labels.py
"""Unit tests for the pure mode/label helpers in epitope_pipeline.run."""

from epitope_pipeline.run import _mode_label
from epitope_pipeline import config


def test_mode_label_whole_ecd_takes_priority():
    # no_distance_filter wins even if a max_distance is somehow present
    assert _mode_label(no_distance_filter=True, max_distance_a=50.0) == (
        "whole ectodomain (no distance filter)"
    )
    assert _mode_label(no_distance_filter=True, max_distance_a=None) == (
        "whole ectodomain (no distance filter)"
    )


def test_mode_label_proximal():
    assert _mode_label(no_distance_filter=False, max_distance_a=50.0) == (
        "proximal (max_distance_a=50.0)"
    )


def test_mode_label_distal_default():
    assert _mode_label(no_distance_filter=False, max_distance_a=None) == (
        "distal (min_distance_a={})".format(config.ECTODOMAIN_MIN_DISTANCE_A)
    )
