# tests/test_run_mode_labels.py
"""Unit tests for the pure mode/label helpers in epitope_pipeline.run."""

from epitope_pipeline.run import _mode_label, _summary_distance_args
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


def test_summary_distance_args_whole_ecd():
    args = _summary_distance_args(no_distance_filter=True, max_distance_a=None)
    assert args == {
        "distance_label": "whole ECD",
        "distance_value": None,
        "distance_mode": "whole_ecd",
    }
    # no_distance_filter wins even if a max_distance is somehow present
    assert _summary_distance_args(
        no_distance_filter=True, max_distance_a=50.0
    )["distance_mode"] == "whole_ecd"


def test_summary_distance_args_proximal():
    args = _summary_distance_args(no_distance_filter=False, max_distance_a=50.0)
    assert args["distance_mode"] == "proximal"
    assert args["distance_value"] == 50.0
    assert args["distance_label"] == "\u2264{:.0f}\u00c5".format(50.0)


def test_summary_distance_args_distal():
    args = _summary_distance_args(no_distance_filter=False, max_distance_a=None)
    assert args["distance_mode"] == "distal"
    assert args["distance_value"] == config.ECTODOMAIN_MIN_DISTANCE_A
    assert args["distance_label"] == "\u2265{:.0f}\u00c5".format(config.ECTODOMAIN_MIN_DISTANCE_A)
