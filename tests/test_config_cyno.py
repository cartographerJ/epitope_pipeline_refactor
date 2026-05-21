# tests/test_config_cyno.py
"""Cyno-mode config defaults."""

from epitope_pipeline import config


def test_cyno_mode_defaults_to_local():
    assert config.CYNO_MODE == "local"


def test_local_threshold_default():
    assert config.MAX_CYNO_MISMATCHES_PER_600A2 == 2


def test_whole_patch_percent_default_kept():
    assert config.MAX_CYNO_MISMATCH_PERCENT == 15.0
