# tests/test_bispecific_cli.py
"""Argparse / CLI tests for epitope_pipeline.bispecific.main.

run_bispecific is monkeypatched to a recording stub so flags are exercised
without running the actual pipeline.
"""

import sys

import pytest

from epitope_pipeline import bispecific as bispecific_module


@pytest.fixture
def captured(monkeypatch):
    calls = []

    def stub(pairs, **kwargs):
        calls.append({"pairs": list(pairs), **kwargs})
        return {"run_dir": "/tmp/fake", "pair_results": [], "zone_results": {}}

    monkeypatch.setattr(bispecific_module, "run_bispecific", stub)
    return calls


def _invoke(monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", ["epitope-bispecific", *argv])
    bispecific_module.main()


def test_pair_parsing(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2:NECTIN4"])
    assert captured[0]["pairs"] == [("ERBB2", "NECTIN4")]


def test_cyno_mode_whole_patch(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2:NECTIN4", "--cyno-mode", "whole-patch"])
    assert captured[0]["cyno_mode"] == "whole_patch"


def test_cyno_mode_default_none(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2:NECTIN4"])
    assert captured[0]["cyno_mode"] is None


def test_cyno_max_window_mismatches(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2:NECTIN4", "--cyno-max-window-mismatches", "1"])
    assert captured[0]["cyno_max_mismatches"] == 1
