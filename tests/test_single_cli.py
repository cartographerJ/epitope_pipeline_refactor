# tests/test_single_cli.py
"""Argparse / CLI tests for epitope_pipeline.single.main.

run_single is monkeypatched to a recording stub so flags are exercised
without running the actual pipeline.
"""

import sys

import pytest

from epitope_pipeline import single as single_module


@pytest.fixture
def captured(monkeypatch):
    calls = []

    def stub(identifiers, **kwargs):
        calls.append({"identifiers": list(identifiers), **kwargs})
        return {"run_dir": "/tmp/fake", "targets": [], "metrics": {}}

    monkeypatch.setattr(single_module, "run_single", stub)
    return calls


def _invoke(monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", ["epitope-single", *argv])
    single_module.main()


def test_positional_targets(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "EGFR"])
    assert captured[0]["identifiers"] == ["ERBB2", "EGFR"]
    assert captured[0]["skip_cyno_gate"] is False
    assert captured[0]["cyno_mismatch_percent"] is None


def test_aliases_and_flags(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2=HER2", "--no-cyno", "--nonspecific-percent", "70"])
    c = captured[0]
    assert c["identifiers"] == ["ERBB2"]
    assert c["aliases"] == {"ERBB2": "HER2"}
    assert c["skip_cyno_gate"] is True
    assert c["nonspecific_percent"] == 70.0


def test_targets_file_merges_and_dedupes(captured, monkeypatch, tmp_path):
    f = tmp_path / "list.txt"
    f.write_text("EGFR\nMSLN\n")
    _invoke(monkeypatch, ["ERBB2", "EGFR", "--targets-file", str(f)])
    assert captured[0]["identifiers"] == ["ERBB2", "EGFR", "MSLN"]


def test_no_targets_errors(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["epitope-single"])
    with pytest.raises(SystemExit):
        single_module.main()
    assert "no targets" in capsys.readouterr().err.lower()


@pytest.mark.parametrize("flag", ["--min-distance", "--max-distance", "--no-distance-filter"])
def test_distance_flags_rejected(monkeypatch, capsys, flag):
    """The single-protein CLI must NOT accept distal/proximal distance flags."""
    argv = ["epitope-single", "ERBB2", flag]
    if flag != "--no-distance-filter":
        argv.append("50")
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SystemExit):
        single_module.main()
    assert "unrecognized arguments" in capsys.readouterr().err.lower()


def test_single_cyno_mode_whole_patch(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--cyno-mode", "whole-patch"])
    assert captured[0]["cyno_mode"] == "whole_patch"


def test_single_cyno_mode_default_none(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2"])
    assert captured[0]["cyno_mode"] is None


def test_single_cyno_max_window_mismatches(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--cyno-max-window-mismatches", "4"])
    assert captured[0]["cyno_max_mismatches"] == 4
