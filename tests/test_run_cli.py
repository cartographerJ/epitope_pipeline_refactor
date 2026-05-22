"""Argparse / CLI tests for epitope_pipeline.run.main.

We monkeypatch run_pipeline to a recording stub so flags are exercised without
running the actual pipeline.
"""

import sys
from types import SimpleNamespace

import pytest

from epitope_pipeline import run as run_module
from epitope_pipeline import config


@pytest.fixture
def captured(monkeypatch):
    """Replace run_pipeline with a stub that records its kwargs."""
    calls = []

    def stub(identifiers, **kwargs):
        calls.append({"identifiers": list(identifiers), **kwargs})
        # Return a minimal result that the print loop in main() can consume
        return {"run_dir": "/tmp/fake", "targets": [], "metrics": {}}

    monkeypatch.setattr(run_module, "run_pipeline", stub)
    return calls


def _invoke(monkeypatch, argv):
    """Run main() with argv (excluding program name)."""
    monkeypatch.setattr(sys, "argv", ["epitope-pipeline", *argv])
    run_module.main()


def test_positional_targets(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "EGFR"])
    assert len(captured) == 1
    assert captured[0]["identifiers"] == ["ERBB2", "EGFR"]
    # Default mode is proximal
    assert captured[0]["max_distance_a"] == config.PROXIMAL_MAX_DISTANCE_A
    assert captured[0]["skip_cyno_gate"] is False
    assert captured[0]["cyno_mismatch_percent"] is None


def test_no_cyno_flag(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--no-cyno"])
    assert captured[0]["skip_cyno_gate"] is True


def test_cyno_mismatch_percent_flag(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--cyno-mismatch-percent", "25"])
    assert captured[0]["cyno_mismatch_percent"] == 25.0


def test_targets_file(captured, monkeypatch, tmp_path):
    targets_file = tmp_path / "list.txt"
    targets_file.write_text("ERBB2\n# a comment\n\nEGFR\nMSLN  # inline comment\n")
    _invoke(monkeypatch, ["--targets-file", str(targets_file)])
    assert captured[0]["identifiers"] == ["ERBB2", "EGFR", "MSLN"]


def test_targets_file_merges_with_positional_and_dedupes(captured, monkeypatch, tmp_path):
    targets_file = tmp_path / "list.txt"
    targets_file.write_text("EGFR\nMSLN\n")
    _invoke(monkeypatch, ["ERBB2", "EGFR", "--targets-file", str(targets_file)])
    # Order: positional first, then file; duplicates dropped on later occurrence
    assert captured[0]["identifiers"] == ["ERBB2", "EGFR", "MSLN"]


def test_targets_file_missing_errors(monkeypatch, tmp_path, capsys):
    bogus = tmp_path / "nope.txt"
    monkeypatch.setattr(sys, "argv", ["epitope-pipeline", "--targets-file", str(bogus)])
    with pytest.raises(SystemExit):
        run_module.main()


def test_no_targets_errors(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["epitope-pipeline"])
    with pytest.raises(SystemExit):
        run_module.main()
    err = capsys.readouterr().err
    assert "no targets" in err.lower()


def test_no_distance_filter_disables_proximal_default(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--no-distance-filter"])
    # When --no-distance-filter is set, the proximal default must NOT be applied
    assert captured[0]["max_distance_a"] is None
    assert captured[0]["no_distance_filter"] is True


def test_distal_mode_via_min_distance(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--min-distance", "80"])
    # Explicit --min-distance should suppress the proximal default
    assert captured[0]["min_distance_a"] == 80.0
    assert captured[0]["max_distance_a"] is None


def test_force_experimental(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--force-experimental"])
    assert captured[0]["force_experimental"] is True


def test_run_name_and_quiet(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--run-name", "my_run", "--quiet"])
    assert captured[0]["run_name"] == "my_run"
    assert captured[0]["verbose"] is False


def test_nonspecific_percent(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--nonspecific-percent", "70"])
    assert captured[0]["nonspecific_percent"] == 70.0


def test_split_alias_helper():
    from epitope_pipeline.run import _split_alias
    assert _split_alias("ERBB2") == ("ERBB2", None)
    assert _split_alias("ERBB2=HER2") == ("ERBB2", "HER2")
    assert _split_alias(" ERBB2 = HER2 ") == ("ERBB2", "HER2")
    assert _split_alias("ERBB2=") == ("ERBB2", None)  # empty alias → None
    assert _split_alias("P04626=Project_X") == ("P04626", "Project_X")


def test_alias_via_positional(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2=HER2", "EGFR=ErbB1", "MSLN"])
    c = captured[0]
    assert c["identifiers"] == ["ERBB2", "EGFR", "MSLN"]
    assert c["aliases"] == {"ERBB2": "HER2", "EGFR": "ErbB1"}


def test_alias_via_targets_file(captured, monkeypatch, tmp_path):
    f = tmp_path / "list.txt"
    f.write_text(
        "ERBB2=HER2\n"
        "# pick a codename for EGFR\n"
        "EGFR=Target_A\n"
        "MSLN\n"
    )
    _invoke(monkeypatch, ["--targets-file", str(f)])
    c = captured[0]
    assert c["identifiers"] == ["ERBB2", "EGFR", "MSLN"]
    assert c["aliases"] == {"ERBB2": "HER2", "EGFR": "Target_A"}


def test_no_aliases_passes_none(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "EGFR"])
    assert captured[0]["aliases"] is None


def test_combined_flags(captured, monkeypatch, tmp_path):
    targets_file = tmp_path / "list.txt"
    targets_file.write_text("ERBB2\nEGFR\n")
    _invoke(monkeypatch, [
        "--targets-file", str(targets_file),
        "--no-cyno",
        "--cyno-mismatch-percent", "30",
        "--run-name", "verify",
    ])
    c = captured[0]
    assert c["identifiers"] == ["ERBB2", "EGFR"]
    assert c["skip_cyno_gate"] is True
    assert c["cyno_mismatch_percent"] == 30.0
    assert c["run_name"] == "verify"


def test_cyno_mode_default_is_none(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2"])
    assert captured[0]["cyno_mode"] is None


def test_cyno_mode_whole_patch(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--cyno-mode", "whole-patch"])
    assert captured[0]["cyno_mode"] == "whole_patch"


def test_cyno_mode_local(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--cyno-mode", "local"])
    assert captured[0]["cyno_mode"] == "local"


def test_cyno_max_window_mismatches(captured, monkeypatch):
    _invoke(monkeypatch, ["ERBB2", "--cyno-max-window-mismatches", "3"])
    assert captured[0]["cyno_max_mismatches"] == 3
