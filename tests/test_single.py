# tests/test_single.py
"""Tests for the single-protein whole-ECD wrapper run_single."""

from epitope_pipeline import single as single_module


def test_run_single_forces_no_distance_filter(monkeypatch):
    """run_single must delegate to run_pipeline with no_distance_filter=True
    and must NOT pass any min/max distance."""
    captured = {}

    def stub(identifiers, **kwargs):
        captured["identifiers"] = list(identifiers)
        captured.update(kwargs)
        return {"run_dir": "/tmp/fake", "targets": [], "metrics": {}}

    monkeypatch.setattr(single_module, "run_pipeline", stub)

    single_module.run_single(["ERBB2"], cyno_mismatch_percent=20.0,
                             skip_cyno_gate=True)

    assert captured["identifiers"] == ["ERBB2"]
    assert captured["no_distance_filter"] is True
    assert captured["cyno_mismatch_percent"] == 20.0
    assert captured["skip_cyno_gate"] is True
    # The distal/proximal knobs must not be forwarded by the single-protein API
    assert "min_distance_a" not in captured
    assert "max_distance_a" not in captured


def test_run_single_default_run_name_marks_single(monkeypatch):
    """When no run_name is given, the auto slug identifies this as a single run."""
    captured = {}

    def stub(identifiers, **kwargs):
        captured.update(kwargs)
        return {"run_dir": "/tmp/fake", "targets": [], "metrics": {}}

    monkeypatch.setattr(single_module, "run_pipeline", stub)
    single_module.run_single(["ERBB2"])
    assert "_single_erbb2" in captured["run_name"]


def test_run_single_explicit_run_name_passes_through(monkeypatch):
    captured = {}

    def stub(identifiers, **kwargs):
        captured.update(kwargs)
        return {"run_dir": "/tmp/fake", "targets": [], "metrics": {}}

    monkeypatch.setattr(single_module, "run_pipeline", stub)
    single_module.run_single(["ERBB2"], run_name="my_run")
    assert captured["run_name"] == "my_run"
