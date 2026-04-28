"""End-to-end smoke test on a small fixture target.

Runs the full pipeline (target → structure → membrane → spatial → surface →
conservation → specificity → scoring → export → viz) on one cached UniProt ID
and asserts the run produced *some* output. Network-dependent on first run;
cached afterwards. Skipped by default with the "slow" marker — run with
`pytest -m slow` or `pytest tests/test_pipeline_smoke.py`.
"""
import pytest

from epitope_pipeline.run import run_pipeline


@pytest.mark.slow
def test_pipeline_runs_end_to_end(tmp_path, monkeypatch):
    """Smoke: ERBB2 produces a populated metrics dict."""
    # Redirect runs to the temporary directory so we don't pollute repo state
    from epitope_pipeline import config as cfg
    monkeypatch.setattr(cfg, "RUNS_DIR", tmp_path / "runs")

    result = run_pipeline(["P04626"], verbose=False)
    assert "targets" in result
    assert len(result["targets"]) == 1
    metrics = result["metrics"][result["targets"][0].uniprot_id]
    # Run completed; we don't care what the score is, just that we got one
    assert "n_patches" in metrics
    assert "best_score" in metrics
