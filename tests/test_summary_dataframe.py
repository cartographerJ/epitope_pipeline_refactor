"""Tests for build_summary_dataframe and write_summary_metrics_csv."""

import math
from types import SimpleNamespace

import pandas as pd
import pytest

from epitope_pipeline.io.export import (
    SUMMARY_DATAFRAME_COLUMNS,
    build_summary_dataframe,
    write_summary_metrics_csv,
)


def _target(uid="P0001", gene="ACME", seq_len=500):
    return SimpleNamespace(
        uniprot_id=uid,
        gene_name=gene,
        sequence_length=seq_len,
    )


def _membrane(topology="single_pass", n_tm=1):
    return SimpleNamespace(
        topology_type=topology,
        tm_segments=[("a", "b")] * n_tm,
    )


def _surface(n_patches=3):
    return SimpleNamespace(patches=list(range(n_patches)))


def _conservation(n_conserved=2, overall=0.95):
    return SimpleNamespace(
        conserved_patches=list(range(n_conserved)),
        overall_identity=overall,
    )


def _specificity(n_specific=1):
    return SimpleNamespace(specific_patches=list(range(n_specific)))


def _score(max_off=0.45):
    return SimpleNamespace(max_off_target_identity=max_off)


def _metric_full():
    return {
        "total_epitope_area_a2": 1200.0,
        "total_ectodomain_area_a2": 8000.0,
        "epitope_fraction": 0.15,
        "n_patches": 2,
        "best_score": 0.72,
        "target_score": 0.42,
        "area_component": 0.24,
        "quality_component": 0.93,
    }


PARAMETERS = {
    "mode": "distal (min_distance_a=80.0)",
    "cyno_mismatch_percent_base": 15.0,
    "cyno_gate_skipped": False,
}


def test_columns_match_spec():
    df = build_summary_dataframe(
        targets=[_target()],
        target_metrics={"P0001": _metric_full()},
        surface_analyses={"P0001": _surface()},
        conservation_results={"P0001": _conservation()},
        specificity_results={"P0001": _specificity()},
        epitope_scores={"P0001": [_score()]},
        membranes={"P0001": _membrane()},
        parameters=PARAMETERS,
        run_name="run_abc",
    )
    assert list(df.columns) == SUMMARY_DATAFRAME_COLUMNS


def test_row_per_target_including_failures():
    """Even targets that failed every step produce a row with NaN/0 sentinels."""
    targets = [_target("P1", "OK"), _target("P2", "FAIL")]
    df = build_summary_dataframe(
        targets=targets,
        target_metrics={"P1": _metric_full()},   # P2 absent → empty_metric defaults
        surface_analyses={"P1": _surface(4)},    # P2 missing
        conservation_results={"P1": _conservation(n_conserved=3, overall=0.91)},
        specificity_results={"P1": _specificity(2)},
        epitope_scores={"P1": [_score(0.5), _score(0.3)]},
        membranes={"P1": _membrane("single_pass", 1)},
        parameters=PARAMETERS,
        run_name="r",
    )
    assert len(df) == 2
    assert df.loc[0, "gene_name"] == "OK"
    assert df.loc[1, "gene_name"] == "FAIL"

    ok = df.loc[0]
    assert ok["n_surface_patches"] == 4
    assert ok["n_conserved_patches"] == 3
    assert ok["n_specific_patches"] == 2
    assert ok["n_scored_patches"] == 2
    assert ok["max_off_target_identity_overall"] == pytest.approx(0.5)
    assert ok["overall_cyno_identity"] == pytest.approx(0.91)
    assert ok["topology"] == "single_pass"

    fail = df.loc[1]
    assert fail["n_surface_patches"] == 0
    assert fail["n_conserved_patches"] == 0
    assert fail["n_specific_patches"] == 0
    assert fail["n_scored_patches"] == 0
    assert fail["n_patches"] == 0
    assert math.isnan(fail["total_epitope_area_a2"])
    assert math.isnan(fail["target_score"])
    assert math.isnan(fail["overall_cyno_identity"])
    assert math.isnan(fail["max_off_target_identity_overall"])
    assert fail["topology"] == ""


def test_empty_metric_keys_become_nan_not_keyerror():
    """target_score / area_component / quality_component are absent in empty_metric()."""
    targets = [_target("P1")]
    minimal_metric = {
        "total_epitope_area_a2": 0.0,
        "total_ectodomain_area_a2": 0.0,
        "epitope_fraction": 0.0,
        "n_patches": 0,
        "best_score": 0.0,
        # NOTE: no target_score / area_component / quality_component
    }
    df = build_summary_dataframe(
        targets=targets,
        target_metrics={"P1": minimal_metric},
        surface_analyses={},
        conservation_results={},
        specificity_results={},
        epitope_scores={},
        membranes={},
        parameters=PARAMETERS,
        run_name="r",
    )
    row = df.loc[0]
    assert math.isnan(row["target_score"])
    assert math.isnan(row["area_component"])
    assert math.isnan(row["quality_component"])
    assert row["best_score"] == 0.0
    assert row["n_patches"] == 0


def test_parameters_are_passed_through():
    params = {
        "mode": "proximal (max_distance_a=40)",
        "cyno_mismatch_percent_base": 25.0,
        "cyno_gate_skipped": True,
    }
    df = build_summary_dataframe(
        targets=[_target("P1")],
        target_metrics={"P1": _metric_full()},
        surface_analyses={"P1": _surface()},
        conservation_results={"P1": _conservation()},
        specificity_results={"P1": _specificity()},
        epitope_scores={"P1": [_score()]},
        membranes={"P1": _membrane()},
        parameters=params,
        run_name="run_xyz",
    )
    row = df.loc[0]
    assert row["run_name"] == "run_xyz"
    assert row["mode"] == "proximal (max_distance_a=40)"
    assert row["cyno_mismatch_percent_used"] == 25.0
    assert bool(row["cyno_gate_skipped"]) is True


def test_csv_round_trip(tmp_path):
    df = build_summary_dataframe(
        targets=[_target("P1", "A"), _target("P2", "B")],
        target_metrics={"P1": _metric_full(), "P2": _metric_full()},
        surface_analyses={"P1": _surface(), "P2": _surface()},
        conservation_results={"P1": _conservation(), "P2": _conservation()},
        specificity_results={"P1": _specificity(), "P2": _specificity()},
        epitope_scores={"P1": [_score()], "P2": [_score()]},
        membranes={"P1": _membrane(), "P2": _membrane()},
        parameters=PARAMETERS,
        run_name="r",
    )
    out = write_summary_metrics_csv(tmp_path, df)
    assert out.exists()
    assert out.name == "summary_metrics.csv"
    loaded = pd.read_csv(out)
    assert list(loaded.columns) == SUMMARY_DATAFRAME_COLUMNS
    assert len(loaded) == 2
    assert loaded.loc[0, "gene_name"] == "A"
    assert loaded.loc[1, "gene_name"] == "B"
