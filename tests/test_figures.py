"""QC gates, connectivity, validation and reporting tests.

The reporting tests exist mainly to pin down one behaviour: a missing analysis
must produce a visible "NOT GENERATED" artifact, never a plausible-looking
number.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from evonexus.connectivity.compute import (
    connectivity_matrix,
    host_lineage_probabilities,
    permutation_null,
)
from evonexus.core.config import load_config
from evonexus.core.io import write_table
from evonexus.evolution.hostphylo import newick_to_covariance, parse_newick
from evonexus.qc.gates import evaluate_gates
from evonexus.reporting.artifacts import ArtifactStore
from evonexus.reporting.figures import build_figures
from evonexus.reporting.legends import build_legends
from evonexus.reporting.tables import build_tables
from evonexus.validation.crossval import calibration_curve, scoring_rules


# ---------------------------------------------------------------------------
# QC gates
# ---------------------------------------------------------------------------
def test_all_gates_evaluate(tables):
    report = evaluate_gates(tables, load_config("qc_gates"))
    assert len(report.results) == 7
    assert {r.gate_id for r in report.results} == {f"G{i}" for i in range(1, 8)}
    for r in report.results:
        assert r.status in {"pass", "fail", "warn"}
        assert isinstance(r.metrics, dict)


def test_gate_failure_blocks_the_primary_model(tables):
    """Remove the habitat labels and G4 must fail, stopping the primary model."""
    broken = {k: v.copy() for k, v in tables.items()}
    broken["sites"]["habitat"] = "unknown"
    report = evaluate_gates(broken, load_config("qc_gates"))
    g4 = next(r for r in report.results if r.gate_id == "G4")
    assert not g4.passed and g4.critical
    assert report.primary_model_allowed is False


def test_gate_g3_fails_without_valid_negatives(tables):
    broken = {k: v.copy() for k, v in tables.items()}
    det = broken["detections"]
    det.loc[det["detection"] == "negative", "detection"] = "unassessable"
    det.loc[det["detection"] == "unassessable", "y"] = np.nan
    report = evaluate_gates(broken, load_config("qc_gates"))
    g3 = next(r for r in report.results if r.gate_id == "G3")
    assert not g3.passed


def test_gates_never_raise_on_empty_tables():
    empty = {k: pd.DataFrame() for k in
             ("sites", "hosts", "libraries", "viruses", "detections", "ecology", "crosswalk")}
    report = evaluate_gates(empty, load_config("qc_gates"))
    assert len(report.results) == 7
    assert report.primary_model_allowed is False


# ---------------------------------------------------------------------------
# connectivity
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def posterior(analysis_frame):
    from evonexus.bayesian.design import build_design
    from evonexus.bayesian.fit import fit_model
    from evonexus.core.config import resolve_variant

    frame = analysis_frame["frame"]
    terms = load_config("model_comparisons")["models"]["M3"]["terms"]
    design = build_design(frame, terms, model_name="M3")
    priors = resolve_variant(load_config("priors"), "regularizing")
    return fit_model(design, priors,
                     sampler={"draws": 100, "tune": 100, "chains": 2, "seed": 4},
                     engine="polya_gamma_gibbs")


def test_connectivity_matrix_is_symmetric_and_bounded(posterior, analysis_frame):
    frame = analysis_frame["frame"]
    pred = host_lineage_probabilities(posterior, frame, habitat="settlement", seed=2)
    C = connectivity_matrix(pred["p"])
    n_hosts = len(pred["hosts"])
    n_lineages = len(pred["lineages"])
    assert C.shape == (posterior["beta"].shape[0] * posterior["beta"].shape[1],
                       n_hosts, n_hosts)
    assert np.allclose(C, np.transpose(C, (0, 2, 1)))
    # C_ij is a sum of products of probabilities, so it cannot exceed n_lineages.
    assert (C >= 0).all() and (C <= n_lineages + 1e-9).all()


def test_marginal_null_differs_from_the_fitted_model(posterior, analysis_frame):
    frame = analysis_frame["frame"]
    full = host_lineage_probabilities(posterior, frame, habitat="settlement", seed=2)
    null = host_lineage_probabilities(posterior, frame, habitat="settlement", seed=2,
                                      drop_interaction=True)
    assert not np.allclose(full["p"], null["p"])


def test_permutation_null_preserves_sampling(analysis_frame):
    frame = analysis_frame["frame"]
    result = permutation_null(frame, n_permutations=25, seed=1)
    assert result["status"] == "ok"
    table = result["table"]
    assert {"observed_co_detection", "null_mean", "z_score", "p_upper_tail"}.issubset(table.columns)
    assert ((table["p_upper_tail"] >= 0) & (table["p_upper_tail"] <= 1)).all()


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------
def test_scoring_rules_reward_a_better_model():
    rng = np.random.default_rng(0)
    y = (rng.random(500) < 0.3).astype(float)
    good = np.tile(np.where(y > 0.5, 0.9, 0.1), (20, 1))
    bad = np.full((20, 500), 0.5)
    g, b = scoring_rules(y, good), scoring_rules(y, bad)
    assert g["elpd"] > b["elpd"]
    assert g["brier"] < b["brier"]
    assert g["auc"] > b["auc"]


def test_calibration_curve_is_well_formed():
    rng = np.random.default_rng(2)
    p = rng.random(1000)
    y = (rng.random(1000) < p).astype(float)
    curve = calibration_curve(y, p, n_bins=5)
    assert len(curve) <= 5
    assert curve["n"].sum() == 1000
    assert (curve["observed_frequency"].between(0, 1)).all()
    assert abs(curve["calibration_error"]).max() < 0.25


# ---------------------------------------------------------------------------
# host phylogeny
# ---------------------------------------------------------------------------
def test_newick_parsing_and_covariance():
    newick = "((Rhinolophus_sinicus:0.1,Rhinolophus_affinis:0.1):0.3,Myotis_davidii:0.4);"
    root = parse_newick(newick)
    assert len(root.children) == 2
    cov = newick_to_covariance(newick)
    assert cov["n_species"] == 3
    matrix = np.asarray(cov["matrix"])
    assert np.allclose(np.diag(matrix), 1.0)
    assert np.allclose(matrix, matrix.T)
    species = cov["species"]
    i, j = species.index("Rhinolophus sinicus"), species.index("Rhinolophus affinis")
    k = species.index("Myotis davidii")
    assert matrix[i, j] > matrix[i, k], "sister taxa must covary more"


def test_malformed_newick_raises():
    with pytest.raises(ValueError):
        parse_newick("((A:0.1,B:0.2)")


def test_missing_species_are_reported_not_imputed():
    cov = newick_to_covariance("(A:0.1,B:0.2);", species=["A", "B", "C"])
    assert cov["missing_species"] == ["C"]
    assert cov["n_species"] == 2


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
def test_missing_analyses_produce_not_generated_artifacts(tmp_path):
    """An empty results tree must still produce labelled, honest outputs."""
    store = ArtifactStore(tmp_path / "results")
    config = load_config("figures")
    out = tmp_path / "out"

    figures = build_figures(store, config, out)
    assert set(figures) == {"Figure1", "Figure2", "Figure3", "Figure4", "Figure5"}
    assert all(f["status"] != "ok" for f in figures.values())
    for result in figures.values():
        assert result["files"], "a placeholder figure must still be written"
        assert result.get("missing")

    tables = build_tables(store, config, out)
    assert all(t["status"] == "not_generated" for t in tables.values())
    assert all(t.get("reason") for t in tables.values())

    legends = build_legends(store, config, figures, out)
    for name, legend in legends.items():
        text = Path(legend["path"]).read_text()
        assert "NOT GENERATED" in text
        assert "No sample size" in text
        # Nothing that could be mistaken for a result.
        assert "95% HDI" not in text


def test_reporting_marks_synthetic_runs(tmp_path, tables, analysis_frame):
    results = tmp_path / "results"
    for name, df in tables.items():
        write_table(df, results / "processed" / name)
    write_table(analysis_frame["frame"], results / "matrix" / "analysis_frame_lineage")
    (results / "matrix").mkdir(parents=True, exist_ok=True)
    (results / "matrix" / "matrix_metadata_lineage.json").write_text(
        json.dumps({"synthetic": True, "summary": analysis_frame["summary"]})
    )
    store = ArtifactStore(results)
    assert store.is_synthetic() is True

    config = load_config("figures")
    figures = build_figures(store, config, tmp_path / "out")
    legends = build_legends(store, config, figures, tmp_path / "out")
    text = Path(legends["Figure1"]["path"]).read_text()
    assert "SYNTHETIC DATA" in text
