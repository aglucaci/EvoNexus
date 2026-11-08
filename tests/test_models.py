"""Model construction, sampling and estimand tests."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from evonexus.bayesian.checks import prior_predictive_check
from evonexus.bayesian.contrast import habitat_contrast, standardized_predictions
from evonexus.bayesian.design import build_design, design_matrix
from evonexus.bayesian.diagnostics import (
    compare_models,
    convergence_report,
    ess,
    hdi,
    psis_loo,
    rhat,
)
from evonexus.bayesian.fit import fit_model
from evonexus.bayesian.predict import PosteriorPredictor, expit
from evonexus.core.config import load_config, resolve_variant
from evonexus.quantify.matrix import apply_restriction

SAMPLER = {"draws": 120, "tune": 120, "chains": 2, "seed": 11}


@pytest.fixture(scope="module")
def priors():
    return resolve_variant(load_config("priors"), "regularizing")


@pytest.fixture(scope="module")
def terms():
    return load_config("model_comparisons")["models"]["M3"]["terms"]


@pytest.fixture(scope="module")
def frame(analysis_frame):
    return analysis_frame["frame"]


@pytest.fixture(scope="module")
def fitted(frame, terms, priors):
    design = build_design(frame, terms, model_name="M3")
    return design, fit_model(design, priors, sampler=SAMPLER, engine="polya_gamma_gibbs")


# ---------------------------------------------------------------------------
# analysis frame
# ---------------------------------------------------------------------------
def test_analysis_frame_has_no_missing_responses(frame):
    assert len(frame) > 0
    assert frame["y"].notna().all()
    assert set(frame["y"].unique()).issubset({0.0, 1.0})


def test_analysis_frame_excludes_mixed_pools_and_unknown_habitat(frame):
    assert set(frame["habitat"]).issubset({"settlement", "wild"})
    assert "mixed" not in set(frame["host_species"])
    assert "unknown" not in set(frame["host_species"])


def test_exclusion_ledger_is_monotone(analysis_frame):
    counts = [row["cells_remaining"] for row in analysis_frame["ledger"]]
    assert counts == sorted(counts, reverse=True), "filters must only remove cells"


def test_standardization_constants_are_recorded(analysis_frame):
    std = analysis_frame["standardization"]
    assert {"log_depth_z", "log_pool_size_z", "time_z"}.issubset(std)
    for spec in std.values():
        assert spec["sd"] > 0


@pytest.mark.parametrize("restriction",
                         ["none", "overlap_only", "no_pools", "drop_rare_lineages",
                          "depth_matched"])
def test_restrictions_return_subsets(frame, tables, restriction):
    sub = apply_restriction(frame, restriction, tables)
    assert len(sub) <= len(frame)
    assert set(sub.columns) == set(frame.columns)


# ---------------------------------------------------------------------------
# design
# ---------------------------------------------------------------------------
def test_design_blocks_match_the_model_specification(frame, terms):
    design = build_design(frame, terms, model_name="M3")
    names = {b.name for b in design.blocks}
    assert {"virus_family", "virus_lineage", "host", "tissue", "site",
            "host_family_ix"}.issubset(names)
    assert "habitat" in design.fixed_names
    assert design.has_intercept


def test_design_matrix_is_consistent(frame, terms):
    design = build_design(frame, terms, model_name="M3")
    X, names, spans = design_matrix(design)
    assert X.shape == (design.n_obs, design.parameter_count())
    assert len(names) == X.shape[1]
    for label, start, stop in spans:
        if label.startswith("block:"):
            block = design.block(label.split(":", 1)[1])
            # Each observation loads on exactly one level of a block.
            assert np.allclose(X[:, start:stop].sum(axis=1), 1.0)
            assert stop - start == block.size


def test_design_refuses_missing_responses(frame, terms):
    bad = frame.copy()
    bad.loc[bad.index[0], "y"] = np.nan
    with pytest.raises(ValueError, match="missing responses"):
        build_design(bad, terms)


def test_model_ladder_is_nested(frame):
    comparisons = load_config("model_comparisons")["models"]
    for lower, upper in (("M0", "M1"), ("M1", "M2"), ("M2", "M3")):
        assert set(comparisons[lower]["terms"]).issubset(set(comparisons[upper]["terms"]))


# ---------------------------------------------------------------------------
# priors
# ---------------------------------------------------------------------------
def test_prior_predictive_rates_are_plausible(frame, terms, priors):
    design = build_design(frame, terms)
    check = prior_predictive_check(design, priors, draws=120, seed=3)
    assert 0.0 <= check["detection_rate_quantiles"]["50%"] <= 1.0
    assert check["fraction_in_plausible_band"] > 0.5


# ---------------------------------------------------------------------------
# sampling
# ---------------------------------------------------------------------------
def test_sampler_returns_the_expected_shapes(fitted):
    design, posterior = fitted
    beta = posterior["beta"]
    assert beta.shape == (SAMPLER["chains"], SAMPLER["draws"], design.parameter_count())
    assert posterior["sigma"].shape[-1] == len(design.blocks)
    assert posterior["log_likelihood"].shape[-1] == design.n_obs
    assert posterior["engine"] == "polya_gamma_gibbs"


def test_posterior_carries_predictor_metadata(fitted):
    _, posterior = fitted
    # Without these, predictions silently drop every fixed effect.
    assert posterior["fixed_names"]
    assert posterior["blocks"]
    assert "habitat" in posterior["fixed_names"]


def test_sampling_is_reproducible(frame, terms, priors):
    design = build_design(frame, terms)
    a = fit_model(design, priors, sampler=SAMPLER, engine="polya_gamma_gibbs")
    b = fit_model(design, priors, sampler=SAMPLER, engine="polya_gamma_gibbs")
    assert np.allclose(a["beta"], b["beta"])


def test_diagnostics_report_runs(fitted):
    _, posterior = fitted
    report = convergence_report(posterior, load_config("priors")["convergence_gates"])
    assert {"max_rhat", "min_ess_bulk", "engine"}.issubset(report["metrics"])
    assert len(report["summary"]) == posterior["beta"].shape[-1] + posterior["sigma"].shape[-1]


def test_rhat_and_ess_behave_on_known_input():
    rng = np.random.default_rng(0)
    independent = rng.standard_normal((4, 500))
    assert rhat(independent) < 1.05
    assert ess(independent) > 400
    # Chains stuck at different values must be flagged.
    stuck = np.stack([np.full(500, -3.0), np.full(500, 3.0)])
    assert rhat(stuck) > 1.5 or np.isnan(rhat(stuck))


def test_hdi_covers_the_bulk():
    rng = np.random.default_rng(1)
    draws = rng.normal(2.0, 1.0, 20000)
    low, high = hdi(draws, 0.95)
    assert low < 2.0 < high
    assert 3.5 < (high - low) < 4.5


# ---------------------------------------------------------------------------
# prediction and the estimand
# ---------------------------------------------------------------------------
def test_habitat_term_changes_predictions(fitted, frame):
    _, posterior = fitted
    predictor = PosteriorPredictor(posterior)
    eta1 = predictor.eta(frame.head(50), integrate=(), habitat_override=1.0,
                         rng=np.random.default_rng(0))
    eta0 = predictor.eta(frame.head(50), integrate=(), habitat_override=0.0,
                         rng=np.random.default_rng(0))
    assert not np.allclose(eta1, eta0), "habitat override had no effect on the predictor"
    assert np.allclose(eta1 - eta0, predictor.fixed("habitat")[:, None])


def test_unseen_levels_are_integrated_over(fitted, frame):
    _, posterior = fitted
    predictor = PosteriorPredictor(posterior)
    profiles = frame.head(20).copy()
    profiles["site_id"] = "__a_site_never_sampled__"
    contribution = predictor.block_contribution(
        "site", profiles["site_id"].tolist(), rng=np.random.default_rng(0)
    )
    assert contribution.std() > 0, "an unseen site must draw from its population"


def test_standardized_contrast_is_coherent(fitted, frame):
    _, posterior = fitted
    std = standardized_predictions(posterior, frame, max_libraries=20)
    assert std["status"] == "ok"
    assert np.allclose(std["difference"], std["p_settlement"] - std["p_wild"])
    assert ((std["p_settlement"] >= 0) & (std["p_settlement"] <= 1)).all()
    # The reference population is recorded so the estimand is reproducible.
    assert std["reference"]["n_library_profiles"] > 0


def test_habitat_contrast_table_carries_the_caveat(fitted, frame):
    _, posterior = fitted
    result = habitat_contrast(posterior, frame, max_libraries=15)
    assert result["status"] == "ok"
    assert "not a causal effect" in result["caveat"].lower()
    table = result["table"]
    assert len(table) >= 3
    row = table[table["quantity"].str.contains("difference")].iloc[0]
    assert row["hdi_95_low"] <= row["mean"] <= row["hdi_95_high"]


# ---------------------------------------------------------------------------
# information criteria
# ---------------------------------------------------------------------------
def test_loo_and_comparison(frame, priors):
    comparisons = load_config("model_comparisons")["models"]
    loos = {}
    for name in ("M0", "M3"):
        design = build_design(frame, comparisons[name]["terms"], model_name=name)
        posterior = fit_model(design, priors, sampler=SAMPLER, engine="polya_gamma_gibbs")
        loos[name] = psis_loo(posterior["log_likelihood"])
        assert loos[name]["elpd_loo"] < 0
    table = compare_models(loos)
    assert list(table["rank"]) == [1, 2]
    assert table.iloc[0]["elpd_diff"] == 0.0


# ---------------------------------------------------------------------------
# parameter recovery on a known generating process
# ---------------------------------------------------------------------------
def test_recovers_a_known_logistic_coefficient():
    """End-to-end check that the sampler is not biased on a simple case."""
    rng = np.random.default_rng(5)
    n = 4000
    x = rng.standard_normal(n)
    group = rng.integers(0, 6, n)
    offsets = rng.normal(0, 0.5, 6)
    eta = -1.0 + 1.25 * x + offsets[group]
    y = (rng.random(n) < expit(eta)).astype(float)

    frame = pd.DataFrame(
        {
            "y": y,
            "U": x,                      # reuse the habitat slot as the known covariate
            "virus_unit": "V",
            "virus_family": "F",
            "host_species": "H",
            "tissue_group": "T",
            "site_id": [f"S{g}" for g in group],
            "sequencing_batch": "B",
            "log_depth_z": 0.0,
            "log_pool_size_z": 0.0,
            "time_z": 0.0,
        }
    )
    design = build_design(frame, ["intercept", "habitat", "site"], model_name="known")
    posterior = fit_model(design, {"intercept": {"dist": "normal", "mu": 0, "sigma": 2},
                                   "habitat_beta": {"dist": "normal", "mu": 0, "sigma": 2},
                                   "site_sd": {"dist": "halfnormal", "sigma": 1}},
                          sampler={"draws": 400, "tune": 400, "chains": 2, "seed": 9},
                          engine="polya_gamma_gibbs")
    predictor = PosteriorPredictor(posterior)
    beta = predictor.fixed("habitat")
    low, high = hdi(beta, 0.95)
    assert low <= 1.25 <= high, f"true 1.25 outside HDI [{low:.3f}, {high:.3f}]"
    assert abs(float(np.mean(beta)) - 1.25) < 0.25
