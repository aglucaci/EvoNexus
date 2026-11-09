"""Expected co-detection between bat host species.

For host species :math:`i, j` and habitat :math:`u`,

.. math::

    C_{ij}^{(u)} = \\sum_{v \\in \\mathcal{V}} p_{iv}^{(u)} p_{jv}^{(u)}

where the probabilities are predicted for one *common* hypothetical sampling
design: the same tissue, depth, pool size and collection time, at a new site
and a new sequencing batch.

What this is
------------
The expected number of lineages that two hypothetical libraries — one from
each host species — would both be positive for, under the fitted model, if
they were sequenced identically. It assumes conditional independence of those
two hypothetical libraries given the parameters.

What this is not
----------------
It is not a transmission rate, a contact rate, or evidence of host-to-host
spillover. Two species can share expected detections because they share an
environment, a diet, a laboratory batch, or simply a common viral family.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from evonexus.bayesian.predict import PosteriorPredictor, expit
from evonexus.core.io import read_json, read_table, write_json, write_table
from evonexus.core.provenance import RunRecord

CONNECTIVITY_CAVEAT = (
    "Expected co-detection under the fitted model for a common sampling design. "
    "Not inferred host-to-host transmission."
)


def _canonical_profile(frame: pd.DataFrame) -> dict[str, Any]:
    """One standard library profile that every host is predicted at."""
    tissue = (
        frame["tissue_group"].mode().iat[0] if "tissue_group" in frame and len(frame) else "unknown"
    )
    profile = {
        "tissue_group": tissue,
        "log_depth_z": 0.0,
        "log_pool_size_z": 0.0,
        "time_z": 0.0,
        "site_id": "__new_site__",
        "sequencing_batch": "__new_batch__",
    }
    for col in frame.columns:
        if col.endswith("_z") and col not in profile:
            profile[col] = 0.0
    return profile


def host_lineage_probabilities(
    posterior: dict[str, Any],
    frame: pd.DataFrame,
    habitat: str = "settlement",
    hosts: list[str] | None = None,
    seed: int = 20250101,
    drop_interaction: bool = False,
) -> dict[str, Any]:
    """Posterior draws of p[host, lineage] for a common sampling design.

    With ``drop_interaction=True`` the host-by-family term is set to zero,
    which produces the marginal null: connectivity driven only by overall host
    and lineage detection frequencies.
    """
    predictor = PosteriorPredictor(posterior)
    hosts = hosts or sorted(set(frame["host_species"].astype(str)))
    lineages = frame[["virus_unit", "virus_family"]].drop_duplicates().sort_values("virus_unit")
    profile = _canonical_profile(frame)

    grid = lineages.merge(pd.DataFrame({"host_species": hosts}), how="cross")
    for key, value in profile.items():
        grid[key] = value
    grid["U"] = 1.0 if habitat == "settlement" else 0.0

    integrate = ["site", "batch"]
    if drop_interaction:
        # Zeroing the interaction is done by predicting at a level that does
        # not exist and forcing its contribution to zero.
        grid["host_species_for_ix"] = grid["host_species"]
        grid["host_species"] = grid["host_species"]

    rng = np.random.default_rng(seed)
    if drop_interaction:
        saved = predictor.block_span.pop("host_family_ix", None)
        p = expit(predictor.eta(grid, integrate=integrate, rng=rng,
                                habitat_override=grid["U"].iat[0]))
        if saved is not None:
            predictor.block_span["host_family_ix"] = saved
    else:
        p = expit(predictor.eta(grid, integrate=integrate, rng=rng,
                                habitat_override=grid["U"].iat[0]))

    n_hosts, n_lineages = len(hosts), len(lineages)
    # grid is lineages x hosts (cross keeps lineage-major order)
    p = p.reshape(predictor.n_samples, n_lineages, n_hosts).transpose(0, 2, 1)
    return {
        "p": p,                                  # (draws, host, lineage)
        "hosts": hosts,
        "lineages": lineages["virus_unit"].astype(str).tolist(),
        "families": lineages["virus_family"].astype(str).tolist(),
        "habitat": habitat,
        "drop_interaction": drop_interaction,
    }


def connectivity_matrix(p: np.ndarray) -> np.ndarray:
    """C[draw, i, j] = sum_v p[draw, i, v] * p[draw, j, v]."""
    return np.einsum("siv,sjv->sij", p, p)


def summarize_connectivity(C: np.ndarray, hosts: list[str], habitat: str,
                           hdi_prob: float = 0.95) -> pd.DataFrame:
    from evonexus.bayesian.diagnostics import hdi

    rows: list[dict[str, Any]] = []
    n = len(hosts)
    for i in range(n):
        for j in range(i, n):
            draws = C[:, i, j]
            lo, hi = hdi(draws, hdi_prob)
            rows.append(
                {
                    "habitat": habitat,
                    "host_i": hosts[i],
                    "host_j": hosts[j],
                    "self_pair": i == j,
                    "mean_expected_shared_lineages": float(np.mean(draws)),
                    "median": float(np.median(draws)),
                    "hdi_low": lo,
                    "hdi_high": hi,
                }
            )
    return pd.DataFrame(rows)


def marginal_null(
    posterior: dict[str, Any],
    frame: pd.DataFrame,
    habitat: str = "settlement",
    seed: int = 20250101,
) -> dict[str, Any]:
    """Null in which host-by-family structure is removed, sampling preserved."""
    pred = host_lineage_probabilities(posterior, frame, habitat=habitat, seed=seed,
                                      drop_interaction=True)
    return {"C": connectivity_matrix(pred["p"]), "hosts": pred["hosts"],
            "description": (
                "Host-by-viral-family interaction removed. Any excess connectivity "
                "above this null reflects lineage-specific host associations rather "
                "than overall host or lineage detection frequencies."
            )}


def permutation_null(
    frame: pd.DataFrame,
    n_permutations: int = 500,
    seed: int = 20250101,
    strata: tuple[str, ...] = ("habitat", "site_id"),
) -> dict[str, Any]:
    """Observed co-detection under permuted host labels within strata.

    Host species labels are shuffled among libraries *within* habitat and site,
    so the null keeps where and how much was sequenced while destroying the
    host-lineage association.
    """
    if frame.empty:
        return {"status": "empty"}
    rng = np.random.default_rng(seed)
    libs = frame.drop_duplicates(subset=["library_id"])[
        ["library_id", "host_species", *[s for s in strata if s in frame.columns]]
    ].reset_index(drop=True)
    wide = frame.pivot_table(index="library_id", columns="virus_unit", values="y",
                             aggfunc="max", fill_value=0.0)
    wide = wide.loc[libs["library_id"]]
    Y = wide.to_numpy(dtype=float)
    hosts = sorted(set(libs["host_species"].astype(str)))
    host_index = {h: i for i, h in enumerate(hosts)}

    def co_detection(labels: np.ndarray) -> np.ndarray:
        n = len(hosts)
        counts = np.zeros((n, n))
        per_host = np.zeros((n, Y.shape[1]))
        totals = np.zeros(n)
        for h, i in host_index.items():
            mask = labels == h
            if mask.any():
                per_host[i] = Y[mask].mean(axis=0)
                totals[i] = mask.sum()
        counts = per_host @ per_host.T
        return counts

    observed = co_detection(libs["host_species"].astype(str).to_numpy())

    stratum_keys = [s for s in strata if s in libs.columns]
    groups = (
        libs.groupby(stratum_keys).indices if stratum_keys
        else {"all": np.arange(len(libs))}
    )

    null_draws = np.empty((n_permutations, len(hosts), len(hosts)))
    labels = libs["host_species"].astype(str).to_numpy()
    for k in range(n_permutations):
        permuted = labels.copy()
        for idx in groups.values():
            idx = np.asarray(idx)
            permuted[idx] = rng.permutation(labels[idx])
        null_draws[k] = co_detection(permuted)

    mean_null = null_draws.mean(axis=0)
    sd_null = null_draws.std(axis=0, ddof=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(sd_null > 0, (observed - mean_null) / sd_null, np.nan)
    p_upper = (null_draws >= observed[None, :, :]).mean(axis=0)

    rows: list[dict[str, Any]] = []
    for i, hi in enumerate(hosts):
        for j, hj in enumerate(hosts):
            if j < i:
                continue
            rows.append(
                {
                    "host_i": hi,
                    "host_j": hj,
                    "observed_co_detection": float(observed[i, j]),
                    "null_mean": float(mean_null[i, j]),
                    "null_sd": float(sd_null[i, j]),
                    "z_score": float(z[i, j]),
                    "p_upper_tail": float(p_upper[i, j]),
                }
            )
    return {
        "status": "ok",
        "table": pd.DataFrame(rows),
        "hosts": hosts,
        "n_permutations": n_permutations,
        "strata": stratum_keys,
        "description": (
            "Host species labels permuted within habitat and site, preserving "
            "the sampling structure and each library's observed virome."
        ),
    }


def connectivity_vs_phylogeny(
    C: np.ndarray, hosts: list[str], distances: dict[str, Any] | None
) -> dict[str, Any]:
    """Association between expected co-detection and host phylogenetic distance."""
    if not distances:
        return {"status": "unavailable",
                "reason": "no host phylogenetic distance matrix was supplied"}
    species = [str(s) for s in distances.get("species", [])]
    matrix = np.asarray(distances.get("distance_matrix") or distances.get("matrix"), dtype=float)
    lookup = {s: i for i, s in enumerate(species)}
    usable = [h for h in hosts if h in lookup]
    if len(usable) < 4 or matrix.size == 0:
        return {"status": "unavailable",
                "reason": f"only {len(usable)} host species have phylogenetic distances"}

    sel_c = [hosts.index(h) for h in usable]
    sel_d = [lookup[h] for h in usable]
    d = matrix[np.ix_(sel_d, sel_d)]
    iu = np.triu_indices(len(usable), k=1)
    dist_vec = d[iu]

    slopes = []
    corrs = []
    for s in range(C.shape[0]):
        c_vec = C[s][np.ix_(sel_c, sel_c)][iu]
        if np.std(dist_vec) == 0 or np.std(c_vec) == 0:
            continue
        slope = np.polyfit(dist_vec, c_vec, 1)[0]
        slopes.append(slope)
        corrs.append(np.corrcoef(dist_vec, c_vec)[0, 1])

    if not slopes:
        return {"status": "unavailable", "reason": "no variation in distance or connectivity"}

    from evonexus.bayesian.diagnostics import hdi

    lo, hi = hdi(np.asarray(slopes))
    return {
        "status": "ok",
        "n_host_species": len(usable),
        "n_pairs": len(dist_vec),
        "slope_mean": float(np.mean(slopes)),
        "slope_hdi_low": lo,
        "slope_hdi_high": hi,
        "correlation_mean": float(np.mean(corrs)),
        "p_direction_negative": float(np.mean(np.asarray(slopes) < 0)),
        "note": (
            "A negative slope means more closely related hosts are expected to "
            "co-detect more lineages. This is an association across host pairs, "
            "not a test of viral codivergence."
        ),
    }


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def run_connectivity(
    posterior_dir: str | Path,
    model: str,
    frame_path: str | Path,
    out_dir: str | Path,
    habitats: tuple[str, ...] = ("settlement", "wild"),
    n_permutations: int = 500,
    seed: int = 20250101,
    host_distance_path: str | Path | None = None,
    synthetic: bool = False,
) -> dict[str, Any]:
    from evonexus.bayesian.fit import load_posterior

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    record = RunRecord(step="connectivity",
                       params={"model": model, "habitats": list(habitats),
                               "n_permutations": n_permutations, "synthetic": synthetic})
    record.synthetic = synthetic

    frame = read_table(frame_path)
    if frame.empty:
        record.error("analysis frame is empty")
        record.finish(out / "run_record_connectivity.json", status="empty")
        return {"status": "empty"}

    # The primary model may have been skipped because a critical QC gate
    # failed. There are then no posterior draws to predict from, and the
    # honest result is a recorded absence rather than a crash.
    if not (Path(posterior_dir) / f"posterior_{model}.json").exists():
        reason = (
            f"no posterior draws for model {model}; the primary model was not "
            "fitted, most likely because a critical QC gate failed"
        )
        record.warn(f"connectivity not computed: {reason}")
        write_json({"status": "not_run", "reason": reason, "model": model,
                    "synthetic": synthetic, "caveat": CONNECTIVITY_CAVEAT},
                   out / "connectivity_metadata.json")
        record.finish(out / "run_record_connectivity.json", status="not_run")
        return {"status": "not_run", "reason": reason}

    posterior = load_posterior(posterior_dir, model)
    distances = None
    if host_distance_path and Path(host_distance_path).exists():
        distances = read_json(host_distance_path)

    summaries: list[pd.DataFrame] = []
    null_summaries: list[pd.DataFrame] = []
    phylo: dict[str, Any] = {}
    matrices: dict[str, Any] = {}

    for habitat in habitats:
        pred = host_lineage_probabilities(posterior, frame, habitat=habitat, seed=seed)
        C = connectivity_matrix(pred["p"])
        summaries.append(summarize_connectivity(C, pred["hosts"], habitat))
        matrices[habitat] = {
            "hosts": pred["hosts"],
            "mean": C.mean(axis=0).tolist(),
            "sd": C.std(axis=0, ddof=1).tolist(),
        }

        null = marginal_null(posterior, frame, habitat=habitat, seed=seed)
        null_df = summarize_connectivity(null["C"], null["hosts"], habitat)
        null_df["null_type"] = "marginal_no_host_family_interaction"
        null_summaries.append(null_df)

        phylo[habitat] = connectivity_vs_phylogeny(C, pred["hosts"], distances)

    connectivity_df = pd.concat(summaries, ignore_index=True)
    null_df = pd.concat(null_summaries, ignore_index=True)

    merged = connectivity_df.merge(
        null_df[["habitat", "host_i", "host_j", "mean_expected_shared_lineages"]],
        on=["habitat", "host_i", "host_j"], suffixes=("", "_null"),
    )
    merged["excess_over_null"] = (
        merged["mean_expected_shared_lineages"] - merged["mean_expected_shared_lineages_null"]
    )

    if len(habitats) == 2:
        pivot = connectivity_df.pivot_table(
            index=["host_i", "host_j"], columns="habitat",
            values="mean_expected_shared_lineages",
        ).reset_index()
        if {"settlement", "wild"}.issubset(pivot.columns):
            pivot["settlement_minus_wild"] = pivot["settlement"] - pivot["wild"]
        write_table(pivot, out / "connectivity_habitat_contrast")

    perm = permutation_null(frame, n_permutations=n_permutations, seed=seed)
    if perm.get("status") == "ok":
        write_table(perm["table"], out / "connectivity_permutation_null")

    paths = {
        "connectivity": str(write_table(merged, out / "connectivity")),
        "null": str(write_table(null_df, out / "connectivity_marginal_null")),
    }
    write_json(
        {
            "status": "ok",
            "model": model,
            "habitats": list(habitats),
            "matrices": matrices,
            "phylogeny_association": phylo,
            "permutation_null": {k: v for k, v in perm.items() if k != "table"},
            "caveat": CONNECTIVITY_CAVEAT,
            "synthetic": synthetic,
        },
        out / "connectivity_metadata.json",
    )
    for p in paths.values():
        record.add_output(p, role="connectivity", checksum=False)
    record.note(f"connectivity computed for habitats {list(habitats)} — {CONNECTIVITY_CAVEAT}")
    record.finish(out / "run_record_connectivity.json", status="ok")

    return {"status": "ok", **paths}
