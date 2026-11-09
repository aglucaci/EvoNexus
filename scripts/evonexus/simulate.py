"""Deterministic synthetic data generator.

Two jobs:

1. Provide a laptop-scale dataset that exercises the *whole* pipeline,
   including metadata harmonization — so the files it writes deliberately use
   messy, source-like column names rather than the canonical schema.
2. Provide ground-truth parameters for parameter recovery, so the Bayesian
   machinery can be checked against a known answer.

Every output is stamped ``synthetic: true``. The reporting layer refuses to
present synthetic results as empirical ones and watermarks every figure.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from evonexus.core.io import write_json

SYNTHETIC_MARKER = "SYNTHETIC"

HOST_SPECIES = [
    "Rhinolophus sinicus",
    "Rhinolophus affinis",
    "Myotis davidii",
    "Hipposideros armiger",
    "Pipistrellus abramus",
    "Miniopterus fuliginosus",
    "Rousettus leschenaultii",
    "Taphozous melanopogon",
]
VIRUS_FAMILIES = [
    "Coronaviridae",
    "Paramyxoviridae",
    "Astroviridae",
    "Picornaviridae",
    "Caliciviridae",
    "Dicistroviridae",   # insect/diet-derived: must be filtered out downstream
]
TISSUES = ["gut", "lung", "oral swab", "liver", "spleen"]


def _truncated_normal(rng: np.random.Generator, mu: float, sigma: float,
                      size: int, lo: float, hi: float) -> np.ndarray:
    draws = rng.normal(mu, sigma, size)
    return np.clip(draws, lo, hi)


def simulate(
    out_dir: str | Path,
    n_sites: int = 12,
    n_libraries: int = 120,
    n_viruses: int = 30,
    seed: int = 20250101,
    unassessable_fraction: float = 0.08,
    mixed_pool_fraction: float = 0.10,
) -> dict[str, Any]:
    """Write source-like synthetic tables plus the generating truth."""
    rng = np.random.default_rng(seed)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ---------------- sites -------------------------------------------
    habitat_raw = rng.choice(
        ["urban", "village", "forest", "cave", "nature reserve", ""],
        size=n_sites,
        p=[0.22, 0.20, 0.22, 0.22, 0.10, 0.04],
    )
    sites = pd.DataFrame(
        {
            "Site_ID": [f"S{i:02d}" for i in range(1, n_sites + 1)],
            "Region": rng.choice(["Yunnan", "Guangdong", "Hainan", "Guangxi"], n_sites),
            "Country": "China",
            "Habitat": habitat_raw,
            "Latitude": np.round(_truncated_normal(rng, 23.0, 3.0, n_sites, 18.0, 29.0), 4),
            "Longitude": np.round(_truncated_normal(rng, 108.0, 5.0, n_sites, 98.0, 118.0), 4),
            "coordinate_precision": "site_centroid",
        }
    )
    # Guarantee at least one site with no recorded habitat: real deposits have
    # metadata gaps, and the "unknown is not a category" rule must be exercised.
    habitat_raw[-1] = ""
    sites.loc[sites.index[-1], "Habitat"] = ""
    settlement = np.isin(habitat_raw, ["urban", "village"]).astype(float)

    # ---------------- hosts / pools ------------------------------------
    n_hosts = n_libraries  # one source pool per library in this simulation
    host_site_idx = rng.integers(0, n_sites, n_hosts)
    # Host species composition depends on habitat, which creates exactly the
    # kind of host/habitat confounding the model has to standardize over.
    species_idx = np.empty(n_hosts, dtype=int)
    for i, s in enumerate(host_site_idx):
        if settlement[s] > 0:
            weights = np.array([0.08, 0.08, 0.10, 0.10, 0.26, 0.14, 0.16, 0.08])
        else:
            weights = np.array([0.20, 0.20, 0.14, 0.16, 0.06, 0.12, 0.06, 0.06])
        species_idx[i] = rng.choice(len(HOST_SPECIES), p=weights / weights.sum())

    is_pool = rng.random(n_hosts) < 0.55
    pool_size = np.where(is_pool, rng.integers(2, 9, n_hosts), 1).astype(float)
    mixed = is_pool & (rng.random(n_hosts) < mixed_pool_fraction)
    species_label = np.array([HOST_SPECIES[i] for i in species_idx], dtype=object)
    species_label[mixed] = "mixed"

    hosts = pd.DataFrame(
        {
            "pool_id": [f"P{i:04d}" for i in range(n_hosts)],
            "site_id": [sites["Site_ID"][s] for s in host_site_idx],
            "species": species_label,
            "pooled": np.where(is_pool, "pool", "individual"),
            "n_individuals": pool_size,
            "sampling_date": pd.to_datetime("2016-01-01")
            + pd.to_timedelta(rng.integers(0, 2200, n_hosts), unit="D"),
        }
    )
    hosts["sampling_date"] = hosts["sampling_date"].dt.strftime("%Y-%m-%d")

    # ---------------- libraries ----------------------------------------
    depth = np.exp(rng.normal(np.log(2.4e7), 0.55, n_libraries))
    libraries = pd.DataFrame(
        {
            "Library": [f"L{i:04d}" for i in range(n_libraries)],
            "pool_id": hosts["pool_id"].to_numpy(),
            "site_id": hosts["site_id"].to_numpy(),
            "Tissue": rng.choice(TISSUES, n_libraries),
            "clean_reads": np.round(depth).astype(np.int64),
            "non_host_non_rrna_reads": np.round(depth * rng.uniform(0.05, 0.28, n_libraries)).astype(np.int64),
            "sequencing_batch": rng.choice([f"B{j}" for j in range(1, 7)], n_libraries),
            "Run": [f"SRR{9000000 + i}" for i in range(n_libraries)],
            "BioProject": "SYNTHETIC-PRJ",
            "collection_date": hosts["sampling_date"].to_numpy(),
            "reference_panel": "panel_v1",
        }
    )

    # ---------------- viruses ------------------------------------------
    fam_idx = rng.integers(0, len(VIRUS_FAMILIES), n_viruses)
    viruses = pd.DataFrame(
        {
            "vOTU": [f"V{i:03d}" for i in range(n_viruses)],
            "Family": [VIRUS_FAMILIES[i] for i in fam_idx],
            "virus_name": [f"synthetic virus {i}" for i in range(n_viruses)],
            "assignment_method": "synthetic_ground_truth",
            "accession": [f"SYN{i:05d}" for i in range(n_viruses)],
        }
    )

    # ---------------- generating parameters ----------------------------
    truth: dict[str, Any] = {
        "alpha": -3.4,
        "beta_habitat": 0.55,          # the estimand the pipeline must recover
        "gamma_depth": 0.45,
        "gamma_poolsize": 0.30,
        "gamma_time": -0.08,
        "sigma_family": 0.8,
        "sigma_lineage": 0.6,
        "sigma_host": 0.5,
        "sigma_tissue": 0.4,
        "sigma_site": 0.35,
        "sigma_interaction": 0.45,
    }
    a_family = rng.normal(0, truth["sigma_family"], len(VIRUS_FAMILIES))
    a_lineage = rng.normal(0, truth["sigma_lineage"], n_viruses)
    b_host = rng.normal(0, truth["sigma_host"], len(HOST_SPECIES))
    c_tissue = rng.normal(0, truth["sigma_tissue"], len(TISSUES))
    g_site = rng.normal(0, truth["sigma_site"], n_sites)
    kappa = rng.normal(0, truth["sigma_interaction"], (len(HOST_SPECIES), len(VIRUS_FAMILIES)))

    eff = libraries["non_host_non_rrna_reads"].to_numpy(dtype=float)
    log_depth_z = (np.log(eff) - np.log(eff).mean()) / np.log(eff).std()
    log_pool_z = np.log(pool_size)
    log_pool_z = (log_pool_z - log_pool_z.mean()) / (log_pool_z.std() or 1.0)
    time_days = (pd.to_datetime(libraries["collection_date"]) - pd.Timestamp("2016-01-01")).dt.days.to_numpy()
    time_z = (time_days - time_days.mean()) / (time_days.std() or 1.0)
    tissue_idx = np.array([TISSUES.index(t) for t in libraries["Tissue"]])
    lib_site_idx = host_site_idx
    U = settlement[lib_site_idx]

    rows = []
    for li in range(n_libraries):
        # A mixed-species pool has no single host effect; the simulation uses
        # the mean, and the pipeline must exclude it from the single-host model.
        if mixed[li]:
            b_h = float(b_host.mean())
            k_h = kappa.mean(axis=0)
        else:
            b_h = float(b_host[species_idx[li]])
            k_h = kappa[species_idx[li]]
        for vi in range(n_viruses):
            f = fam_idx[vi]
            eta = (
                truth["alpha"]
                + a_family[f]
                + a_lineage[vi]
                + truth["beta_habitat"] * U[li]
                + b_h
                + k_h[f]
                + c_tissue[tissue_idx[li]]
                + g_site[lib_site_idx[li]]
                + truth["gamma_depth"] * log_depth_z[li]
                + truth["gamma_poolsize"] * log_pool_z[li]
                + truth["gamma_time"] * time_z[li]
            )
            p = 1.0 / (1.0 + np.exp(-eta))
            y = int(rng.random() < p)
            if y:
                mapped = int(rng.integers(12, 4000))
                breadth = float(np.clip(rng.beta(3, 2), 0.1, 1.0))
            else:
                mapped = int(rng.integers(0, 9))
                breadth = float(np.clip(rng.beta(1, 25), 0.0, 0.12))
            rows.append((libraries["Library"][li], viruses["vOTU"][vi], mapped, breadth,
                         mapped / eff[li] * 1e6))

    detections = pd.DataFrame(
        rows, columns=["Library", "vOTU", "reads_mapped", "coverage", "RPKM"]
    )

    # A realistic deposit does not report every cell. Drop a fraction of rows
    # and mark their libraries as not screened for those lineages, so the
    # "missing is not a zero" logic has something to act on.
    n_drop = int(len(detections) * unassessable_fraction)
    drop_idx = rng.choice(len(detections), n_drop, replace=False)
    not_screened = detections.iloc[drop_idx][["Library", "vOTU"]].copy()
    detections = detections.drop(index=detections.index[drop_idx]).reset_index(drop=True)

    # ---------------- ecology ------------------------------------------
    ecology = pd.DataFrame(
        {
            "site_id": sites["Site_ID"],
            "forest_cover": np.round(np.clip(0.75 - 0.45 * settlement + rng.normal(0, 0.10, n_sites), 0, 1), 3),
            "relative_humidity": np.round(_truncated_normal(rng, 76, 6, n_sites, 50, 95), 2),
            "mean_temp": np.round(_truncated_normal(rng, 21.5, 2.5, n_sites, 12, 30), 2),
            "precipitation": np.round(_truncated_normal(rng, 1450, 320, n_sites, 500, 2600), 1),
        }
    )

    # ---------------- write source-like files --------------------------
    paths = {
        "sites": out / "metadata_sites.csv",
        "hosts": out / "metadata_pools.csv",
        "libraries": out / "metadata_libraries.csv",
        "viruses": out / "viral_catalog.csv",
        "detections": out / "detection_table.csv",
        "ecology": out / "site_environment.csv",
        "not_screened": out / "not_screened_cells.csv",
    }
    sites.to_csv(paths["sites"], index=False)
    hosts.to_csv(paths["hosts"], index=False)
    libraries.to_csv(paths["libraries"], index=False)
    viruses.to_csv(paths["viruses"], index=False)
    detections.to_csv(paths["detections"], index=False)
    ecology.to_csv(paths["ecology"], index=False)
    not_screened.to_csv(paths["not_screened"], index=False)

    truth_payload = {
        "synthetic": True,
        "marker": SYNTHETIC_MARKER,
        "seed": seed,
        "dimensions": {
            "n_sites": n_sites,
            "n_libraries": n_libraries,
            "n_viruses": n_viruses,
            "n_host_species": len(HOST_SPECIES),
            "n_cells_written": len(detections),
            "n_cells_not_screened": len(not_screened),
            "n_mixed_species_pools": int(mixed.sum()),
        },
        "parameters": truth,
        "random_effects": {
            "a_family": dict(zip(VIRUS_FAMILIES, map(float, a_family))),
            "b_host": dict(zip(HOST_SPECIES, map(float, b_host))),
            "c_tissue": dict(zip(TISSUES, map(float, c_tissue))),
        },
        "note": (
            "Synthetic data generated for software validation only. "
            "These numbers describe no real bat, virus, place or study."
        ),
    }
    write_json(truth_payload, out / "synthetic_truth.json")

    # A mapping file that is already approved for the synthetic source tables,
    # so the test profile can run unattended while the real mapping still
    # requires human approval.
    write_json(
        {
            "approved": True,
            "approved_by": "synthetic-fixture",
            "entities": {
                "sites": {"source_table": paths["sites"].name},
                "hosts": {"source_table": paths["hosts"].name},
                "libraries": {"source_table": paths["libraries"].name},
                "viruses": {"source_table": paths["viruses"].name},
                "detections": {"source_table": paths["detections"].name},
                "ecology": {"source_table": paths["ecology"].name},
            },
            "extras": {"not_screened_table": paths["not_screened"].name},
        },
        out / "mapping_overrides.json",
    )

    return {"paths": {k: str(v) for k, v in paths.items()}, "truth": truth_payload}
