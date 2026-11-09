"""The seven mandatory data gates.

Each gate answers one question that, if answered wrongly, would invalidate the
habitat estimand. Gates report numbers, not opinions; the thresholds live in
``config/qc_gates.yaml`` so that a reviewer can see and change them.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import pandas as pd

Tables = dict[str, pd.DataFrame]


@dataclass
class GateResult:
    gate_id: str
    name: str
    condition: str
    critical: bool
    passed: bool
    metrics: dict[str, Any] = field(default_factory=dict)
    failures: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def status(self) -> str:
        if self.passed:
            return "pass"
        return "fail" if self.critical else "warn"

    def to_dict(self) -> dict[str, Any]:
        return {
            "gate_id": self.gate_id,
            "name": self.name,
            "condition": self.condition,
            "critical": self.critical,
            "status": self.status,
            "passed": self.passed,
            "metrics": self.metrics,
            "failures": self.failures,
            "notes": self.notes,
        }


@dataclass
class GateReport:
    results: list[GateResult]
    on_critical_failure: dict[str, Any] = field(default_factory=dict)

    @property
    def critical_failures(self) -> list[GateResult]:
        return [r for r in self.results if r.critical and not r.passed]

    @property
    def warnings(self) -> list[GateResult]:
        return [r for r in self.results if not r.critical and not r.passed]

    @property
    def primary_model_allowed(self) -> bool:
        return not self.critical_failures

    def to_dict(self) -> dict[str, Any]:
        return {
            "primary_model_allowed": self.primary_model_allowed,
            "n_critical_failures": len(self.critical_failures),
            "n_warnings": len(self.warnings),
            "gates": [r.to_dict() for r in self.results],
            "on_critical_failure": self.on_critical_failure,
        }


def _frac(numerator: int, denominator: int) -> float:
    return float(numerator) / float(denominator) if denominator else 0.0


# ---------------------------------------------------------------------------
# individual gates
# ---------------------------------------------------------------------------
def _gate_g1(t: Tables, th: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    det, lib, vir, host = t["detections"], t["libraries"], t["viruses"], t["hosts"]
    known_libs = set(lib["library_id"].astype(str)) if len(lib) else set()
    known_vir = set(vir["virus_id"].astype(str)) if len(vir) else set()
    known_hosts = set(host["host_id"].astype(str)) if len(host) else set()
    known_sites = set(t["sites"]["site_id"].astype(str)) if len(t["sites"]) else set()

    orphan_det = 0
    if len(det):
        bad_lib = ~det["library_id"].astype(str).isin(known_libs)
        bad_vir = ~det["virus_id"].astype(str).isin(known_vir)
        orphan_det = int((bad_lib | bad_vir).sum())
    lib_no_host = int((~lib["host_id"].astype(str).isin(known_hosts)).sum()) if len(lib) else 0
    host_no_site = int((~host["site_id"].astype(str).isin(known_sites)).sum()) if len(host) else 0

    metrics = {
        "n_detection_cells": len(det),
        "orphan_detection_cells": orphan_det,
        "orphan_detection_fraction": _frac(orphan_det, len(det)),
        "libraries_without_host": lib_no_host,
        "libraries_without_host_fraction": _frac(lib_no_host, len(lib)),
        "hosts_without_site": host_no_site,
        "hosts_without_site_fraction": _frac(host_no_site, len(host)),
    }
    failures = []
    if metrics["orphan_detection_fraction"] > th.get("max_orphan_detection_fraction", 0.0):
        failures.append(
            f"{orphan_det} detection cells reference an unknown library or lineage"
        )
    if metrics["libraries_without_host_fraction"] > th.get("max_library_without_host_fraction", 0.02):
        failures.append(f"{lib_no_host} libraries have no resolvable source pool or animal")
    if metrics["hosts_without_site_fraction"] > th.get("max_host_without_site_fraction", 0.05):
        failures.append(f"{host_no_site} host records have no resolvable sampling site")
    return metrics, failures


def _gate_g2(t: Tables, th: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    host = t["hosts"]
    n = len(host)
    unknown_pool = int((host["pool_status"].astype(str) == "unknown").sum()) if n else 0
    unknown_species = int((host["host_species_resolved"].astype(str) == "unknown").sum()) if n else 0
    mixed = int(host["mixed_species_pool"].fillna(False).sum()) if n else 0
    sized = int(host["pool_size"].notna().sum()) if n else 0
    metrics = {
        "n_hosts": n,
        "unknown_pool_status": unknown_pool,
        "unknown_pool_status_fraction": _frac(unknown_pool, n),
        "unknown_host_species": unknown_species,
        "unknown_host_species_fraction": _frac(unknown_species, n),
        "mixed_species_pools": mixed,
        "pool_size_documented": sized,
        "pool_size_documented_fraction": _frac(sized, n),
    }
    failures = []
    if metrics["unknown_pool_status_fraction"] > th.get("max_unknown_pool_status_fraction", 0.2):
        failures.append(
            f"{unknown_pool}/{n} host records do not state whether they are a pool or one animal"
        )
    if metrics["unknown_host_species_fraction"] > th.get("max_unknown_host_species_fraction", 0.1):
        failures.append(f"{unknown_species}/{n} host records have no resolvable species")
    return metrics, failures


def _gate_g3(t: Tables, th: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    det = t["detections"]
    n = len(det)
    counts = det["detection"].value_counts().to_dict() if n else {}
    n_pos = int(counts.get("positive", 0))
    n_neg = int(counts.get("negative", 0))
    n_un = int(counts.get("unassessable", 0))
    metrics = {
        "n_cells": n,
        "n_positive": n_pos,
        "n_valid_negative": n_neg,
        "n_unassessable": n_un,
        "valid_negative_fraction": _frac(n_neg, n),
        "unassessable_fraction": _frac(n_un, n),
        "marginal_detection_rate_over_assessable": _frac(n_pos, n_pos + n_neg),
    }
    failures = []
    if n_pos < th.get("min_positive_cells", 20):
        failures.append(f"only {n_pos} positive cells; too few to model lineage detection")
    if metrics["valid_negative_fraction"] < th.get("min_valid_negative_fraction", 0.5):
        failures.append(
            f"only {metrics['valid_negative_fraction']:.1%} of cells are assay-tested "
            "negatives; absence cannot be treated as a zero for the rest"
        )
    if metrics["unassessable_fraction"] > th.get("max_unassessable_fraction", 0.4):
        failures.append(
            f"{metrics['unassessable_fraction']:.1%} of cells are unassessable "
            "(not screened against the common panel, or below the depth floor)"
        )
    return metrics, failures


def _library_habitat(t: Tables) -> pd.Series:
    lib, sites = t["libraries"], t["sites"]
    if not len(lib) or not len(sites):
        return pd.Series(dtype=object)
    lookup = sites.set_index("site_id")["habitat"]
    return lib["site_id"].astype(str).map(lookup).fillna("unknown")


def _gate_g4(t: Tables, th: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    habitat = _library_habitat(t)
    n = len(habitat)
    known = int((habitat != "unknown").sum())
    per_habitat = habitat.value_counts().to_dict()
    sites = t["sites"]
    sites_per_habitat = sites["habitat"].value_counts().to_dict() if len(sites) else {}
    metrics = {
        "n_libraries": n,
        "libraries_with_known_habitat": known,
        "known_habitat_fraction": _frac(known, n),
        "libraries_per_habitat": {str(k): int(v) for k, v in per_habitat.items()},
        "sites_per_habitat": {str(k): int(v) for k, v in sites_per_habitat.items()},
        "geo_resolution_values": (
            sorted(set(sites["geo_resolution"].dropna().astype(str))) if len(sites) else []
        ),
    }
    failures = []
    if metrics["known_habitat_fraction"] < th.get("min_known_habitat_library_fraction", 0.7):
        failures.append(
            f"only {metrics['known_habitat_fraction']:.1%} of libraries have a classified habitat"
        )
    for level in ("settlement", "wild"):
        if per_habitat.get(level, 0) < th.get("min_libraries_per_habitat", 15):
            failures.append(
                f"{per_habitat.get(level, 0)} libraries in habitat '{level}' "
                f"(minimum {th.get('min_libraries_per_habitat', 15)})"
            )
        if sites_per_habitat.get(level, 0) < th.get("min_sites_per_habitat", 3):
            failures.append(
                f"{sites_per_habitat.get(level, 0)} sites in habitat '{level}'; "
                "site and habitat would be collinear"
            )
    return metrics, failures


def _gate_g5(t: Tables, th: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    lib, hosts = t["libraries"], t["hosts"]
    habitat = _library_habitat(t)
    if not len(lib) or not len(hosts):
        return {"n_overlapping_species": 0}, ["no libraries or hosts to assess overlap"]
    species = lib["host_id"].astype(str).map(
        hosts.set_index("host_id")["host_species_resolved"]
    ).fillna("unknown")
    frame = pd.DataFrame({"species": species.to_numpy(), "habitat": habitat.to_numpy()})
    frame = frame[(frame["habitat"] != "unknown") & (~frame["species"].isin(["unknown", "mixed"]))]
    table = pd.crosstab(frame["species"], frame["habitat"]) if len(frame) else pd.DataFrame()

    min_per_cell = th.get("min_libraries_per_overlapping_species_per_habitat", 2)
    overlapping: list[str] = []
    if {"settlement", "wild"}.issubset(set(table.columns)):
        mask = (table["settlement"] >= min_per_cell) & (table["wild"] >= min_per_cell)
        overlapping = sorted(table.index[mask].astype(str))
    n_overlap_libraries = int(frame["species"].isin(overlapping).sum())

    metrics = {
        "overlap_table": table.astype(int).to_dict() if len(table) else {},
        "overlapping_species": overlapping,
        "n_overlapping_species": len(overlapping),
        "n_libraries_in_overlapping_species": n_overlap_libraries,
        "overlap_library_fraction": _frac(n_overlap_libraries, len(frame)),
    }
    failures = []
    if len(overlapping) < th.get("min_overlapping_host_species", 2):
        failures.append(
            f"only {len(overlapping)} host species are sampled in both habitats with "
            f">={min_per_cell} libraries each; the habitat contrast would rest on "
            "between-species differences rather than a within-species comparison"
        )
    if metrics["overlap_library_fraction"] < th.get("min_overlap_library_fraction", 0.3):
        failures.append(
            f"only {metrics['overlap_library_fraction']:.1%} of habitat-classified libraries "
            "belong to overlap-supported host species"
        )
    return metrics, failures


def _gate_g6(t: Tables, th: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    lib = t["libraries"]
    n = len(lib)
    per_host = lib.groupby("host_id").size() if n else pd.Series(dtype=int)
    repeated = int((per_host > 1).sum())
    has_batch = bool(n and lib["sequencing_batch"].notna().any())
    metrics = {
        "n_libraries": n,
        "n_source_pools": int(lib["host_id"].nunique()) if n else 0,
        "pools_with_multiple_libraries": repeated,
        "max_libraries_per_pool": int(per_host.max()) if len(per_host) else 0,
        "sequencing_batch_available": has_batch,
        "n_batches": int(lib["sequencing_batch"].nunique()) if has_batch else 0,
    }
    failures = []
    if repeated and metrics["max_libraries_per_pool"] > th.get(
        "max_libraries_per_host_without_grouping", 1
    ):
        failures.append(
            f"{repeated} source pools contributed more than one library; the model must "
            "include a pool-level grouping term or these libraries are not independent"
        )
    if th.get("require_batch_variable", True) and not has_batch:
        failures.append("no sequencing batch variable is available to absorb run-level effects")
    return metrics, failures


def _gate_g7(t: Tables, th: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    vir = t["viruses"]
    n = len(vir)
    eligible = int(vir["eligible"].fillna(False).sum()) if n else 0
    unclassified = (
        int(vir["virus_family"].astype(str).str.lower().isin(["unknown", "unclassified", "nan"]).sum())
        if n else 0
    )
    det = t["detections"]
    detected_eligible = 0
    if n and len(det):
        elig_ids = set(vir.loc[vir["eligible"].fillna(False), "virus_id"].astype(str))
        pos = det[det["detection"] == "positive"]
        detected_eligible = int(pos.loc[pos["virus_id"].astype(str).isin(elig_ids), "virus_id"].nunique())
    metrics = {
        "n_lineages": n,
        "n_eligible": eligible,
        "n_eligible_detected": detected_eligible,
        "unclassified_fraction": _frac(unclassified, n),
        "excluded_reasons": (
            vir.loc[~vir["eligible"].fillna(False), "eligibility_reason"].value_counts().to_dict()
            if n else {}
        ),
    }
    failures = []
    if eligible < th.get("min_eligible_lineages", 5):
        failures.append(f"only {eligible} eligible mammal-associated lineages after curation")
    if metrics["unclassified_fraction"] > th.get("max_unclassified_fraction", 0.5):
        failures.append(
            f"{metrics['unclassified_fraction']:.1%} of lineages have no family assignment"
        )
    return metrics, failures


_GATE_FUNCS: dict[str, Callable[[Tables, dict[str, Any]], tuple[dict[str, Any], list[str]]]] = {
    "G1": _gate_g1,
    "G2": _gate_g2,
    "G3": _gate_g3,
    "G4": _gate_g4,
    "G5": _gate_g5,
    "G6": _gate_g6,
    "G7": _gate_g7,
}


def evaluate_gates(tables: Tables, config: dict[str, Any]) -> GateReport:
    """Run every configured gate against the harmonized tables."""
    results: list[GateResult] = []
    for gate_id, spec in (config.get("gates") or {}).items():
        func = _GATE_FUNCS.get(gate_id)
        if func is None:
            continue
        thresholds = spec.get("thresholds", {}) or {}
        try:
            metrics, failures = func(tables, thresholds)
        except Exception as exc:  # a gate must never crash the report
            metrics, failures = {}, [f"gate evaluation failed: {type(exc).__name__}: {exc}"]
        results.append(
            GateResult(
                gate_id=gate_id,
                name=spec.get("name", gate_id),
                condition=spec.get("condition", ""),
                critical=bool(spec.get("critical", False)),
                passed=not failures,
                metrics=metrics,
                failures=failures,
                notes=[f"thresholds: {thresholds}"] if thresholds else [],
            )
        )
    results.sort(key=lambda r: r.gate_id)
    return GateReport(results=results, on_critical_failure=config.get("on_critical_failure", {}))
