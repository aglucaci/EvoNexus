"""Driver for the QC gate step, including the human-readable report."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from evonexus.core.config import load_config
from evonexus.core.io import read_table, write_json, write_table, write_text
from evonexus.core.provenance import RunRecord
from evonexus.qc.gates import GateReport, evaluate_gates

TABLE_NAMES = ("sites", "hosts", "libraries", "viruses", "detections", "ecology", "crosswalk")


def load_canonical(processed_dir: str | Path) -> dict[str, pd.DataFrame]:
    base = Path(processed_dir)
    tables: dict[str, pd.DataFrame] = {}
    for name in TABLE_NAMES:
        try:
            tables[name] = read_table(base / name)
        except FileNotFoundError:
            tables[name] = pd.DataFrame()
    return tables


def _exclusion_waterfall(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """How many observations survive each eligibility rule, in order."""
    lib, hosts, sites, vir, det = (
        tables["libraries"], tables["hosts"], tables["sites"], tables["viruses"], tables["detections"]
    )
    rows = []

    def add(stage: str, libraries: int, lineages: int, cells: int, note: str) -> None:
        rows.append({"stage": stage, "libraries": libraries, "lineages": lineages,
                     "cells": cells, "note": note})

    add("Source records", len(lib), len(vir), len(det), "After harmonization")

    if not len(lib) or not len(det):
        return pd.DataFrame(rows)

    habitat = lib["site_id"].astype(str).map(
        sites.set_index("site_id")["habitat"] if len(sites) else pd.Series(dtype=object)
    ).fillna("unknown")
    keep_lib = set(lib.loc[habitat != "unknown", "library_id"].astype(str))
    sub = det[det["library_id"].astype(str).isin(keep_lib)]
    add("Known habitat", len(keep_lib), vir["virus_id"].nunique(), len(sub),
        "Libraries whose site habitat is settlement or wild")

    species = lib["host_id"].astype(str).map(
        hosts.set_index("host_id")["host_species_resolved"] if len(hosts) else pd.Series(dtype=object)
    ).fillna("unknown")
    keep_lib &= set(lib.loc[~species.isin(["unknown", "mixed"]), "library_id"].astype(str))
    sub = det[det["library_id"].astype(str).isin(keep_lib)]
    add("Single documented host species", len(keep_lib), vir["virus_id"].nunique(), len(sub),
        "Mixed-species and unidentified pools excluded from the single-host model")

    elig = set(vir.loc[vir["eligible"].fillna(False), "virus_id"].astype(str)) if len(vir) else set()
    sub = sub[sub["virus_id"].astype(str).isin(elig)]
    add("Eligible viral lineages", len(keep_lib), len(elig), len(sub),
        "Mammal-associated RNA viral families only")

    assessable = sub[sub["detection"].isin(["positive", "negative"])]
    add("Assessable cells", assessable["library_id"].nunique(), assessable["virus_id"].nunique(),
        len(assessable), "Screened against the common panel at adequate depth")

    return pd.DataFrame(rows)


def _render_report(report: GateReport, tables: dict[str, pd.DataFrame],
                   waterfall: pd.DataFrame, synthetic: bool) -> str:
    lines = ["# EvoNexus QC report", ""]
    if synthetic:
        lines += ["> **SYNTHETIC DATA — software validation only. Not a scientific result.**", ""]

    verdict = "PASS" if report.primary_model_allowed else "STOP"
    lines += [
        f"## Verdict: **{verdict}**",
        "",
        (
            "All critical gates passed. The primary Bayesian habitat analysis may run."
            if report.primary_model_allowed
            else report.on_critical_failure.get(
                "message",
                "A critical gate failed; the primary habitat analysis will not run.",
            )
        ),
        "",
        "## Gate summary",
        "",
        "| Gate | Name | Critical | Status |",
        "| --- | --- | --- | --- |",
    ]
    for r in report.results:
        lines.append(
            f"| {r.gate_id} | {r.name} | {'yes' if r.critical else 'no'} | "
            f"**{r.status.upper()}** |"
        )
    lines.append("")

    for r in report.results:
        lines += [f"### {r.gate_id} — {r.name}", "", f"*{r.condition}*", ""]
        lines.append(f"Status: **{r.status.upper()}**")
        lines.append("")
        if r.failures:
            lines += ["Unmet conditions:", ""] + [f"- {f}" for f in r.failures] + [""]
        if r.metrics:
            lines += ["| Metric | Value |", "| --- | --- |"]
            for k, v in r.metrics.items():
                if isinstance(v, float):
                    v = f"{v:.4f}"
                elif isinstance(v, (dict, list)):
                    v = str(v)[:240]
                lines.append(f"| `{k}` | {v} |")
            lines.append("")

    if len(waterfall):
        lines += ["## Exclusion waterfall", "",
                  "| Stage | Libraries | Lineages | Cells | Note |",
                  "| --- | ---: | ---: | ---: | --- |"]
        for _, row in waterfall.iterrows():
            lines.append(
                f"| {row['stage']} | {row['libraries']} | {row['lineages']} | "
                f"{row['cells']} | {row['note']} |"
            )
        lines.append("")

    lines += [
        "## Interpretation limits carried forward",
        "",
        "- Detections are library-level viral RNA detections. They are not confirmed",
        "  infections of individual animals, and pooled libraries cannot be read as",
        "  one independently sampled animal per constituent species.",
        "- Cells marked `unassessable` are not negatives. They carry no response and",
        "  are excluded from the likelihood rather than imputed as zeros.",
        "- Any habitat difference estimated downstream is an observational association",
        "  standardized to a common sampling distribution, not a causal effect of",
        "  urbanization.",
        "",
    ]
    return "\n".join(lines)


def run_qc(
    processed_dir: str | Path,
    out_dir: str | Path,
    gates_config: str = "qc_gates",
    synthetic: bool = False,
    fail_on_critical: bool = False,
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    record = RunRecord(step="qc", params={"gates_config": gates_config, "synthetic": synthetic})
    record.synthetic = synthetic

    tables = load_canonical(processed_dir)
    for name, df in tables.items():
        if len(df):
            record.add_input(Path(processed_dir) / name, role=f"canonical:{name}", checksum=False)

    config = load_config(gates_config)
    report = evaluate_gates(tables, config)
    waterfall = _exclusion_waterfall(tables)

    payload = report.to_dict()
    payload["synthetic"] = synthetic
    json_path = write_json(payload, out / "qc_gates.json")
    md_path = write_text(_render_report(report, tables, waterfall, synthetic), out / "qc_report.md")
    # A one-word verdict file. The workflow reads this to route to the primary
    # model or the descriptive branch; a plain file works on every Nextflow
    # version, whereas an `env` output does not.
    verdict_path = write_text(
        f"{report.primary_model_allowed}\n", out / "primary_model_allowed.txt"
    )
    wf_path = write_table(waterfall, out / "exclusion_waterfall")
    gate_table = write_table(
        pd.DataFrame(
            [
                {"gate_id": r.gate_id, "name": r.name, "critical": r.critical,
                 "status": r.status, "failures": "; ".join(r.failures)}
                for r in report.results
            ]
        ),
        out / "qc_gate_summary",
    )
    for p in (json_path, md_path, wf_path, gate_table, verdict_path):
        record.add_output(p, role="qc")

    for r in report.critical_failures:
        record.error(f"{r.gate_id} ({r.name}) failed: {'; '.join(r.failures)}")
    for r in report.warnings:
        record.warn(f"{r.gate_id} ({r.name}): {'; '.join(r.failures)}")

    status = "ok" if report.primary_model_allowed else "gate_failure"
    record.finish(out / "run_record_qc.json", status=status)

    if fail_on_critical and not report.primary_model_allowed:
        raise SystemExit(
            "Critical QC gate failure; the primary habitat analysis was not run. "
            f"See {md_path}"
        )

    return {
        "qc_json": str(json_path),
        "qc_report": str(md_path),
        "waterfall": str(wf_path),
        "primary_model_allowed": report.primary_model_allowed,
        "status": status,
    }
