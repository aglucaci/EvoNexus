"""Driver for manuscript production, plus the overall provenance report."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

from evonexus.core.config import load_config
from evonexus.core.io import write_json, write_text
from evonexus.core.provenance import RunRecord, record_environment
from evonexus.reporting.artifacts import ArtifactStore
from evonexus.reporting.figures import build_figures
from evonexus.reporting.legends import build_legends
from evonexus.reporting.tables import build_tables


def _provenance_report(store: ArtifactStore, figures: Dict[str, Any],
                       tables: Dict[str, Any]) -> str:
    synthetic = store.is_synthetic()
    status = store.bayesian_status() or {}
    qc = store.qc() or {}
    lines = ["# EvoNexus analysis report", ""]
    if synthetic:
        lines += ["> **SYNTHETIC DATA — software validation only. "
                  "Nothing below is a scientific result.**", ""]

    lines += ["## Provenance", ""]
    env = record_environment()
    lines += [
        f"- Generated: {env['timestamp_utc']}",
        f"- Git commit: `{env.get('git_commit') or 'unknown'}`"
        + (" (working tree dirty)" if env.get("git_dirty") else ""),
        f"- Python {env['python']} on {env['platform']}",
        f"- Sampler engine: {status.get('engine', 'not run')}",
        "",
    ]
    packages = env.get("packages", {})
    lines += ["| Package | Version |", "| --- | --- |"]
    lines += [f"| {k} | {v} |" for k, v in sorted(packages.items())]
    lines.append("")

    lines += ["## Quality gates", ""]
    if qc:
        lines.append(
            f"Primary habitat model allowed: **{qc.get('primary_model_allowed')}** "
            f"({qc.get('n_critical_failures', 0)} critical failures, "
            f"{qc.get('n_warnings', 0)} warnings)."
        )
        lines.append("")
        lines += ["| Gate | Name | Status |", "| --- | --- | --- |"]
        for gate in qc.get("gates", []):
            lines.append(f"| {gate['gate_id']} | {gate['name']} | {gate['status'].upper()} |")
    else:
        lines.append("No QC report was found.")
    lines.append("")

    lines += ["## Figures", "", "| Figure | Status | Files |", "| --- | --- | --- |"]
    for name, result in figures.items():
        files = ", ".join(Path(f).name for f in result.get("files", []))
        lines.append(f"| {name} | {result.get('status')} | {files} |")
    lines.append("")

    lines += ["## Tables", "", "| Table | Status | Detail |", "| --- | --- | --- |"]
    for name, result in tables.items():
        detail = (
            f"{result.get('rows')} rows" if result.get("status") == "ok"
            else result.get("reason", "")
        )
        lines.append(f"| {name} | {result.get('status')} | {detail} |")
    lines.append("")

    not_generated = [n for n, r in figures.items() if r.get("status") != "ok"]
    not_generated += [n for n, r in tables.items() if r.get("status") != "ok"]
    if not_generated:
        lines += [
            "## Not generated",
            "",
            "These outputs were not produced because their inputs were unavailable. "
            "No substitute values were used.",
            "",
        ]
        lines += [f"- {n}" for n in not_generated]
        lines.append("")

    lines += [
        "## Interpretation limits",
        "",
        "- Observations are library-level viral RNA detections, not confirmed infections",
        "  of individual animals.",
        "- The habitat contrast is an observational association standardized to a common",
        "  sampling distribution. It is not a causal effect of urbanization.",
        "- Connectivity is expected co-detection under the fitted model, not inferred",
        "  host-to-host transmission.",
        "- The independent individual-bat dataset is a secondary evaluation of whether",
        "  host-sharing summaries remain interpretable without pooled sequencing. It is",
        "  not a replication of the habitat contrast.",
        "",
    ]
    return "\n".join(lines)


def run_reporting(
    results_dir: str | Path,
    out_dir: str | Path | None = None,
    figures_config: str = "figures",
    synthetic: bool | None = None,
) -> Dict[str, Any]:
    results_dir = Path(results_dir)
    out = Path(out_dir) if out_dir else results_dir
    store = ArtifactStore(results_dir)
    config = load_config(figures_config)

    is_synthetic = store.is_synthetic() if synthetic is None else synthetic
    record = RunRecord(step="reporting", params={"results_dir": str(results_dir)})
    record.synthetic = is_synthetic

    figures = build_figures(store, config, out)
    tables = build_tables(store, config, out)
    legends = build_legends(store, config, figures, out)

    report_path = write_text(_provenance_report(store, figures, tables),
                             out / "manuscript" / "analysis_report.md")
    manifest = {
        "synthetic": is_synthetic,
        "figures": {k: {kk: vv for kk, vv in v.items() if kk != "source_tables"}
                    for k, v in figures.items()},
        "tables": tables,
        "legends": legends,
        "missing_artifacts": sorted(set(store.missing)),
        "report": str(report_path),
    }
    manifest_path = write_json(manifest, out / "manuscript" / "manuscript_manifest.json")
    record.add_output(report_path, role="report", checksum=False)
    record.add_output(manifest_path, role="manifest", checksum=False)

    n_missing = sum(1 for v in figures.values() if v.get("status") != "ok")
    if n_missing:
        record.warn(f"{n_missing} figure(s) could not be generated; marked as NOT GENERATED")
    record.note(f"wrote {len(figures)} figures, {len(tables)} tables, {len(legends)} legends")
    record.finish(out / "manuscript" / "run_record_reporting.json", status="ok")

    return {"status": "ok", "report": str(report_path), "manifest": str(manifest_path),
            "figures": len(figures), "tables": len(tables)}
