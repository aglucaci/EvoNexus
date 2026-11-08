"""Driver for the HARMONIZE step."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

from evonexus.core.config import load_config
from evonexus.core.io import read_json, write_json, write_table
from evonexus.core.provenance import RunRecord
from evonexus.core.schema import check_foreign_keys, validate_table
from evonexus.harmonize.canonical import build_canonical
from evonexus.harmonize.mapping import MappingError, resolve_mapping, write_mapping_report


def run_harmonize(
    inventory_path: str | Path,
    raw_dir: str | Path,
    out_dir: str | Path,
    mapping_config: str = "metadata_mapping",
    filters_config: str = "viral_filters",
    overrides: str | Path | None = None,
    allow_role_guess: bool = False,
    require_approval: bool = True,
    synthetic: bool = False,
) -> Dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    record = RunRecord(
        step="harmonize",
        params={
            "mapping_config": mapping_config,
            "filters_config": filters_config,
            "allow_role_guess": allow_role_guess,
            "require_approval": require_approval,
            "synthetic": synthetic,
        },
    )
    record.synthetic = synthetic
    record.add_input(inventory_path, role="inventory")

    inventory = read_json(inventory_path)
    mapping_cfg = load_config(mapping_config)
    filters = load_config(filters_config)

    if overrides is None:
        candidate = Path(raw_dir) / "synthetic" / "mapping_overrides.json"
        overrides = candidate if candidate.exists() else None

    mapping = resolve_mapping(
        inventory,
        mapping_cfg,
        overrides_path=overrides,
        allow_role_guess=allow_role_guess,
        require_approval=require_approval,
    )
    report_path = write_mapping_report(mapping, out / "mapping_report.md")
    record.add_output(report_path, role="mapping_report")

    problems = mapping.blocking_problems
    if require_approval and not mapping.approved:
        problems.insert(0, "metadata mapping is not approved by a human reviewer")

    if problems:
        for p in problems:
            record.error(p)
        write_json(
            {
                "status": "blocked",
                "problems": problems,
                "next_step": (
                    "Read results/qc/inventory.md and mapping_report.md, complete "
                    "config/metadata_mapping.yaml, then set `approved: true`."
                ),
            },
            out / "harmonize_status.json",
        )
        record.finish(out / "run_record_harmonize.json", status="blocked")
        raise MappingError(
            "harmonization blocked:\n  - " + "\n  - ".join(problems)
        )

    built = build_canonical(raw_dir, mapping, filters)
    tables = built["tables"]

    schema_problems: Dict[str, Any] = {}
    written: Dict[str, str] = {}
    for name, df in tables.items():
        issues = validate_table(df, name, strict=False)
        if issues:
            schema_problems[name] = issues
            for issue in issues:
                record.warn(issue)
        path = write_table(df, out / name)
        written[name] = str(path)
        record.add_output(path, role=f"canonical:{name}")

    fk_problems = check_foreign_keys(tables, strict=False)
    for p in fk_problems:
        record.warn(f"foreign key: {p}")

    diagnostics = built["diagnostics"]
    diagnostics["schema_problems"] = schema_problems
    diagnostics["foreign_key_problems"] = fk_problems
    diagnostics["mapping_approved"] = mapping.approved
    diagnostics["synthetic"] = synthetic
    diag_path = write_json(diagnostics, out / "harmonize_diagnostics.json")
    record.add_output(diag_path, role="diagnostics")

    hard_failures = {k: v for k, v in schema_problems.items() if k in {"detections", "libraries"}}
    status = "ok" if not hard_failures else "warning"
    record.finish(out / "run_record_harmonize.json", status=status)

    return {"tables": written, "diagnostics": str(diag_path),
            "mapping_report": str(report_path), "status": status}
