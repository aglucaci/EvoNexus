"""Command-line interface. Every Nextflow process calls one of these."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


def _print(result: Any) -> None:
    print(json.dumps(result, indent=2, default=str))


# ---------------------------------------------------------------------------
# subcommands
# ---------------------------------------------------------------------------
def cmd_simulate(args: argparse.Namespace) -> int:
    from evonexus.simulate import simulate

    out = simulate(args.out_dir, n_sites=args.sites, n_libraries=args.libraries,
                   n_viruses=args.viruses, seed=args.seed)
    _print({"status": "ok", "paths": out["paths"],
            "dimensions": out["truth"]["dimensions"]})
    return 0


def cmd_acquire(args: argparse.Namespace) -> int:
    from evonexus.acquire.run import run_acquire

    _print(run_acquire(
        out_dir=args.out_dir, mode=args.mode, offline=args.offline,
        synthetic=args.synthetic, synthetic_seed=args.seed,
        synthetic_libraries=args.libraries, synthetic_viruses=args.viruses,
    ))
    return 0


def cmd_inventory(args: argparse.Namespace) -> int:
    from evonexus.acquire.inventory import inventory_directory, write_inventory

    inv = inventory_directory(args.raw_dir)
    paths = write_inventory(inv, args.out_dir)
    _print({"status": "ok", "n_files": inv["n_files"],
            "expectations": inv["expectations"],
            "paths": {k: str(v) for k, v in paths.items()}})
    return 0


def cmd_harmonize(args: argparse.Namespace) -> int:
    from evonexus.harmonize.mapping import MappingError
    from evonexus.harmonize.run import run_harmonize

    try:
        _print(run_harmonize(
            inventory_path=args.inventory, raw_dir=args.raw_dir, out_dir=args.out_dir,
            overrides=args.overrides, allow_role_guess=args.allow_role_guess,
            require_approval=not args.skip_approval, synthetic=args.synthetic,
        ))
    except MappingError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    return 0


def cmd_qc(args: argparse.Namespace) -> int:
    from evonexus.qc.run import run_qc

    _print(run_qc(processed_dir=args.processed_dir, out_dir=args.out_dir,
                  synthetic=args.synthetic, fail_on_critical=args.fail_on_critical))
    return 0


def cmd_matrix(args: argparse.Namespace) -> int:
    from evonexus.quantify.matrix import run_matrix

    _print(run_matrix(processed_dir=args.processed_dir, out_dir=args.out_dir,
                      resolution=args.resolution, synthetic=args.synthetic))
    return 0


def cmd_bayesian(args: argparse.Namespace) -> int:
    from evonexus.bayesian.run import run_bayesian

    _print(run_bayesian(
        frame_path=args.frame, out_dir=args.out_dir, models=args.models,
        prior_variant=args.prior_variant, engine=args.engine, draws=args.draws,
        tune=args.tune, chains=args.chains, seed=args.seed,
        host_covariance_path=args.host_covariance, truth_path=args.truth,
        synthetic=args.synthetic, run_sensitivity=not args.no_sensitivity,
        skip_primary=args.skip_primary, processed_dir=args.processed_dir,
    ))
    return 0


def cmd_connectivity(args: argparse.Namespace) -> int:
    from evonexus.connectivity.compute import run_connectivity

    _print(run_connectivity(
        posterior_dir=args.posterior_dir, model=args.model, frame_path=args.frame,
        out_dir=args.out_dir, n_permutations=args.permutations, seed=args.seed,
        host_distance_path=args.host_covariance, synthetic=args.synthetic,
    ))
    return 0


def cmd_evolution(args: argparse.Namespace) -> int:
    from evonexus.evolution.pipeline import run_evolution

    _print(run_evolution(
        processed_dir=args.processed_dir, sequences_path=args.sequences,
        out_dir=args.out_dir, host_newick=args.host_newick,
        connectivity_path=args.connectivity, synthetic=args.synthetic,
    ))
    return 0


def cmd_validate(args: argparse.Namespace) -> int:
    from evonexus.validation.crossval import run_validation

    _print(run_validation(
        frame_path=args.frame, out_dir=args.out_dir, posterior_dir=args.posterior_dir,
        model=args.model, group_columns=args.group_columns, n_folds=args.folds,
        engine=args.engine, draws=args.draws, tune=args.tune, chains=args.chains,
        synthetic=args.synthetic,
    ))
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    from evonexus.reporting.run import run_reporting

    _print(run_reporting(results_dir=args.results_dir, out_dir=args.out_dir))
    return 0


def cmd_schema(args: argparse.Namespace) -> int:
    from evonexus.core.schema import schema_as_dict

    schema = schema_as_dict()
    if args.markdown:
        lines = ["# Data dictionary", "",
                 "Generated from `evonexus.core.schema`. Do not edit by hand.", ""]
        for name, spec in schema.items():
            lines += [
                f"## `{name}.parquet`", "",
                f"- Unit of observation: {spec['unit']}",
                f"- Primary key: {', '.join('`' + k + '`' for k in spec['primary_key'])}",
            ]
            if spec["foreign_keys"]:
                lines.append("- Foreign keys: " + ", ".join(
                    f"`{k}` → `{v}`" for k, v in spec["foreign_keys"].items()))
            lines += ["", "| Column | Type | Required | Nullable | Allowed values | Description |",
                      "| --- | --- | --- | --- | --- | --- |"]
            for col in spec["columns"]:
                allowed = ", ".join(f"`{a}`" for a in col["allowed"]) if col["allowed"] else "—"
                lines.append(
                    f"| `{col['name']}` | {col['dtype']} | "
                    f"{'yes' if col['required'] else 'no'} | "
                    f"{'yes' if col['nullable'] else 'no'} | {allowed} | "
                    f"{col['description'] or '—'} |"
                )
            lines.append("")
        Path(args.markdown).parent.mkdir(parents=True, exist_ok=True)
        Path(args.markdown).write_text("\n".join(lines), encoding="utf-8")
        _print({"status": "ok", "written": args.markdown})
    else:
        _print(schema)
    return 0


def cmd_engines(args: argparse.Namespace) -> int:
    from evonexus.bayesian.fit import available_engines
    from evonexus.core.provenance import record_environment

    _print({"engines": available_engines(), "environment": record_environment()})
    return 0


# ---------------------------------------------------------------------------
# parser
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="evonexus",
        description="EvoNexus — reproducible bat virome ecology and evolution pipeline",
    )
    parser.add_argument("--version", action="store_true", help="print the version and exit")
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("simulate", help="generate deterministic synthetic source data")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sites", type=int, default=12)
    p.add_argument("--libraries", type=int, default=120)
    p.add_argument("--viruses", type=int, default=30)
    p.add_argument("--seed", type=int, default=20250101)
    p.set_defaults(func=cmd_simulate)

    p = sub.add_parser("acquire", help="download sources and write the manifest and inventory")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--mode", choices=["processed_only", "raw_reads"], default="processed_only")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--synthetic", action="store_true")
    p.add_argument("--seed", type=int, default=20250101)
    p.add_argument("--libraries", type=int, default=120)
    p.add_argument("--viruses", type=int, default=30)
    p.set_defaults(func=cmd_acquire)

    p = sub.add_parser("inventory", help="inventory a raw-data directory")
    p.add_argument("--raw-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.set_defaults(func=cmd_inventory)

    p = sub.add_parser("harmonize", help="build the canonical tables")
    p.add_argument("--inventory", required=True)
    p.add_argument("--raw-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--overrides", default=None)
    p.add_argument("--allow-role-guess", action="store_true",
                   help="allow table roles to be guessed from the inventory (test use)")
    p.add_argument("--skip-approval", action="store_true",
                   help="do not require human approval of the mapping (test use)")
    p.add_argument("--synthetic", action="store_true")
    p.set_defaults(func=cmd_harmonize)

    p = sub.add_parser("qc", help="evaluate the mandatory data gates")
    p.add_argument("--processed-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--synthetic", action="store_true")
    p.add_argument("--fail-on-critical", action="store_true")
    p.set_defaults(func=cmd_qc)

    p = sub.add_parser("matrix", help="build the analysis frame and detection matrix")
    p.add_argument("--processed-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--resolution", choices=["lineage", "family", "strict"], default="lineage")
    p.add_argument("--synthetic", action="store_true")
    p.set_defaults(func=cmd_matrix)

    p = sub.add_parser("bayesian", help="fit the hierarchical detection models")
    p.add_argument("--frame", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--models", nargs="*", default=None)
    p.add_argument("--prior-variant", default="regularizing")
    p.add_argument("--engine", default="auto", choices=["auto", "pymc_nuts", "polya_gamma_gibbs"])
    p.add_argument("--draws", type=int, default=None)
    p.add_argument("--tune", type=int, default=None)
    p.add_argument("--chains", type=int, default=None)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--host-covariance", default=None)
    p.add_argument("--truth", default=None, help="synthetic ground-truth JSON for recovery")
    p.add_argument("--processed-dir", default=None)
    p.add_argument("--no-sensitivity", action="store_true")
    p.add_argument("--skip-primary", action="store_true",
                   help="descriptive output only (used when a critical QC gate failed)")
    p.add_argument("--synthetic", action="store_true")
    p.set_defaults(func=cmd_bayesian)

    p = sub.add_parser("connectivity", help="posterior cross-host connectivity")
    p.add_argument("--posterior-dir", required=True)
    p.add_argument("--model", default="M3")
    p.add_argument("--frame", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--permutations", type=int, default=500)
    p.add_argument("--seed", type=int, default=20250101)
    p.add_argument("--host-covariance", default=None)
    p.add_argument("--synthetic", action="store_true")
    p.set_defaults(func=cmd_connectivity)

    p = sub.add_parser("evolution", help="alignment, phylogeny and selection analyses")
    p.add_argument("--processed-dir", required=True)
    p.add_argument("--sequences", default=None)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--host-newick", default=None)
    p.add_argument("--connectivity", default=None)
    p.add_argument("--synthetic", action="store_true")
    p.set_defaults(func=cmd_evolution)

    p = sub.add_parser("validate", help="held-out cross-validation and calibration")
    p.add_argument("--frame", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--posterior-dir", default=None)
    p.add_argument("--model", default="M3")
    p.add_argument("--group-columns", nargs="*", default=["site_id"])
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--engine", default="auto")
    p.add_argument("--draws", type=int, default=None)
    p.add_argument("--tune", type=int, default=None)
    p.add_argument("--chains", type=int, default=None)
    p.add_argument("--synthetic", action="store_true")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("report", help="generate figures, tables and legends")
    p.add_argument("--results-dir", required=True)
    p.add_argument("--out-dir", default=None)
    p.set_defaults(func=cmd_report)

    p = sub.add_parser("schema", help="print or render the canonical data dictionary")
    p.add_argument("--markdown", default=None, help="write Markdown to this path")
    p.set_defaults(func=cmd_schema)

    p = sub.add_parser("engines", help="report the available samplers and environment")
    p.set_defaults(func=cmd_engines)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "version", False):
        import evonexus

        print(evonexus.__version__)
        return 0
    if not getattr(args, "func", None):
        parser.print_help()
        return 1
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
