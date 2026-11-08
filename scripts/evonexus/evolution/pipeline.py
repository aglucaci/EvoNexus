"""Driver for the evolutionary subworkflow.

Order: select sequences, align, check recombination, infer a tree, test for
selection, relate the supported results to connectivity. Every stage can end
in a documented exclusion, and an excluded lineage produces no selection
numbers at all rather than numbers with a caveat attached.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd

from evonexus.core.config import load_config
from evonexus.core.io import write_json, write_table
from evonexus.core.provenance import RunRecord
from evonexus.evolution import tools
from evonexus.evolution.hostphylo import covariance_from_file
from evonexus.evolution.sequences import (
    gap_fractions,
    mean_pairwise_divergence,
    read_fasta,
    select_lineages,
    write_fasta,
)


def _alignment_quality(records: Dict[str, str], cfg: Dict[str, Any]) -> Dict[str, Any]:
    per_seq, per_col = gap_fractions(records)
    columns = len(next(iter(records.values()))) if records else 0
    max_seq_gap = float(cfg.get("max_gap_fraction_per_sequence", 0.5))
    max_col_gap = float(cfg.get("max_gap_fraction_per_column", 0.5))
    bad_seqs = [n for n, f in per_seq.items() if f > max_seq_gap]
    good_columns = int((per_col <= max_col_gap).sum()) if per_col.size else 0
    return {
        "n_sequences": len(records),
        "n_columns": columns,
        "n_columns_below_gap_threshold": good_columns,
        "sequences_over_gap_threshold": bad_seqs,
        "mean_pairwise_divergence": mean_pairwise_divergence(records),
    }


def _selection_eligible(quality: Dict[str, Any], cfg: Dict[str, Any],
                        recombination: Dict[str, Any]) -> tuple[bool, List[str]]:
    rules = cfg.get("eligibility", {}) or {}
    reasons: List[str] = []
    if quality["n_sequences"] < int(rules.get("min_sequences", 10)):
        reasons.append(
            f"{quality['n_sequences']} sequences (minimum {rules.get('min_sequences', 10)})"
        )
    if quality["n_columns_below_gap_threshold"] < int(rules.get("min_codons", 100)) * 3:
        reasons.append(
            f"{quality['n_columns_below_gap_threshold']} well-covered alignment columns "
            f"(minimum {int(rules.get('min_codons', 100)) * 3})"
        )
    div = quality.get("mean_pairwise_divergence")
    if div is None or not np.isfinite(div):
        reasons.append("mean pairwise divergence could not be computed")
    else:
        lo = float(rules.get("min_mean_pairwise_divergence", 0.01))
        hi = float(rules.get("max_mean_pairwise_divergence", 0.60))
        if div < lo:
            reasons.append(f"mean pairwise divergence {div:.4f} below {lo}; too little signal")
        if div > hi:
            reasons.append(f"mean pairwise divergence {div:.4f} above {hi}; likely saturated")
    if rules.get("require_no_recombination", True):
        if recombination.get("status") == "detected":
            reasons.append("recombination detected; a single-tree selection test is not valid")
        elif recombination.get("status") in {"tool_missing", "not_run"}:
            reasons.append(
                "recombination could not be assessed "
                f"({recombination.get('status')}); selection tests withheld"
            )
    return (not reasons), reasons


def run_evolution(
    processed_dir: str | Path,
    sequences_path: str | Path | None,
    out_dir: str | Path,
    config_name: str = "evolution",
    host_newick: str | Path | None = None,
    connectivity_path: str | Path | None = None,
    synthetic: bool = False,
) -> Dict[str, Any]:
    from evonexus.qc.run import load_canonical

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    record = RunRecord(step="evolution", params={"synthetic": synthetic,
                                                 "sequences": str(sequences_path)})
    record.synthetic = synthetic
    cfg = load_config(config_name)
    tables = load_canonical(processed_dir)

    # ---- host phylogeny -------------------------------------------------
    host_cov_path = None
    host_cfg = cfg.get("host_phylogeny", {}) or {}
    newick = host_newick or host_cfg.get("newick_path")
    if newick and Path(newick).exists():
        species = sorted(
            set(tables["hosts"]["host_species_resolved"].astype(str))
            - {"unknown", "mixed"}
        ) if len(tables["hosts"]) else None
        try:
            cov = covariance_from_file(
                newick, species, bool(host_cfg.get("scale_to_correlation", True))
            )
            host_cov_path = write_json(cov, out / "host_covariance.json")
            record.note(
                f"host phylogeny: {cov['n_species']} species matched, "
                f"{len(cov['missing_species'])} missing"
            )
        except Exception as exc:
            record.warn(f"host phylogeny unusable: {type(exc).__name__}: {exc}")
    else:
        record.note("no host phylogeny supplied; model M5 will not be fitted")

    # ---- viral sequences -------------------------------------------------
    if not sequences_path or not Path(sequences_path).exists():
        payload = {
            "status": "not_run",
            "reason": "no viral sequence file was supplied",
            "host_covariance": str(host_cov_path) if host_cov_path else None,
            "message": (
                "Figures and tables depending on viral phylogenetics or selection "
                "will be marked as not generated."
            ),
        }
        write_json(payload, out / "evolution_status.json")
        record.warn("evolutionary analyses not run: no sequences supplied")
        record.finish(out / "run_record_evolution.json", status="not_run")
        return payload

    sequences = read_fasta(sequences_path)
    selection = select_lineages(tables["viruses"], tables["detections"], sequences, cfg)
    write_table(selection["table"], out / "lineage_selection")
    write_table(selection["excluded"], out / "lineage_exclusions")
    record.note(
        f"{selection['n_selected']} lineages selected for phylogenetics, "
        f"{selection['n_excluded']} excluded"
    )

    align_cfg = cfg.get("alignment", {}) or {}
    phylo_cfg = cfg.get("phylogeny", {}) or {}
    sel_cfg = cfg.get("selection", {}) or {}
    recomb_cfg = cfg.get("recombination", {}) or {}

    results: List[Dict[str, Any]] = []
    site_rows: List[Dict[str, Any]] = []

    for lineage, records in selection["selected"].items():
        lin_dir = out / "lineages" / lineage
        lin_dir.mkdir(parents=True, exist_ok=True)
        raw_fa = write_fasta(records, lin_dir / f"{lineage}.fasta")
        entry: Dict[str, Any] = {"lineage": lineage, "n_sequences": len(records)}

        aln_path = lin_dir / f"{lineage}.aln.fasta"
        aln = tools.run_mafft(raw_fa, aln_path, align_cfg.get("args", ["--auto"]))
        entry["alignment_status"] = aln["status"]
        if aln["status"] != "ok":
            entry["excluded_reason"] = f"alignment not produced ({aln['status']})"
            entry["remedy"] = aln.get("remedy")
            results.append(entry)
            continue

        aligned = read_fasta(aln_path)
        quality = _alignment_quality(aligned, align_cfg)
        entry.update({f"alignment_{k}": v for k, v in quality.items()
                      if k != "sequences_over_gap_threshold"})
        if quality["n_columns_below_gap_threshold"] < int(align_cfg.get("min_aligned_columns", 300)):
            entry["excluded_reason"] = (
                f"{quality['n_columns_below_gap_threshold']} well-covered columns "
                f"(minimum {align_cfg.get('min_aligned_columns', 300)})"
            )
            results.append(entry)
            continue

        recombination: Dict[str, Any] = {"status": "not_run"}
        if recomb_cfg.get("enabled", True):
            gard = tools.run_hyphy("gard", aln_path, aln_path, lin_dir / "gard.json")
            if gard["status"] == "tool_missing":
                recombination = {"status": "tool_missing", "tool": "hyphy"}
            elif gard["status"] == "ok":
                parsed = tools.parse_hyphy_sites(lin_dir / "gard.json", "GARD")
                detected = bool((parsed.get("test_result") or {}).get("breakpoints", 0))
                recombination = {"status": "detected" if detected else "none_detected",
                                 "detail": parsed}
            else:
                recombination = {"status": "failed", "detail": gard.get("stderr_tail")}
        entry["recombination_status"] = recombination["status"]

        tree = tools.run_iqtree(
            aln_path, lin_dir / lineage,
            model=phylo_cfg.get("model", "MFP"),
            bootstrap=int(phylo_cfg.get("bootstrap", 1000)),
            extra=phylo_cfg.get("args", []),
        )
        entry["tree_status"] = tree["status"]
        entry["treefile"] = tree.get("treefile")
        if tree["status"] != "ok" or not tree.get("treefile"):
            entry["excluded_reason"] = f"phylogeny not inferred ({tree['status']})"
            entry["remedy"] = tree.get("remedy")
            results.append(entry)
            continue

        eligible, reasons = _selection_eligible(quality, sel_cfg, recombination)
        entry["selection_eligible"] = eligible
        entry["selection_exclusion_reasons"] = "; ".join(reasons)
        if not eligible:
            record.note(f"{lineage}: selection analyses withheld — {'; '.join(reasons)}")
            results.append(entry)
            continue

        for analysis in sel_cfg.get("analyses", []):
            if not analysis.get("enabled", True):
                continue
            method = analysis["name"]
            res = tools.run_hyphy(method, aln_path, tree["treefile"],
                                  lin_dir / f"{method}.json")
            parsed = tools.parse_hyphy_sites(
                lin_dir / f"{method}.json", method,
                float((sel_cfg.get("significance") or {}).get("p_value", 0.05)),
            )
            site_rows.append(
                {
                    "lineage": lineage,
                    "method": method,
                    "question": analysis.get("question"),
                    "status": parsed.get("status", res["status"]),
                    "n_sites": parsed.get("n_sites"),
                    "n_significant": parsed.get("n_significant"),
                    "p_threshold": parsed.get("p_threshold"),
                }
            )
        results.append(entry)

    results_df = pd.DataFrame(results)
    sites_df = pd.DataFrame(site_rows)
    write_table(results_df, out / "evolution_lineage_results")
    write_table(sites_df, out / "selection_results")

    # ---- integration with connectivity ------------------------------------
    integration: Dict[str, Any] = {"status": "not_run"}
    if (cfg.get("integration", {}) or {}).get("relate_to_connectivity") and connectivity_path:
        try:
            pos = tables["detections"]
            breadth = (
                pos[pos["detection"] == "positive"]
                .merge(tables["libraries"][["library_id", "host_id"]], on="library_id")
                .merge(tables["hosts"][["host_id", "host_species_resolved"]], on="host_id")
                .groupby("virus_id")["host_species_resolved"].nunique()
                .rename("host_breadth").reset_index()
            )
            merged = results_df.merge(breadth, left_on="lineage", right_on="virus_id", how="left")
            write_table(merged, out / "evolution_host_breadth")
            integration = {
                "status": "ok",
                "n_lineages": int(len(merged)),
                "note": (cfg.get("integration", {}) or {}).get("note"),
            }
        except Exception as exc:
            integration = {"status": "failed", "error": f"{type(exc).__name__}: {exc}"}

    missing_tools = sorted({
        t for t in ("mafft", "iqtree2", "hyphy") if not tools.tool_available(t)
    })
    payload = {
        "status": "ok" if not results_df.empty else "no_lineages",
        "n_lineages_attempted": int(len(results_df)),
        "n_with_tree": int((results_df.get("tree_status") == "ok").sum()) if len(results_df) else 0,
        "n_selection_eligible": int(results_df.get("selection_eligible", pd.Series(dtype=bool)).fillna(False).sum()) if len(results_df) else 0,
        "missing_tools": missing_tools,
        "host_covariance": str(host_cov_path) if host_cov_path else None,
        "integration": integration,
        "synthetic": synthetic,
    }
    if missing_tools:
        record.warn(
            "external tools not found: " + ", ".join(missing_tools)
            + " — affected analyses were excluded and recorded, not approximated"
        )
    write_json(payload, out / "evolution_status.json")
    record.finish(out / "run_record_evolution.json", status=payload["status"])
    return payload
