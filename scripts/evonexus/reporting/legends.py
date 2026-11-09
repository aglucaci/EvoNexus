"""Figure legends generated from validated analysis outputs.

A legend has five parts: the question, the data and the real number of
independent observations, the panels, the statistics, and the caveats. Every
number comes from an output file. If a figure was not generated, the legend
says so and states nothing further.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from evonexus.core.io import write_text
from evonexus.reporting.artifacts import ArtifactStore


def _observations(store: ArtifactStore, figure: str) -> list[str]:
    """The honest denominator for each figure."""
    lines: list[str] = []
    meta = store.matrix_metadata() or {}
    summary = meta.get("summary", {})
    libs = store.canonical("libraries")
    qc = store.qc() or {}

    if figure == "Figure1":
        if libs is not None:
            lines.append(f"{len(libs)} sequencing libraries in the harmonized dataset.")
        wf = store.waterfall()
        if wf is not None and len(wf):
            last = wf.iloc[-1]
            lines.append(
                f"{int(last['cells']):,} library × lineage cells remain after eligibility "
                f"filtering, from {int(last['libraries'])} libraries."
            )
    elif figure in {"Figure2", "Figure5"}:
        if summary:
            lines.append(
                f"{summary.get('n_cells', 0):,} assessed library × lineage cells from "
                f"{summary.get('n_libraries', 0)} libraries, "
                f"{summary.get('n_virus_units', 0)} viral lineage units, "
                f"{summary.get('n_host_species', 0)} host species and "
                f"{summary.get('n_sites', 0)} sites."
            )
            by_habitat = summary.get("libraries_by_habitat") or {}
            if by_habitat:
                lines.append(
                    "Libraries by habitat: "
                    + ", ".join(f"{k} {v}" for k, v in sorted(by_habitat.items()))
                    + "."
                )
        ref = (store.contrast_metadata() or {}).get("reference") or {}
        if ref:
            lines.append(
                f"The contrast is standardized over {ref.get('n_library_profiles', 0)} library "
                f"profiles and {ref.get('n_lineages', 0)} lineages"
                + (" restricted to host species observed in both habitats."
                   if ref.get("overlap_only") else ".")
            )
    elif figure == "Figure3":
        conn_meta = store.connectivity_metadata() or {}
        matrices = conn_meta.get("matrices") or {}
        if matrices:
            first = next(iter(matrices.values()))
            lines.append(
                f"{len(first.get('hosts', []))} host species, "
                f"{summary.get('n_virus_units', 0)} viral lineage units, "
                f"habitats: {', '.join(conn_meta.get('habitats', []))}."
            )
        perm = conn_meta.get("permutation_null") or {}
        if perm.get("n_permutations"):
            lines.append(
                f"The permutation null uses {perm['n_permutations']} shuffles of host labels "
                f"within {', '.join(perm.get('strata', []))}."
            )
    elif figure == "Figure4":
        status = store.evolution_status() or {}
        if status.get("n_lineages_attempted") is not None:
            lines.append(
                f"{status.get('n_lineages_attempted', 0)} lineages entered the evolutionary "
                f"workflow; {status.get('n_with_tree', 0)} yielded a phylogeny and "
                f"{status.get('n_selection_eligible', 0)} met the selection-analysis gates."
            )

    if qc and not qc.get("primary_model_allowed", True):
        lines.append(
            "A critical QC gate failed, so the primary habitat model was not fitted."
        )
    return lines


def _statistics(store: ArtifactStore, figure: str) -> list[str]:
    lines: list[str] = []
    status = store.bayesian_status() or {}
    engine = status.get("engine")
    engine_text = {
        "pymc_nuts": "PyMC with NUTS",
        "polya_gamma_gibbs": (
            "the built-in Polya-Gamma Gibbs sampler (PyMC was not installed; "
            "see REPRODUCIBILITY.md)"
        ),
    }.get(engine, engine)

    if figure in {"Figure2", "Figure5"} and status:
        sampler = status.get("sampler", {})
        lines.append(
            f"Model {status.get('primary_model')} fitted with {engine_text}: "
            f"{sampler.get('chains')} chains, {sampler.get('draws')} post-warmup draws each, "
            f"{sampler.get('tune')} warmup iterations, seed {sampler.get('seed')}."
        )
        conv = store.convergence(status.get("primary_model", "")) or {}
        metrics = conv.get("metrics", {})
        if metrics:
            lines.append(
                f"Convergence: maximum R-hat {metrics.get('max_rhat'):.3f}, "
                f"minimum bulk ESS {metrics.get('min_ess_bulk'):.0f}, "
                f"{metrics.get('divergences', 0)} divergent transitions."
                + ("" if conv.get("passed") else " Convergence gates were NOT met.")
            )
        lines.append("Intervals are 95% highest-density intervals of the posterior.")

    if figure == "Figure2":
        contrast = store.contrast()
        if contrast is not None and len(contrast):
            row = contrast[contrast["quantity"].str.contains("standardized difference")]
            if len(row):
                r = row.iloc[0]
                lines.append(
                    f"Standardized difference (settlement − wild): {r['mean']:.4f} "
                    f"(95% HDI {r['hdi_95_low']:.4f} to {r['hdi_95_high']:.4f}); "
                    f"posterior probability of a positive difference "
                    f"{r['p_direction_positive']:.2f}."
                )
    if figure == "Figure3":
        lines.append(
            "C(i,j) is the posterior expected number of viral lineages that two "
            "hypothetical libraries, one from each host species, would both be "
            "positive for under a common sampling design."
        )
    if figure == "Figure5":
        val = (store.validation_metadata() or {}).get("cross_validation") or {}
        for group, res in val.items():
            overall = (res or {}).get("overall") or {}
            if overall:
                lines.append(
                    f"Held-out {group}: mean log predictive density "
                    f"{overall.get('mean_lpd'):.4f}, Brier score {overall.get('brier'):.4f}, "
                    f"AUC {overall.get('auc'):.3f} over {overall.get('n'):,} held-out cells."
                )
        lines.append(
            "Random effects of held-out sites and batches were integrated over their "
            "estimated population distributions, not estimated from held-out data."
        )
    return lines


def build_legends(store: ArtifactStore, config: dict[str, Any],
                  figure_results: dict[str, Any], out_dir: str | Path) -> dict[str, Any]:
    out = Path(out_dir) / "figure_legends"
    template = config.get("legend_template", {}) or {}
    boilerplate = template.get("caveat_boilerplate", []) or []
    synthetic = store.is_synthetic()
    written: dict[str, Any] = {}

    for name, spec in (config.get("figures") or {}).items():
        result = figure_results.get(name, {})
        status = result.get("status", "not_generated")
        lines: list[str] = [f"# {name} legend", ""]

        if synthetic:
            lines += ["> **SYNTHETIC DATA — software validation only. "
                      "No value below describes a real bat, virus or place.**", ""]

        if status != "ok":
            lines += [
                f"**{name}: {spec.get('title', name)} — NOT GENERATED.**",
                "",
                "Required analysis outputs were unavailable:",
                "",
            ]
            lines += [f"- {m}" for m in (result.get("missing") or [result.get("error", "unknown")])]
            lines += [
                "",
                "No sample size, credible interval or conclusion is reported for this "
                "figure because none was computed.",
                "",
            ]
            written[name] = {"status": "not_generated",
                             "path": str(write_text("\n".join(lines), out / f"{name}_legend.md"))}
            continue

        lines += [f"**{spec.get('title', name)}.**", ""]
        lines += ["## Question", "", spec.get("question", "—"), ""]

        observations = _observations(store, name)
        lines += ["## Data and number of observations", ""]
        lines += observations or ["Observation counts were not available."]
        lines.append("")

        lines += ["## Panels", ""]
        for panel in spec.get("panels", []):
            lines.append(f"**({panel['id']})** {panel['name']}.")
        lines.append("")

        statistics = _statistics(store, name)
        lines += ["## Statistics", ""]
        lines += statistics or ["No model statistics apply to this figure."]
        lines.append("")

        lines += ["## Caveats and interpretation limits", ""]
        lines += [f"- {c}" for c in boilerplate]
        qc = store.qc() or {}
        failed = [g for g in qc.get("gates", []) if g.get("status") != "pass"]
        if failed:
            lines.append(
                "- Unmet QC conditions at the time of generation: "
                + "; ".join(f"{g['gate_id']} ({g['status']})" for g in failed)
                + "."
            )
        if synthetic:
            lines.append("- These data are synthetic; the figure validates software only.")
        lines.append("")

        written[name] = {"status": "ok",
                         "path": str(write_text("\n".join(lines), out / f"{name}_legend.md"))}
    return written
