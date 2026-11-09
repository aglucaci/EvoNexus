"""Figure generation.

Each figure declares its required artifacts. When one is missing, the figure
is still written — as a clearly labelled placeholder naming what is absent —
so a reader can never mistake a gap for a result. Synthetic runs are
watermarked across every panel.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from evonexus.core.io import write_json, write_table
from evonexus.reporting.artifacts import ArtifactStore

MM = 1 / 25.4


def _apply_style(style: dict[str, Any]) -> None:
    plt.rcParams.update(
        {
            "font.size": style.get("font_size_pt", 7),
            "font.family": "sans-serif",
            "font.sans-serif": style.get("font_family", ["DejaVu Sans"]),
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "savefig.bbox": "tight",
            "figure.dpi": 150,
        }
    )


def _new_figure(style: dict[str, Any]) -> plt.Figure:
    return plt.figure(
        figsize=(style.get("width_mm", 180) * MM, style.get("height_mm", 200) * MM)
    )


def _panel_label(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.12, 1.08, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left")


def _watermark(fig: plt.Figure, text: str) -> None:
    fig.text(0.5, 0.5, text, fontsize=22, color="#D62828", alpha=0.16,
             ha="center", va="center", rotation=30, zorder=1000)


def _not_generated(fig: plt.Figure, title: str, missing: list[str]) -> None:
    fig.clf()
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.5, 0.72, title, ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(0.5, 0.58, "NOT GENERATED", ha="center", va="center", fontsize=16,
            color="#B4503C", fontweight="bold")
    ax.text(
        0.5, 0.40,
        "Required analysis outputs are unavailable:\n" + "\n".join(f"• {m}" for m in missing),
        ha="center", va="top", fontsize=8,
    )
    ax.text(0.5, 0.06,
            "No sample size, interval or conclusion is shown because none was computed.",
            ha="center", va="center", fontsize=7, style="italic", color="#555555")


def _save(fig: plt.Figure, out_dir: Path, name: str, style: dict[str, Any]) -> list[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    for fmt in style.get("formats", ["pdf"]):
        path = out_dir / f"{name}.{fmt}"
        fig.savefig(path, format=fmt, dpi=style.get("dpi", 400))
        written.append(str(path))
    plt.close(fig)
    return written


# ---------------------------------------------------------------------------
# Figure 1 — sampling structure
# ---------------------------------------------------------------------------
def figure1(store: ArtifactStore, style: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    sites = store.canonical("sites")
    libs = store.canonical("libraries")
    hosts = store.canonical("hosts")
    waterfall = store.waterfall()
    frame = store.analysis_frame()

    missing = [n for n, v in (("sites", sites), ("libraries", libs), ("hosts", hosts),
                              ("exclusion waterfall", waterfall)) if v is None or not len(v)]
    fig = _new_figure(style)
    if missing:
        _not_generated(fig, "Figure 1: Sampling structure", missing)
        return {"figure": "Figure1", "status": "not_generated", "missing": missing,
                "files": _save(fig, out_dir, "Figure1", style), "source_tables": {}}

    palette = style.get("palette", {})
    gs = GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)
    sources: dict[str, pd.DataFrame] = {}

    # (a) geography
    ax = fig.add_subplot(gs[0, 0])
    _panel_label(ax, "a")
    per_site = libs.groupby("site_id").size().rename("n_libraries").reset_index()
    geo = sites.merge(per_site, on="site_id", how="left").fillna({"n_libraries": 0})
    for habitat, colour in (("settlement", palette.get("settlement", "#B4503C")),
                            ("wild", palette.get("wild", "#2F6B5E")),
                            ("unknown", palette.get("light", "#C9C9C9"))):
        sub = geo[geo["habitat"] == habitat]
        if not len(sub):
            continue
        ax.scatter(sub["longitude"], sub["latitude"], s=6 + 2.2 * sub["n_libraries"],
                   c=colour, alpha=0.8, edgecolors="white", linewidths=0.4, label=habitat)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Sampling sites (point area = libraries)", fontsize=7)
    ax.legend(frameon=False, fontsize=6, loc="best")
    sources["Figure1a_sites"] = geo

    # (b) host x habitat overlap
    ax = fig.add_subplot(gs[0, 1])
    _panel_label(ax, "b")
    merged = libs.merge(hosts[["host_id", "host_species_resolved"]], on="host_id", how="left")
    merged = merged.merge(sites[["site_id", "habitat"]], on="site_id", how="left")
    overlap = pd.crosstab(
        merged["host_species_resolved"].fillna("unknown"), merged["habitat"].fillna("unknown")
    )
    for col in ("settlement", "wild"):
        if col not in overlap.columns:
            overlap[col] = 0
    overlap = overlap.sort_values(["settlement", "wild"], ascending=False).head(14)
    display = overlap[[c for c in ("settlement", "wild", "unknown") if c in overlap.columns]]
    ax.imshow(display.to_numpy(), aspect="auto", cmap="Greys")
    ax.set_xticks(range(display.shape[1]))
    ax.set_xticklabels(display.columns, rotation=30, ha="right")
    ax.set_yticks(range(display.shape[0]))
    ax.set_yticklabels(display.index, fontsize=5.5, style="italic")
    for i in range(display.shape[0]):
        for j in range(display.shape[1]):
            v = int(display.iat[i, j])
            ax.text(j, i, str(v), ha="center", va="center", fontsize=5,
                    color="white" if v > display.to_numpy().max() * 0.55 else "black")
    ax.set_title("Libraries per host species and habitat", fontsize=7)
    sources["Figure1b_overlap"] = overlap.reset_index()

    # (c) sequencing effort
    ax = fig.add_subplot(gs[1, 0])
    _panel_label(ax, "c")
    depth = merged.copy()
    depth["log10_effective_depth"] = np.log10(
        pd.to_numeric(depth["effective_depth"], errors="coerce").clip(lower=1)
    )
    groups, labels, colours = [], [], []
    for habitat, colour in (("settlement", palette.get("settlement", "#B4503C")),
                            ("wild", palette.get("wild", "#2F6B5E"))):
        values = depth.loc[depth["habitat"] == habitat, "log10_effective_depth"].dropna()
        if len(values):
            groups.append(values.to_numpy())
            labels.append(habitat)
            colours.append(colour)
    if groups:
        parts = ax.violinplot(groups, showextrema=False, widths=0.75)
        for body, colour in zip(parts["bodies"], colours):
            body.set_facecolor(colour)
            body.set_alpha(0.55)
        for i, values in enumerate(groups, start=1):
            ax.scatter(np.full(len(values), i) + np.random.default_rng(0).normal(0, 0.035, len(values)),
                       values, s=2, color="black", alpha=0.35, zorder=3)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels)
    ax.set_ylabel("log₁₀ effective depth (reads)")
    ax.set_title("Sequencing effort by habitat", fontsize=7)
    sources["Figure1c_depth"] = depth[["library_id", "habitat", "tissue_group",
                                       "effective_depth", "log10_effective_depth"]]

    # (d) exclusion waterfall
    ax = fig.add_subplot(gs[1, 1])
    _panel_label(ax, "d")
    wf = waterfall.copy()
    ax.barh(range(len(wf)), wf["cells"], color=palette.get("neutral", "#4A4A4A"), alpha=0.85)
    ax.set_yticks(range(len(wf)))
    ax.set_yticklabels(wf["stage"], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("Library × lineage cells")
    for i, (cells, libs_n) in enumerate(zip(wf["cells"], wf["libraries"])):
        ax.text(cells, i, f" {int(cells):,} ({int(libs_n)} lib.)", va="center", fontsize=5.5)
    ax.set_title("From source records to modelled cells", fontsize=7)
    sources["Figure1d_exclusions"] = wf

    n_cells = len(frame) if frame is not None else None
    fig.suptitle("Figure 1 · Sampling structure of the analysed libraries", fontsize=9, y=0.98)
    if store.is_synthetic():
        _watermark(fig, style.get("synthetic_watermark", "SYNTHETIC DATA"))

    return {
        "figure": "Figure1",
        "status": "ok",
        "files": _save(fig, out_dir, "Figure1", style),
        "source_tables": sources,
        "n": {"libraries": len(libs), "sites": len(sites),
              "modelled_cells": n_cells},
    }


# ---------------------------------------------------------------------------
# Figure 2 — habitat inference
# ---------------------------------------------------------------------------
def figure2(store: ArtifactStore, style: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    contrast = store.contrast()
    draws = store.contrast_draws()
    by_family = store.contrast_by_family()
    raw = store.raw_rates()
    sensitivity = store.sensitivity()

    missing = [n for n, v in (("habitat contrast", contrast), ("contrast draws", draws),
                              ("raw detection rates", raw)) if v is None or not len(v)]
    fig = _new_figure(style)
    if missing:
        _not_generated(fig, "Figure 2: Bayesian habitat inference", missing)
        return {"figure": "Figure2", "status": "not_generated", "missing": missing,
                "files": _save(fig, out_dir, "Figure2", style), "source_tables": {}}

    palette = style.get("palette", {})
    gs = GridSpec(2, 2, figure=fig, hspace=0.5, wspace=0.55)
    sources: dict[str, pd.DataFrame] = {}

    # (a) raw rates
    ax = fig.add_subplot(gs[0, 0])
    _panel_label(ax, "a")
    for i, (habitat, colour) in enumerate(
        (("settlement", palette.get("settlement", "#B4503C")),
         ("wild", palette.get("wild", "#2F6B5E")))
    ):
        values = raw.loc[raw["habitat"] == habitat, "detection_rate"].dropna()
        if not len(values):
            continue
        rng = np.random.default_rng(1)
        ax.scatter(np.full(len(values), i) + rng.normal(0, 0.06, len(values)), values,
                   s=4, color=colour, alpha=0.6)
        ax.hlines(values.mean(), i - 0.25, i + 0.25, color="black", linewidth=1.2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["settlement", "wild"])
    ax.set_ylabel("Per-library detection rate")
    ax.set_title("Unadjusted detections", fontsize=7)
    sources["Figure2a_raw"] = raw

    # (b) standardized predictions
    ax = fig.add_subplot(gs[0, 1])
    _panel_label(ax, "b")
    for colour, column, label in (
        (palette.get("settlement", "#B4503C"), "p_settlement", "settlement"),
        (palette.get("wild", "#2F6B5E"), "p_wild", "wild"),
    ):
        if column in draws:
            ax.hist(draws[column], bins=40, color=colour, alpha=0.6, label=label, density=True)
    ax.set_xlabel("Standardized detection probability")
    ax.set_ylabel("Posterior density")
    ax.legend(frameon=False, fontsize=6)
    ax.set_title("Standardized to a common sampling design", fontsize=7)
    sources["Figure2b_standardized"] = draws

    # (c) contrast overall and by family
    ax = fig.add_subplot(gs[1, 0])
    _panel_label(ax, "c")
    rows: list[dict[str, Any]] = []
    overall = contrast[contrast["quantity"].str.contains("standardized difference")]
    if len(overall):
        r = overall.iloc[0]
        rows.append({"label": "Overall", "mean": r["mean"],
                     "low": r["hdi_95_low"], "high": r["hdi_95_high"]})
    if by_family is not None and len(by_family):
        for _, r in by_family.sort_values("mean").iterrows():
            rows.append({"label": r.get("virus_family", r["quantity"]), "mean": r["mean"],
                         "low": r["hdi_95_low"], "high": r["hdi_95_high"]})
    forest = pd.DataFrame(rows)
    y = np.arange(len(forest))
    ax.errorbar(forest["mean"], y,
                xerr=[forest["mean"] - forest["low"], forest["high"] - forest["mean"]],
                fmt="o", markersize=3, linewidth=0.9, color=palette.get("neutral", "#4A4A4A"))
    ax.axvline(0.0, color="#999999", linewidth=0.7, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels(forest["label"], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("Settlement − wild (probability difference)")
    ax.set_title("Posterior contrast, 95% HDI", fontsize=7)
    sources["Figure2c_contrast"] = forest

    # (d) sensitivity
    ax = fig.add_subplot(gs[1, 1])
    _panel_label(ax, "d")
    if sensitivity is not None and len(sensitivity):
        sens = sensitivity[sensitivity["status"] == "ok"].copy()
        if len(sens):
            y = np.arange(len(sens))
            ax.errorbar(sens["mean_difference"], y,
                        xerr=[sens["mean_difference"] - sens["hdi_low"],
                              sens["hdi_high"] - sens["mean_difference"]],
                        fmt="s", markersize=3, linewidth=0.9,
                        color=palette.get("accent", "#3B6EA5"))
            ax.axvline(0.0, color="#999999", linewidth=0.7, linestyle="--")
            ax.set_yticks(y)
            ax.set_yticklabels(
                [f"{str(d)[:5]} · {v}" for d, v in zip(sens["dimension"], sens["variant"])],
                fontsize=5.5,
            )
            ax.invert_yaxis()
        not_estimable = sensitivity[sensitivity["status"] != "ok"]
        if len(not_estimable):
            ax.text(0.02, -0.22,
                    "Not estimable: " + ", ".join(not_estimable["variant"].astype(str)),
                    transform=ax.transAxes, fontsize=5.5, color="#B4503C")
        sources["Figure2d_sensitivity"] = sensitivity
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "Sensitivity analyses not generated", ha="center", va="center",
                fontsize=7, color="#B4503C")
    ax.set_xlabel("Settlement − wild")
    ax.set_title("Sensitivity", fontsize=7)

    fig.suptitle("Figure 2 · Bayesian inference of habitat-associated detection",
                 fontsize=9, y=0.98)
    fig.text(0.5, 0.005,
             "Observational association standardized to a common sampling distribution; "
             "not a causal effect of urbanization.",
             ha="center", fontsize=5.5, style="italic", color="#555555")
    if store.is_synthetic():
        _watermark(fig, style.get("synthetic_watermark", "SYNTHETIC DATA"))

    return {"figure": "Figure2", "status": "ok",
            "files": _save(fig, out_dir, "Figure2", style), "source_tables": sources,
            "n": {"posterior_draws": len(draws)}}


# ---------------------------------------------------------------------------
# Figure 3 — connectivity
# ---------------------------------------------------------------------------
def figure3(store: ArtifactStore, style: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    conn = store.connectivity()
    meta = store.connectivity_metadata()
    perm = store.connectivity_permutation()
    frame = store.analysis_frame()

    missing = [n for n, v in (("connectivity", conn), ("connectivity metadata", meta))
               if v is None or (hasattr(v, "__len__") and not len(v))]
    fig = _new_figure(style)
    if missing:
        _not_generated(fig, "Figure 3: Cross-host connectivity", missing)
        return {"figure": "Figure3", "status": "not_generated", "missing": missing,
                "files": _save(fig, out_dir, "Figure3", style), "source_tables": {}}

    palette = style.get("palette", {})
    gs = GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.4)
    sources: dict[str, pd.DataFrame] = {}

    # (a) host-virus bipartite structure
    ax = fig.add_subplot(gs[0, 0])
    _panel_label(ax, "a")
    if frame is not None and len(frame):
        counts = (
            frame[frame["y"] > 0]
            .groupby(["host_species", "virus_family"]).size()
            .rename("n").reset_index()
        )
        pivot = counts.pivot(index="host_species", columns="virus_family", values="n").fillna(0)
        im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="Greys")
        ax.set_xticks(range(pivot.shape[1]))
        ax.set_xticklabels(pivot.columns, rotation=40, ha="right", fontsize=5)
        ax.set_yticks(range(pivot.shape[0]))
        ax.set_yticklabels(pivot.index, fontsize=5, style="italic")
        fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02).set_label("detections", fontsize=5)
        sources["Figure3a_bipartite"] = pivot.reset_index()
    ax.set_title("Observed host × viral family detections", fontsize=7)

    # (b) connectivity matrix
    ax = fig.add_subplot(gs[0, 1])
    _panel_label(ax, "b")
    habitat = "settlement" if "settlement" in (meta.get("matrices") or {}) else next(
        iter(meta.get("matrices") or {"none": None})
    )
    matrices = (meta.get("matrices") or {}).get(habitat)
    if matrices:
        hosts = matrices["hosts"]
        C = np.asarray(matrices["mean"], dtype=float)
        im = ax.imshow(C, cmap="magma")
        ax.set_xticks(range(len(hosts)))
        ax.set_xticklabels(hosts, rotation=40, ha="right", fontsize=5, style="italic")
        ax.set_yticks(range(len(hosts)))
        ax.set_yticklabels(hosts, fontsize=5, style="italic")
        fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02).set_label(
            "expected shared lineages", fontsize=5
        )
        sources["Figure3b_matrix"] = pd.DataFrame(C, index=hosts, columns=hosts).reset_index()
    ax.set_title(f"Posterior C(i,j) — {habitat}", fontsize=7)

    # (c) observed vs null
    ax = fig.add_subplot(gs[1, 0])
    _panel_label(ax, "c")
    pairs = conn[~conn["self_pair"]] if "self_pair" in conn else conn
    if "mean_expected_shared_lineages_null" in pairs:
        ax.scatter(pairs["mean_expected_shared_lineages_null"],
                   pairs["mean_expected_shared_lineages"],
                   s=6, alpha=0.7, color=palette.get("accent", "#3B6EA5"))
        lim = [0, float(max(pairs["mean_expected_shared_lineages"].max(),
                            pairs["mean_expected_shared_lineages_null"].max())) * 1.05]
        ax.plot(lim, lim, color="#999999", linewidth=0.7, linestyle="--")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel("Marginal null (no host × family term)")
        ax.set_ylabel("Fitted model")
    sources["Figure3c_null"] = pairs
    if perm is not None and len(perm):
        sources["Figure3c_permutation_null"] = perm
    ax.set_title("Connectivity against its null", fontsize=7)

    # (d) connectivity vs phylogenetic distance
    ax = fig.add_subplot(gs[1, 1])
    _panel_label(ax, "d")
    phylo = (meta.get("phylogeny_association") or {}).get(habitat, {})
    if phylo.get("status") == "ok":
        ax.text(0.5, 0.6,
                f"slope = {phylo['slope_mean']:.3g}\n"
                f"95% HDI [{phylo['slope_hdi_low']:.3g}, {phylo['slope_hdi_high']:.3g}]\n"
                f"{phylo['n_host_species']} species, {phylo['n_pairs']} pairs",
                ha="center", va="center", fontsize=7, transform=ax.transAxes)
        ax.axis("off")
        sources["Figure3d_phylo"] = pd.DataFrame([phylo])
    else:
        ax.axis("off")
        ax.text(0.5, 0.5,
                "Host phylogenetic distances not available\n"
                f"({phylo.get('reason', 'no host phylogeny supplied')})",
                ha="center", va="center", fontsize=6.5, color="#B4503C")
    ax.set_title("Connectivity vs host relatedness", fontsize=7)

    fig.suptitle("Figure 3 · Posterior cross-host connectivity", fontsize=9, y=0.98)
    fig.text(0.5, 0.005,
             "Expected co-detection under a common sampling design; "
             "not inferred host-to-host transmission.",
             ha="center", fontsize=5.5, style="italic", color="#555555")
    if store.is_synthetic():
        _watermark(fig, style.get("synthetic_watermark", "SYNTHETIC DATA"))

    return {"figure": "Figure3", "status": "ok",
            "files": _save(fig, out_dir, "Figure3", style), "source_tables": sources,
            "n": {"host_pairs": len(conn)}}


# ---------------------------------------------------------------------------
# Figure 4 — evolution
# ---------------------------------------------------------------------------
def figure4(store: ArtifactStore, style: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    selection = store.selection_results()
    lineages = store.lineage_results()
    status = store.evolution_status()

    missing: list[str] = []
    if status is None or status.get("status") in {None, "not_run"}:
        missing.append("evolutionary analyses (not run)")
    if lineages is None or not len(lineages):
        missing.append("per-lineage phylogenetic results")

    fig = _new_figure(style)
    if missing:
        if status and status.get("missing_tools"):
            missing.append("external tools: " + ", ".join(status["missing_tools"]))
        _not_generated(fig, "Figure 4: Viral evolutionary structure", missing)
        return {"figure": "Figure4", "status": "not_generated", "missing": missing,
                "files": _save(fig, out_dir, "Figure4", style), "source_tables": {}}

    palette = style.get("palette", {})
    gs = GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)
    sources: dict[str, pd.DataFrame] = {"Figure4a_tree": lineages}

    ax = fig.add_subplot(gs[0, 0])
    _panel_label(ax, "a")
    ax.axis("off")
    trees = lineages[lineages.get("tree_status") == "ok"] if "tree_status" in lineages else lineages
    ax.text(0.5, 0.5,
            f"{len(trees)} lineage phylogenies inferred\n"
            "Newick files are published alongside this figure",
            ha="center", va="center", fontsize=7)
    ax.set_title("Focal phylogenies", fontsize=7)

    ax = fig.add_subplot(gs[0, 1])
    _panel_label(ax, "b")
    col = "alignment_mean_pairwise_divergence"
    if col in lineages and lineages[col].notna().any():
        ax.bar(range(len(lineages)), lineages[col].fillna(0),
               color=palette.get("neutral", "#4A4A4A"))
        ax.set_xticks(range(len(lineages)))
        ax.set_xticklabels(lineages["lineage"], rotation=60, ha="right", fontsize=5)
        ax.set_ylabel("Mean pairwise divergence")
        sources["Figure4b_distance"] = lineages[["lineage", col]]
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "Divergence not computed", ha="center", va="center",
                fontsize=7, color="#B4503C")
    ax.set_title("Sequence divergence", fontsize=7)

    ax = fig.add_subplot(gs[1, 0])
    _panel_label(ax, "c")
    if selection is not None and len(selection):
        ok = selection[selection["status"] == "ok"]
        if len(ok):
            ax.bar(range(len(ok)), ok["n_significant"].fillna(0),
                   color=palette.get("accent", "#3B6EA5"))
            ax.set_xticks(range(len(ok)))
            ax.set_xticklabels([f"{lin}\n{meth}" for lin, meth in zip(ok["lineage"], ok["method"])],
                               rotation=60, ha="right", fontsize=4.5)
            ax.set_ylabel("Significant sites")
        withheld = selection[selection["status"] != "ok"]
        if len(withheld):
            ax.text(0.02, -0.35,
                    f"{len(withheld)} analyses withheld (ineligible alignment or missing tool)",
                    transform=ax.transAxes, fontsize=5, color="#B4503C")
        sources["Figure4c_selection"] = selection
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "Selection analyses withheld", ha="center", va="center",
                fontsize=7, color="#B4503C")
    ax.set_title("Selection analyses", fontsize=7)

    ax = fig.add_subplot(gs[1, 1])
    _panel_label(ax, "d")
    breadth = store.table("evolution", "evolution_host_breadth")
    if breadth is not None and len(breadth) and "host_breadth" in breadth:
        ax.scatter(breadth.get(col, pd.Series(np.nan, index=breadth.index)),
                   breadth["host_breadth"], s=10, color=palette.get("neutral", "#4A4A4A"))
        ax.set_xlabel("Mean pairwise divergence")
        ax.set_ylabel("Host species detected in")
        sources["Figure4d_breadth"] = breadth
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "Host breadth integration not generated", ha="center",
                va="center", fontsize=7, color="#B4503C")
    ax.set_title("Host breadth vs divergence", fontsize=7)

    fig.suptitle("Figure 4 · Viral evolutionary structure across hosts", fontsize=9, y=0.98)
    if store.is_synthetic():
        _watermark(fig, style.get("synthetic_watermark", "SYNTHETIC DATA"))
    return {"figure": "Figure4", "status": "ok",
            "files": _save(fig, out_dir, "Figure4", style), "source_tables": sources,
            "n": {"lineages": len(lineages)}}


# ---------------------------------------------------------------------------
# Figure 5 — validation
# ---------------------------------------------------------------------------
def figure5(store: ArtifactStore, style: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    comparison = store.model_comparison()
    folds = store.cv_folds()
    calibration = store.calibration("site_id")
    ppc = store.posterior_predictive()

    missing = [n for n, v in (("model comparison", comparison),
                              ("cross-validation folds", folds),
                              ("posterior predictive check", ppc))
               if v is None or (hasattr(v, "__len__") and not len(v))]
    fig = _new_figure(style)
    if missing:
        _not_generated(fig, "Figure 5: Validation", missing)
        return {"figure": "Figure5", "status": "not_generated", "missing": missing,
                "files": _save(fig, out_dir, "Figure5", style), "source_tables": {}}

    palette = style.get("palette", {})
    gs = GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)
    sources: dict[str, pd.DataFrame] = {}

    ax = fig.add_subplot(gs[0, 0])
    _panel_label(ax, "a")
    comp = comparison.sort_values("elpd_loo")
    ax.errorbar(comp["elpd_loo"], range(len(comp)),
                xerr=comp["se"].fillna(0), fmt="o", markersize=3, linewidth=0.9,
                color=palette.get("neutral", "#4A4A4A"))
    ax.set_yticks(range(len(comp)))
    ax.set_yticklabels(comp["model"], fontsize=6)
    ax.set_xlabel("ELPD (LOO)")
    ax.set_title("Model comparison", fontsize=7)
    sources["Figure5a_comparison"] = comparison

    ax = fig.add_subplot(gs[0, 1])
    _panel_label(ax, "b")
    if calibration is not None and len(calibration):
        ax.plot([0, calibration["predicted_mean"].max()], [0, calibration["predicted_mean"].max()],
                color="#999999", linewidth=0.7, linestyle="--")
        ax.errorbar(calibration["predicted_mean"], calibration["observed_frequency"],
                    yerr=[calibration["observed_frequency"] - calibration["lower"],
                          calibration["upper"] - calibration["observed_frequency"]],
                    fmt="o", markersize=3, linewidth=0.8,
                    color=palette.get("accent", "#3B6EA5"))
        ax.set_xlabel("Predicted probability")
        ax.set_ylabel("Observed frequency")
        sources["Figure5b_calibration"] = calibration
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "Calibration not generated", ha="center", va="center",
                fontsize=7, color="#B4503C")
    ax.set_title("Held-out calibration", fontsize=7)

    ax = fig.add_subplot(gs[1, 0])
    _panel_label(ax, "c")
    rep = ppc.get("replicate_detection_rate") or []
    if rep:
        ax.hist(rep, bins=30, color=palette.get("light", "#C9C9C9"),
                edgecolor="white", linewidth=0.3)
        ax.axvline(ppc.get("observed_detection_rate", np.nan), color="#B4503C", linewidth=1.2)
        ax.set_xlabel("Replicate detection rate")
        ax.set_ylabel("Count")
        sources["Figure5c_ppc"] = pd.DataFrame({"replicate_detection_rate": rep})
    ax.set_title("Posterior predictive check", fontsize=7)

    ax = fig.add_subplot(gs[1, 1])
    _panel_label(ax, "d")
    ax.scatter(folds["brier"], folds["mean_lpd"], s=12,
               color=palette.get("neutral", "#4A4A4A"))
    for _, r in folds.iterrows():
        ax.annotate(str(r["fold"]), (r["brier"], r["mean_lpd"]), fontsize=5,
                    xytext=(2, 2), textcoords="offset points")
    ax.set_xlabel("Brier score (lower is better)")
    ax.set_ylabel("Mean log predictive density")
    ax.set_title("Held-out site folds", fontsize=7)
    sources["Figure5d_heldout"] = folds

    fig.suptitle("Figure 5 · Model validation and predictive performance", fontsize=9, y=0.98)
    if store.is_synthetic():
        _watermark(fig, style.get("synthetic_watermark", "SYNTHETIC DATA"))
    return {"figure": "Figure5", "status": "ok",
            "files": _save(fig, out_dir, "Figure5", style), "source_tables": sources,
            "n": {"folds": len(folds)}}


FIGURE_BUILDERS: dict[str, Callable[..., dict[str, Any]]] = {
    "Figure1": figure1,
    "Figure2": figure2,
    "Figure3": figure3,
    "Figure4": figure4,
    "Figure5": figure5,
}


def build_figures(store: ArtifactStore, config: dict[str, Any],
                  out_dir: str | Path) -> dict[str, Any]:
    style = config.get("style", {}) or {}
    _apply_style(style)
    figures_dir = Path(out_dir) / "figures"
    tables_dir = Path(out_dir) / "tables"
    results: dict[str, Any] = {}

    for name in (config.get("figures") or {}):
        builder = FIGURE_BUILDERS.get(name)
        if builder is None:
            continue
        try:
            result = builder(store, style, figures_dir)
        except Exception as exc:  # a broken panel must not lose the other figures
            fig = _new_figure(style)
            _not_generated(fig, f"{name}: generation failed",
                           [f"{type(exc).__name__}: {exc}"])
            result = {"figure": name, "status": "failed", "error": str(exc),
                      "files": _save(fig, figures_dir, name, style), "source_tables": {}}
        for table_name, df in (result.pop("source_tables", {}) or {}).items():
            if df is not None and len(df):
                write_table(df, tables_dir / table_name)
        results[name] = result

    write_json(results, Path(out_dir) / "figure_manifest.json")
    return results
