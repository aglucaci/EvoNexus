# Manuscript outline

**Working title:** Beyond viral sharing — Bayesian inference of urban ecology
and cross-host viral detection in bats

This outline is a scaffold for writing. Numbers belong in the generated
tables and legends, not here; every claim below should be checkable against
an artifact the pipeline produced.

## Abstract
- What is known: bat metatranscriptomic surveys report large viral
  compendia and apparent host sharing.
- The gap: detection is confounded with sampling. Host composition, tissue,
  sequencing depth, pooling and geography differ systematically between
  habitats, and absence is often not distinguishable from "not assessed".
- What we do: a hierarchical Bayesian model of library-level viral RNA
  detection, with a habitat contrast standardized to a common sampling
  distribution, and posterior cross-host connectivity with
  sampling-preserving nulls.
- What we find: `results/*/tables/Table2` and `Table3`.
- What it means, and does not: an observational association, not a causal
  effect of urbanization; expected co-detection, not transmission.

## Introduction
1. Bats, RNA viromes and the interest in human-modified habitats.
2. Why library-level metatranscriptomic detection is not individual-animal
   infection, and why pooled libraries break the usual denominator.
3. Why unadjusted host-sharing counts mostly measure sampling effort.
4. The estimand this paper actually targets, stated before any result.

## Results
1. **What was sampled, and where a contrast is estimable** — Figure 1.
   Geography, host × habitat overlap, sequencing effort, and the exclusion
   waterfall from source records to modelled cells.
2. **Habitat-associated detection** — Figure 2, Table 2. Raw rates, the
   standardized posterior contrast overall and by viral family, and its
   sensitivity to priors, lineage resolution and sample restrictions.
3. **Cross-host connectivity** — Figure 3, Table 3. Posterior expected
   co-detection, the marginal null without the host × family interaction, the
   permutation null, and the relationship to host relatedness.
4. **Viral evolutionary structure** — Figure 4, Supplementary 4. Only for
   lineages that passed the alignment and selection eligibility gates.
5. **Validation** — Figure 5, Table 4. Model comparison, calibration,
   posterior predictive checks, and held-out-site prediction.

## Discussion
- What the standardized contrast does and does not license.
- Why connectivity is co-detection: shared environment, diet, laboratory
  batch and shared viral family all produce it without transmission.
- Pooling, and what it costs in interpretability.
- The unassessable fraction as a measurable limit on the whole analysis.
- Generalization: the same framework applied to urban metatranscriptomic
  surveys, wastewater viromics and other host–environment surveillance.

## Limitations
Written from `qc_report.md` and the exclusion waterfall, not from memory:
- Observational design; habitat is not randomized and sites are not
  interchangeable.
- Library-level detection, not infection.
- Pooled and mixed-species libraries.
- Cells that could not be assessed, and what that does to the denominator.
- Reference-panel dependence of what counts as a lineage.
- Any gate that failed or warned.

## Methods
Condensed from [`docs/STATISTICAL_METHODS.md`](../docs/STATISTICAL_METHODS.md)
and [`methods.md`](methods.md).

## Data and code availability
- Source datasets and accessions: [`docs/DATA_SOURCES.md`](../docs/DATA_SOURCES.md)
- Code, configuration and the release manifest:
  [`docs/REPRODUCIBILITY.md`](../docs/REPRODUCIBILITY.md)

## Display items
| Item | Generated at | Legend |
| --- | --- | --- |
| Figure 1 | `results/*/figures/Figure1.*` | `figure_legends/Figure1_legend.md` |
| Figure 2 | `results/*/figures/Figure2.*` | `figure_legends/Figure2_legend.md` |
| Figure 3 | `results/*/figures/Figure3.*` | `figure_legends/Figure3_legend.md` |
| Figure 4 | `results/*/figures/Figure4.*` | `figure_legends/Figure4_legend.md` |
| Figure 5 | `results/*/figures/Figure5.*` | `figure_legends/Figure5_legend.md` |
| Tables 1–4, S1–S5 | `results/*/tables/`, `results/*/manuscript/` | in the table files |
