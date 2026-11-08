# Methods (manuscript draft)

Draft prose for the manuscript. Every number is left as a placeholder of the
form `{{artifact}}`; fill them from the generated tables rather than by hand,
and keep the qualifications — they are what make the claims defensible.

## Data sources

Processed viral catalog, ecological metadata and sequencing-library metadata
were obtained from the BtCN-Virome deposit (Yan et al. 2025; Figshare
10.6084/m9.figshare.24515056). Sequencing run metadata were retrieved from
NCBI SRA BioProjects PRJNA862267 and PRJNA997461. An independent
individual-bat dataset (PRJNA929070) was used for secondary evaluation only.
Every file was recorded in an acquisition manifest with its URL, retrieval
date and SHA-256 checksum (Supplementary Table 1).

## Metadata harmonization

Source tables were mapped onto a canonical schema of sites, hosts, libraries,
viral lineages, detections and ecological covariates, with a crosswalk
recording the source table and column behind every canonical identifier. The
mapping was reviewed and approved against the real table headers before any
transformation was applied.

Habitat was reduced to a controlled vocabulary (settlement, wild, unknown).
Values not matching the vocabulary were retained as `unknown` and excluded
from the habitat contrast rather than assigned to a category. Libraries
derived from mixed-species pools were flagged and excluded from the
single-host model; they were not represented as one independently sampled
animal per constituent species.

## Detection calls

A library × lineage cell was called positive when mapping support met the
preregistered thresholds (≥ {{min_mapped_reads}} mapped reads and ≥
{{min_coverage_breadth}} genome coverage breadth). A cell was treated as a
valid negative only when that library had been screened against the common
reference panel at an effective depth of at least {{min_effective_depth}}
reads. All other cells were recorded as unassessable and excluded from the
likelihood; they were not imputed as zeros. Detections observed in libraries
that did not meet the screening criteria were recorded but excluded, because
including them while discarding the same libraries' negatives would bias
detection upward where sequencing effort was lowest. Of {{n_cells_total}}
cells, {{n_positive}} were positive, {{n_negative}} were assay-tested
negatives and {{n_unassessable}} were unassessable.

## Quality gates

Seven gates were evaluated before modelling: identifier integrity, documented
biological structure, detection validity, habitat exposure, host overlap
across habitats, sample dependence, and viral classification. The primary
analysis proceeded only because all critical gates were satisfied
(Supplementary Table 1).

## Statistical model

We modelled the probability that viral lineage *v* is detected in library *l*
as a Bernoulli outcome with a logit link, including varying intercepts for
viral family and lineage, bat host species, tissue group, sampling site and
sequencing batch; a shrinkage-regularized host-by-viral-family interaction;
and fixed effects for habitat, log effective sequencing depth, log pool size,
standardized collection time and supported site-level ecological covariates.
Continuous predictors were standardized. Priors were weakly informative and
regularizing (Supplementary Table 2), and prior predictive simulation
confirmed they imply plausible marginal detection rates before seeing data.

Posterior inference used {{n_chains}} NUTS chains of {{n_draws}} post-warmup
draws after {{n_tune}} warmup iterations. Convergence was assessed by
split rank-normalized R̂ (maximum {{max_rhat}}), bulk and tail effective
sample size (minimum {{min_ess}}) and divergent transitions ({{n_divergences}});
posterior predictive checks are shown in Figure 5 and Supplementary Table 3.

## Primary estimand

The primary estimand is the difference in average lineage-level detection
probability between settlement and wild habitats, standardized to a common
sampling distribution: the same host species — restricted to those observed
in both habitats — the same tissues, sequencing depths, pool sizes and
collection times, at a new site and new sequencing batch integrated over
their estimated population distributions. The estimated difference was
{{contrast_mean}} (95% highest-density interval {{contrast_hdi}}).

**This is an observational association, not a causal effect of
urbanization.** Habitat was not randomized, sites are not interchangeable,
and standardization removes only measured and modelled differences in
sampling. Robustness across prior specifications, lineage resolutions and
sample restrictions is reported in Figure 2d and Supplementary Table 5.

## Cross-host connectivity

For each pair of host species we computed the posterior expected number of
lineages co-detected under a common sampling design, summing the product of
predicted detection probabilities across lineages. This quantity assumes
conditional independence of two hypothetical libraries given the parameters
and measures expected co-detection, **not host-to-host transmission**. It was
compared with two nulls: a marginal null removing the host-by-family
interaction, and a permutation null shuffling host labels within habitat and
site, which preserves the sampling structure.

## Predictive validation

Model performance was assessed by leave-one-out cross-validation (PSIS-LOO)
across the model ladder, and by held-out-site and held-out-region
cross-validation scored with log predictive density, Brier score and
calibration. Random effects of held-out sites and batches were integrated
over their estimated population distributions rather than estimated from the
held-out observations.

## Evolutionary analyses

Lineages with sufficient homologous coding sequence were aligned with MAFFT,
screened for recombination, and used to infer maximum-likelihood phylogenies
with IQ-TREE 2. Selection was tested with HyPhy only for alignments meeting
preregistered eligibility criteria for sequence number, alignment length,
divergence and absence of recombination. Lineages failing any criterion are
reported as documented exclusions in Supplementary Table 4 rather than
assigned uncertain estimates.

## Reproducibility

The analysis was executed as a Nextflow DSL2 pipeline. Software versions,
random seeds, configuration and input checksums for the release that produced
these results are given in Supplementary Table 1 and in the repository's
reproducibility documentation.
