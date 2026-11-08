# Statistical methods

## 1. Observation unit and response

The unit of observation is a **sequencing library × viral lineage unit** cell.
A library is one sequenced nucleic-acid preparation, derived from a single
animal or a documented pool. A lineage unit is a published viral lineage
(vOTU), or a viral family when lineages too rare to support a lineage-specific
effect are pooled (`min_libraries_detected` in
[`config/viral_filters.yaml`](../config/viral_filters.yaml)).

Each cell takes one of three states, and the distinction is enforced by the
schema layer:

| State | Meaning | Response |
| --- | --- | --- |
| `positive` | Detection passing the mapping-support thresholds, in a screened library | `y = 1` |
| `negative` | Screened against the common panel at adequate depth, not detected | `y = 0` |
| `unassessable` | Not screened for that lineage, or below the effective-depth floor | missing |

**Unassessable cells never enter the likelihood and are never imputed as
zeros.** A positive observed in an unscreened library is recorded as
`positive_outside_assessed_universe` and excluded, because counting it while
discarding that library's zeros would bias detection upward precisely where
sequencing effort was lowest.

## 2. Likelihood

For library *l* and lineage unit *v*:

$$Y_{lv} \sim \mathrm{Bernoulli}(p_{lv})$$

$$
\mathrm{logit}(p_{lv}) = \alpha
+ a_{f(v)} + a_v
+ \beta_U U_l
+ b_{h(l)} + \kappa_{h(l),f(v)} + c_{t(l)}
+ g_{r(l)} + s_{\mathrm{batch}(l)}
+ \gamma_D \log D_l + \gamma_N \log N_l + \gamma_T T_l
+ \boldsymbol{\theta}^\top E_l
$$

| Term | Meaning |
| --- | --- |
| $\alpha$ | Intercept |
| $a_{f(v)}, a_v$ | Viral family and lineage varying intercepts |
| $\beta_U U_l$ | Habitat: settlement (1) versus wild (0) |
| $b_{h(l)}$ | Bat host species |
| $\kappa_{h(l),f(v)}$ | Shrinkage-regularized host × viral family interaction |
| $c_{t(l)}$ | Tissue group |
| $g_{r(l)}, s_{\mathrm{batch}(l)}$ | Sampling site and sequencing batch |
| $\gamma_D \log D_l$ | Effective sequencing depth |
| $\gamma_N \log N_l$ | Pool size (animals contributing to the library) |
| $\gamma_T T_l$ | Standardized collection time |
| $\boldsymbol{\theta}^\top E_l$ | Supported site-level ecological covariates |

The **host × family interaction** is not optional decoration. Without it,
host-sharing predictions are dominated by overall differences in host and
lineage detection frequency rather than by lineage-specific host
associations, and the connectivity analysis measures sampling rather than
biology.

The default model uses a **common habitat coefficient**. Lineage-varying
habitat effects are added only if the data support them — they are not in the
default ladder.

### Effective depth

The depth denominator must be defined identically for every library
(`non_host_non_rrna_reads` by default). If the published tables report only
total reads, the fallback is used and the substitution is recorded in the QC
report and the model metadata, because it changes what the depth coefficient
means.

## 3. Model ladder

All members share the same eligible-observation set, so predictive
comparisons are well defined.

| Model | Terms added | Question |
| --- | --- | --- |
| M0 | Viral identity, depth, pool size, time, batch | What does sequencing alone predict? |
| M1 | + host species, tissue | Does sample composition explain it? |
| M2 | + site | Is habitat confounded with place? |
| M3 | + habitat, host × family | **Primary habitat comparison** |
| M4 | + ecological covariates | Does measured ecology explain the habitat term? |
| M5 | + host phylogenetic covariance | Does host relatedness explain host effects? |

M4 is fitted only when ecological covariates cover enough modelled cells; M5
only when a bat phylogeny is supplied. An unmet requirement is recorded as a
skipped model, not silently dropped.

## 4. Priors

See [`config/priors.yaml`](../config/priors.yaml). All continuous predictors
are standardized on the modelled subset, so priors are on the standardized
scale.

| Parameter | Prior |
| --- | --- |
| $\alpha$ | Normal(−3, 1.5) |
| $\beta_U$, $\gamma_D$, $\gamma_N$ | Normal(0, 1) |
| $\gamma_T$, $\boldsymbol{\theta}$ | Normal(0, 0.5) |
| Family / lineage / host / tissue / site scales | HalfNormal(1) |
| Batch scale | HalfNormal(0.75) |
| Host × family interaction scale | HalfNormal(0.5) |

Variants `tight` and `wide` halve and roughly double the key scales for the
sensitivity analysis.

**Prior predictive check.** Detections are simulated under the priors alone
before any data is used. If the implied marginal detection rate falls outside
a plausible band (0.05%–35%), the run warns: priors that imply implausible
virome prevalence are not weakly informative, they are wrong.

## 5. Inference

PyMC with NUTS, 4 chains, 1000 warmup + 1000 post-warmup draws,
`target_accept = 0.9`. Varying effects use a non-centred parameterization.

Where PyMC is unavailable, a self-contained engine runs the same
`ModelDesign`: Pólya-Gamma data augmentation (Polson, Scott & Windle 2013)
gives a closed-form Gaussian conditional for the coefficients, and each
hierarchical scale is updated by slice sampling on $\log\sigma$, interweaved
between the centered and non-centered parameterizations (ASIS). Its documented
limitations — a truncated Pólya-Gamma representation and no divergence
diagnostic — are recorded in the output, and the engine is named in every
legend. Published results should use PyMC/NUTS.

### Convergence gates

Results are not released unless: $\hat{R} \le 1.01$ for every parameter, bulk
and tail ESS $\ge 400$, divergent transitions $\le 0.5\%$ of draws, and
maximum-tree-depth hits $\le 5\%$. R-hat and ESS follow the split-chain
rank-normalized definitions of Vehtari et al. (2021).

## 6. Primary estimand

> The difference in **average lineage-level detection probability** between
> settlement and wild habitats, standardized to one common sampling
> distribution.

Concretely: posterior predictions are formed over a reference population of
library profiles — restricted to host species observed in **both** habitats —
crossed with every modelled lineage, holding tissue, depth, pool size, time
and ecology at their observed values, at a **new** site and a **new**
sequencing batch drawn from their estimated population distributions. The
habitat indicator is then set to 1 and to 0 using the *same* new-site draws,
so the contrast is paired.

$$\Delta = \overline{p}(U=1) - \overline{p}(U=0)$$

Reported with its posterior mean, median, 95% HDI and the posterior
probability of a positive difference; also as $\beta_U$ on the log-odds scale
and as an odds ratio, and stratified by viral family and host species.

**This is an observational association, not a causal effect of urbanization.**
Habitat is not randomized, sites are not interchangeable, and the
standardization removes only measured and modelled differences in sampling.
That sentence is attached to every artifact the contrast produces.

### Sensitivity

Across prior variants (`regularizing`, `tight`, `wide`), lineage resolutions
(`lineage`, `family`, `strict`) and sample restrictions (`overlap_only`,
`no_pools`, `drop_rare_lineages`, `depth_matched`). A restriction that leaves
too few cells or one habitat is reported as `not_estimable`, not as a null
result.

## 7. Cross-host connectivity

For host species *i*, *j* and habitat *u*:

$$C_{ij}^{(u)} = \sum_{v \in \mathcal{V}} p_{iv}^{(u)} p_{jv}^{(u)}$$

with probabilities predicted for one common sampling design. This is the
expected number of lineages that two hypothetical libraries, one from each
host, would both be positive for — assuming conditional independence given
the parameters.

It is **expected co-detection, not transmission**. Two species can co-detect
because they share an environment, a diet, a laboratory batch, or simply a
common viral family.

Two nulls:

1. **Marginal null** — the same model with the host × family interaction
   removed. Excess over this null reflects lineage-specific host association
   rather than overall host or lineage detection frequency.
2. **Permutation null** — host species labels shuffled among libraries
   *within* habitat and site, preserving where and how much was sequenced
   while destroying the host–lineage association.

Connectivity is also related to host phylogenetic distance when a bat
phylogeny is supplied; otherwise that panel is marked not generated.

## 8. Predictive validation

**Model validity:** prior predictive simulation, parameter recovery against a
known generating process, MCMC diagnostics, posterior predictive checks, and
sensitivity to rare lineages, pooling, classification thresholds and
geographic overlap.

**Predictive validity:** site- and region-held-out cross-validation scored by
log predictive density, Brier score, AUC and calibration curves. **Held-out
site and batch random effects are integrated over their estimated population
distributions, never estimated from the held-out observations** — estimating
them would leak the answer and overstate generalization.

**Independent comparison:** the individual-bat dataset (PRJNA929070) is used
to ask whether host-sharing summaries remain interpretable without pooled
sequencing. It is not presented as an independent validation of the habitat
contrast.

## 9. Model comparison

PSIS-LOO via ArviZ, reported as ELPD with standard errors, pairwise
differences with their standard errors, and Pareto-k warnings. Without ArviZ
the pipeline falls back to plain importance-sampled LOO and says so
(`method: is_loo`), with no Pareto-k diagnostics claimed.

## References

- Polson, Scott & Windle (2013). Bayesian inference for logistic models using
  Pólya-Gamma latent variables. *JASA* 108(504):1339–1349.
- Vehtari, Gelman, Simpson, Carpenter & Bürkner (2021). Rank-normalization,
  folding, and localization: an improved R-hat. *Bayesian Analysis* 16(2).
- Vehtari, Gelman & Gabry (2017). Practical Bayesian model evaluation using
  leave-one-out cross-validation and WAIC. *Statistics and Computing* 27.
- Yu & Meng (2011). To center or not to center: that is not the question —
  an ancillarity-sufficiency interweaving strategy. *JCGS* 20(3).
