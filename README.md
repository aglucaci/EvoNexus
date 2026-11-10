# EvoNexus

**Evolutionary and Ecological Networks of Viral Host Sharing**

A provenance-aware Nextflow DSL2 pipeline that downloads public bat virome
datasets, harmonizes their metadata, runs Bayesian ecological and evolutionary
analyses, and generates manuscript figures, tables and legends from validated
analysis outputs.

---

## The question

Across bat sequencing libraries, does the probability of detecting a given
viral RNA lineage differ between human-settlement and wild habitats, once the
comparison is standardized to the same host species, tissues, sequencing
effort and sampling design — and which host species are expected to co-detect
the same lineages under that common design?

## What EvoNexus does not claim

These limits are enforced in code, not only in prose. They appear in the QC
report, in every figure legend, and in the model metadata.

- **Detections are library-level viral RNA detections.** They are not
  confirmed infections of individual animals. A pooled library is not one
  independently sampled animal per constituent species, and mixed-species
  pools are excluded from the single-host model rather than split.
- **Absence is not evidence of absence.** A library × lineage cell counts as a
  negative only when that library was screened against the common reference
  panel at adequate sequencing depth. Everything else stays `unassessable`
  and carries no response. Positives from unscreened libraries are recorded
  but kept out of the likelihood, because counting them while discarding the
  same library's zeros would bias detection upward exactly where effort was
  lowest.
- **The habitat contrast is an observational association**, standardized to a
  common sampling distribution. It is not a causal effect of urbanization.
- **Connectivity is expected co-detection** under the fitted model, not
  inferred host-to-host transmission.
- **Synthetic results are never presented as empirical ones.** They are
  watermarked across every panel and labelled in every table and legend.

## Design principle

Data provenance and sample harmonization are validated *before* any Bayesian
model runs. An incorrect relationship between an SRA run, a sequencing
library, a bat species and a viral lineage would invalidate every downstream
estimate, so the pipeline refuses to proceed past that point on trust.

```
discovery → download → harmonize → QC gates ─┬─ pass → inference → manuscript
                                             └─ fail → descriptive only + report
```

---

## Quick start

```bash
git clone https://github.com/aglucaci/EvoNexus.git
cd EvoNexus
conda env create -f environment.yml && conda activate evonexus
```

Conda is recommended because the evolution workflow needs external tools
(MAFFT, IQ-TREE 2, HyPhy) that pip cannot install. For the Python side
alone:

```bash
python -m pip install -r requirements.txt
```

Run the deterministic synthetic test — no downloads, a few minutes on a laptop:

```bash
nextflow run . -profile test
```

Nextflow 24.04 or newer is required. Note that Nextflow bundles its own Groovy
runtime, so it must also match your JDK: Java 25 needs Nextflow 26+, while
Nextflow 24.x and 25.x run on Java 17–24. A mismatch fails with
`Unsupported class file major version` before any pipeline code is read.

Or drive the same steps directly, without Nextflow:

```bash
PYTHONPATH=scripts python -m evonexus.cli acquire --out-dir results/test/acquisition --synthetic --offline
```

Check which sampler is available in your environment:

```bash
PYTHONPATH=scripts python -m evonexus.cli engines
```

---

## Running on real data

The real analysis is intentionally **not** a single command, because step 2
requires a human decision.

**1. Acquire and inventory.** Nothing is interpreted yet; the published files
are downloaded unmodified and described.

```bash
nextflow run . -profile processed --acquisition_mode processed_only
```

This writes `results/processed/acquisition/manifest.json` (identifiers, URLs,
retrieval dates, checksums, versions) and `results/processed/qc/inventory.md`
(every real table header, dtype, missingness and example values).

**2. Approve the column mapping.** Read `inventory.md`, complete
`config/metadata_mapping.yaml`, and set `approved: true` with your name and
the inventory checksum you reviewed. Until you do, harmonization stops with an
explanatory report. This is the single most consequential step in the project.

**3. Run the analysis.**

```bash
nextflow run . -profile processed
```

**4. Read the QC report before any result.**
`results/processed/qc/qc_report.md` states whether the primary habitat
analysis was allowed to run at all.

---

## Execution profiles

| Profile | Data | Computation | Purpose |
| --- | --- | --- | --- |
| `test` | Deterministic synthetic | Very small | Unit and integration testing |
| `laptop` | Small verified public subset | Limited chains and viral panel | Exercise metadata and modelling |
| `processed` | Full published processed dataset | Moderate | Main analysis |
| `full` | Full published processed dataset | SLURM/HPC | Main analysis at full sampling, M0–M5 |
| `evolution` | Selected coding lineages | Moderate/HPC | Molecular evolutionary analyses |
| `manuscript` | Existing validated results | Low | Regenerate figures, tables, legends |

`full` is `processed` with a SLURM executor, the complete model ladder
including M5, and a larger sampling budget (4 chains × 2000 draws, 1000
permutations). Set `EVONEXUS_QUEUE` and, if your site needs them,
`EVONEXUS_CLUSTER_OPTS` before submitting, and give Nextflow a work directory
with quota:

```bash
EVONEXUS_QUEUE=your_partition nextflow run . -profile full,conda -w /scratch/$USER/evonexus-work -resume
```

M5 additionally needs `--host_newick`, a Newick tree over the host species
names as they appear in `hosts.host_species_resolved`. Without it M5 is
skipped and the connectivity-versus-relatedness panel of Figure 3 is marked
as not generated.

Combine with `docker`, `singularity`, `conda` or `hpc`. The container profiles
expect an image built from the repository `Dockerfile`; **no image is
published**, so build and push it first, or point `--container` at your own
tag:

```bash
docker build -t ghcr.io/aglucaci/evonexus:0.1.0 . && docker push ghcr.io/aglucaci/evonexus:0.1.0
```

### What `raw_reads` does and does not do

`--acquisition_mode raw_reads` retrieves FASTQ from SRA with per-run
checksums, and **nothing downstream consumes it**. EvoNexus has no
read-mapping, assembly or detection-reconstruction stage: assessability and
the library × lineage response both come from the published deposit's own
screening metadata. A `raw_reads` run therefore produces the same results as
`processed_only` after a far larger download, so the pipeline refuses the mode
unless you pass `--allow_unused_read_download true` to stage the files for
your own use.

---

## Statistical model

For library *l* and viral lineage unit *v*:

$$Y_{lv} \sim \mathrm{Bernoulli}(p_{lv})$$

$$\mathrm{logit}(p_{lv}) = \alpha + a_{f(v)} + a_v + \beta_U U_l + b_{h(l)}
+ \kappa_{h(l),f(v)} + c_{t(l)} + g_{r(l)}
+ \gamma_D \log D_l + \gamma_N \log N_l + \gamma_T T_l + \boldsymbol{\theta}^\top E_l$$

Viral family and lineage effects, host and tissue effects, a
shrinkage-regularized host-by-family interaction, site and batch effects,
sequencing depth, pool size, time, and supported ecological covariates. The
host-by-family interaction matters: without it, host-sharing predictions are
dominated by overall differences in host and lineage detection frequency
rather than by lineage-specific host associations.

| Model | Adds | Purpose |
| --- | --- | --- |
| M0 | Viral identity + technical covariates | Detection baseline |
| M1 | Host and tissue | Sample composition |
| M2 | Site | Geographic confounding |
| M3 | Habitat + host × family | **Primary habitat comparison** |
| M4 | Supported ecology | Ecological explanation |
| M5 | Host phylogenetic covariance | Evolutionary host relatedness |

**Primary estimand:** the difference in average lineage-level detection
probability between settlement and wild habitats, standardized to the same
overlap-supported host and sampling distribution, at a *new* site and batch
drawn from their estimated population distributions.

### Samplers

PyMC with NUTS is the intended engine for published results. If PyMC is not
installed, EvoNexus falls back to a self-contained Pólya-Gamma Gibbs sampler
with interweaved centered/non-centered scale updates, so the whole pipeline
runs on a bare laptop and in CI. The engine is recorded in every artifact and
named in every legend — a fallback fit can never be mistaken for a NUTS fit.
The same R-hat, ESS and posterior-predictive gates apply to both.

---

## The seven data gates

A critical failure stops the primary habitat analysis and emits an explanatory
QC report. Descriptive outputs still run, labelled as descriptive-only.

| Gate | Condition | Critical |
| --- | --- | --- |
| G1 | Every detection maps to a valid library and catalogued lineage | yes |
| G2 | Host species and pool relationships are documented | no |
| G3 | Positives and assay-tested negatives are distinguishable | yes |
| G4 | Habitat is known at a supported observational level | yes |
| G5 | Host comparisons across habitats are estimable | yes |
| G6 | Repeated libraries and shared pools are representable | no |
| G7 | Eligible mammal-associated RNA viral lineages are curated | yes |

---

## Outputs and the manuscript

| Figure | File | Content |
| --- | --- | --- |
| 1 | `results/*/figures/Figure1.*` | Sampling structure, host × habitat overlap, effort, exclusions |
| 2 | `results/*/figures/Figure2.*` | Raw rates, standardized predictions, posterior contrast, sensitivity |
| 3 | `results/*/figures/Figure3.*` | Host–virus structure, connectivity, nulls, host relatedness |
| 4 | `results/*/figures/Figure4.*` | Phylogenies, divergence, selection, host breadth |
| 5 | `results/*/figures/Figure5.*` | Model comparison, calibration, posterior checks, held-out scores |

Each figure writes its source tables and a `FigureN_legend.md` generated from
the analysis outputs. Tables 1–4 and Supplementary 1–5 are emitted as both
data files and formatted Markdown. **If an analysis is unavailable, the figure,
table and legend say "NOT GENERATED" and name what is missing.** No sample
size, interval or conclusion is ever fabricated.

---

## Repository layout

```
main.nf  nextflow.config  conf/      Nextflow entry point and profiles
workflows/  modules/                 DSL2 subworkflows and processes
scripts/evonexus/                    The Python package (all analysis logic)
config/                              Datasets, mapping, filters, priors, gates
docs/                                Data sources, dictionary, methods, reproducibility
tests/                               Unit and integration tests
results/                             Generated outputs (not committed)
```

Large raw data, read-mapping intermediates and MCMC traces are not committed.
The repository keeps manifests, checksums, code, configuration and lightweight
publication outputs.

## Documentation

- [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) — accessions, publications, licensing
- [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md) — canonical schema (generated)
- [`docs/STATISTICAL_METHODS.md`](docs/STATISTICAL_METHODS.md) — likelihood, priors, estimands, diagnostics
- [`docs/WORKFLOW.md`](docs/WORKFLOW.md) — DAG, processes, parameters, profiles
- [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) — versions, seeds, checksums, rerunning a release

## Testing

```bash
PYTHONPATH=scripts python -m pytest tests -q
nextflow run . -profile test
```

CI runs the unit tests, the metadata-schema tests and the small deterministic
integration test. Full-data analyses and long MCMC runs are not CI jobs.

## Citation

See [`CITATION.cff`](CITATION.cff). Cite the underlying data sources
independently — see [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md).

## License

MIT for the pipeline code. The public datasets carry their own terms; verify
them in the source deposits before redistribution.
