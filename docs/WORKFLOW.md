# Workflow

## Pipeline DAG

```
                    ┌──────────────────────────┐
                    │ 01 ACQUIRE_DATA          │  manifest + raw files
                    │                          │  + inventory of real headers
                    └────────────┬─────────────┘
                                 ▼
                    ┌──────────────────────────┐
                    │ 02 HARMONIZE_METADATA    │  canonical tables + crosswalk
                    │    (blocks until the     │
                    │     mapping is approved) │
                    └────────────┬─────────────┘
                                 ▼
                                 ▼
                                      ┌────────────────────┐
                                      │ 04 CREATE_MATRIX   │
                                      │    library×lineage │
                                      └─────────┬──────────┘
                                                ▼
                                            ┌────────────────────┐
                                            │ 05 RUN_QC_GATES    │
                                            │    G1 … G7         │
                                            └─────────┬──────────┘
                                   pass ◄─────────────┴─────────────► fail
                                    ▼                                  ▼
                        ┌───────────────────────┐        ┌────────────────────────┐
                        │ 06 FIT_MODEL M0…M5    │        │ descriptive branch,    │
                        │    prior check,       │        │ labelled, no habitat   │
                        │    diagnostics,       │        │ effect estimated       │
                        │    contrast, sensitivity│      └───────────┬────────────┘
                        └───────┬───────────────┘                    │
             ┌──────────────────┼──────────────────┐                 │
             ▼                  ▼                  ▼                 │
┌─────────────────────┐ ┌───────────────┐ ┌─────────────────┐        │
│ 07 RUN_EVOLUTION    │ │ 09 CONNECTIV. │ │ 10 CROSS_VALID. │        │
│    align/tree/select│ │    + nulls    │ │    held-out site│        │
└──────────┬──────────┘ └───────┬───────┘ └────────┬────────┘        │
           └────────────────────┴──────────────────┴─────────────────┘
                                 ▼
                    ┌──────────────────────────┐
                    │ 11-12 GENERATE_REPORT    │  Figures 1-5, Tables, legends
                    └──────────────────────────┘
```

Render the real DAG with `nextflow run . -profile test -with-dag dag.html`.

The curated viral catalog has no process of its own: harmonization produces
`processed/viruses.*` with eligibility flags and exclusion reasons, and a
process that only copied it under another name would add a failure mode
without adding information.

### Publishing contract

Each process writes into a stage-named directory (`acquisition/`,
`processed/`, `qc/`, `matrix/`, `bayesian/`, …) and publishes to
`params.outdir` itself, so results land at `results/<profile>/<stage>/`.
Naming the stage in both the output path and `publishDir` would produce
`results/acquisition/acquisition/`.

### Nextflow version

`nextflow.config` declares `nextflowVersion '>=23.10.0'`, and CI runs the test
profile on both that minimum and the current stable release. Nextflow's strict
syntax (25.10 onwards) rejects constructs the older parser accepted —
`include` inside a workflow body, `return` in a workflow, and
`workflow.onComplete` as a script-level statement — so the pipeline is written
to satisfy both.

## Processes

| Process | Module | Label | Does |
| --- | --- | --- | --- |
| `ACQUIRE_DATA` | acquisition | network | Manifest, raw download, inventory |
| `DOWNLOAD_READS` | acquisition | network | SRA FASTQ (`raw_reads` only) |
| `HARMONIZE_METADATA` | metadata | light | Canonical tables + crosswalk |
| `RUN_QC_GATES` | metadata | light | G1–G7, QC report, waterfall |
| `CREATE_MATRIX` | viromics | light | Analysis frame per resolution |
| `FIT_MODEL` | statistics | heavy | M0–M5, contrast, sensitivity |
| `COMPUTE_CONNECTIVITY` | statistics | standard | C(i,j) and nulls |
| `CROSS_VALIDATE` | statistics | heavy | Held-out site/region scores |
| `BUILD_HOST_PHYLOGENY` | phylogenetics | light | Host covariance for M5 |
| `RUN_EVOLUTION` | phylogenetics | heavy | MAFFT, IQ-TREE 2, HyPhy |
| `GENERATE_REPORT` | reporting | light | Figures, tables, legends |
| `RENDER_DATA_DICTIONARY` | reporting | tiny | Generated schema docs |

## Parameters

See [`params.yaml`](../params.yaml) for a commented template and
[`nextflow.config`](../nextflow.config) for defaults.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `acquisition_mode` | `processed_only` | `processed_only` or `raw_reads` |
| `require_mapping_approval` | `true` | Harmonization stops until a human approves the mapping |
| `fail_on_critical_gate` | `false` | Skip the primary model vs. abort the run |
| `models` | `[M0..M4]` | Ladder members to fit |
| `primary_model` | `M3` | Model behind the primary estimand |
| `prior_variant` | `regularizing` | `regularizing`, `tight` or `wide` |
| `engine` | `auto` | `pymc_nuts` preferred when installed |
| `draws` / `tune` / `chains` / `seed` | 1000/1000/4/20250101 | Sampler settings |
| `resolutions` | `[lineage, family, strict]` | Lineage resolutions; the first drives modelling |
| `permutations` | `500` | Connectivity permutation null |
| `cv_group_columns` / `cv_folds` | `[site_id, region]` / 5 | Held-out grouping |
| `sequences` / `host_newick` | `null` | Optional evolutionary inputs |

## Running without Nextflow

Every process is a thin wrapper around one CLI subcommand, so the whole
pipeline can be driven directly — useful for debugging a single step:

```bash
export PYTHONPATH=scripts
python -m evonexus.cli acquire      --out-dir results/x/acquisition --synthetic --offline
python -m evonexus.cli harmonize    --inventory results/x/acquisition/inventory.json \
                                    --raw-dir  results/x/acquisition/raw \
                                    --out-dir  results/x/processed --skip-approval --synthetic
python -m evonexus.cli qc           --processed-dir results/x/processed --out-dir results/x/qc --synthetic
python -m evonexus.cli matrix       --processed-dir results/x/processed --out-dir results/x/matrix --synthetic
python -m evonexus.cli bayesian     --frame results/x/matrix/analysis_frame_lineage \
                                    --out-dir results/x/bayesian --models M0 M3 --synthetic
python -m evonexus.cli connectivity --posterior-dir results/x/bayesian --model M3 \
                                    --frame results/x/matrix/analysis_frame_lineage \
                                    --out-dir results/x/connectivity --synthetic
python -m evonexus.cli validate     --frame results/x/matrix/analysis_frame_lineage \
                                    --out-dir results/x/validation --posterior-dir results/x/bayesian --synthetic
python -m evonexus.cli evolution    --processed-dir results/x/processed --out-dir results/x/evolution --synthetic
python -m evonexus.cli report       --results-dir results/x
```

`python -m evonexus.cli --help` lists every subcommand.

## Resources

`conf/base.config` defines `tiny`, `light`, `standard`, `heavy` and `network`
labels; each profile overrides them. On SLURM, set `EVONEXUS_QUEUE` and
`EVONEXUS_CLUSTER_OPTS`. Network processes are rate-limited (`maxForks`) to
stay inside NCBI's unauthenticated limits.
