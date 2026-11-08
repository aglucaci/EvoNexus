# Reproducibility

## What is recorded, automatically

Every pipeline step writes a `run_record_<step>.json` next to its outputs:

- the step, its parameters, its status and its duration;
- every input and output path with size and SHA-256 checksum;
- every info/warning/error message the step emitted;
- the environment: Python version, platform, package versions, git commit,
  git branch and whether the working tree was dirty.

Acquisition additionally writes a permanent `manifest.json` (repository
identifiers, filenames, URLs, UTC retrieval dates, checksums, deposit
versions and licenses), and the reporting step writes
`manuscript/analysis_report.md`, which carries the provenance table, the gate
summary and an explicit list of anything that was not generated.

## Seeds

A single seed (`params.seed`, default `20250101`) drives the synthetic data
generator, MCMC chain initialization, the posterior standardization draws for
new sites and batches, the connectivity permutation null and the
cross-validation fold assignment. Chain *k* uses `seed + 1000k`, so chains are
independent but the whole run is deterministic.

Determinism is tested: `test_sampling_is_reproducible` asserts that two fits
with the same seed produce identical draws.

## Software versions

Pin with `environment.yml`. The exact versions used by a run are in every run
record and in the report:

```bash
PYTHONPATH=scripts python -m evonexus.cli engines
```

Containers (`-profile docker` or `-profile singularity`) pin external tools:
MAFFT 7.520, IQ-TREE 2.2.2.7, HyPhy 2.5.60.

## Sampling engines

`engine: auto` prefers PyMC/NUTS and falls back to the built-in Pólya-Gamma
Gibbs sampler. **The two engines do not produce bit-identical draws**, and the
engine is recorded in `bayesian_status.json`, in every posterior sidecar and
in every figure legend. For a release, pin `engine: pymc_nuts` so a rerun
cannot silently change engines.

## Reproducing a release

A release identifies five things. Record all of them:

| Item | Where |
| --- | --- |
| Source-data manifest | `results/<profile>/acquisition/manifest.json` |
| Code commit | `git_commit` in any run record |
| Configuration | `config/*.yaml` at that commit, plus the `params.yaml` used |
| Environment | `environment.yml` and the `packages` block in the run records |
| Posterior results | `results/<profile>/bayesian/posterior_*.npz` and `.json` |

To rerun:

```bash
git checkout <release-tag>
conda env create -f environment.yml && conda activate evonexus
nextflow run . -profile processed -params-file <the params.yaml from the release>
```

Then compare the new `manifest.json` checksums against the released one. A
checksum mismatch means the public deposit changed — investigate before
comparing results, because the inputs are not the same data.

## Resuming

Nextflow caches completed tasks:

```bash
nextflow run . -profile processed -resume
```

SRA downloads use `storeDir`, so they survive `work/` being cleaned.

## What is not committed

Raw data, read-mapping intermediates and MCMC traces are excluded by
`.gitignore`. The repository keeps acquisition manifests, checksums, scripts,
configuration and lightweight publication outputs. A release may additionally
attach the posterior `.npz` files as release assets.

## Known sources of non-determinism

- **Figshare and NCBI content can change.** The manifest's checksums are the
  only way to know. Run acquisition and diff the manifest before assuming two
  runs saw the same data.
- **BLAS threading** can reorder floating-point reductions. Set
  `OMP_NUM_THREADS=1` for bit-level reproducibility at some cost in speed.
- **External phylogenetic tools** (IQ-TREE, HyPhy) have their own seeds and
  heuristics; pin their versions through the containers.
- **Parquet versus CSV.** Without `pyarrow`, tables are written as `.csv.gz`.
  Values round-trip, but dtypes may differ subtly. Install `pyarrow` for
  published runs; the run record states which format was used.
