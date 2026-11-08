# Contributing

## Ground rules

EvoNexus is a research pipeline whose conclusions depend on its data handling
being correct. Three rules are not negotiable, and the tests enforce them:

1. **Unknown is not zero.** Habitat, host species and detection status each
   carry an explicit unknown/unassessable state. Never coerce one into a
   modelled category or a numeric zero.
2. **A missing analysis is reported as missing.** If an input is unavailable,
   the figure, table or legend says "NOT GENERATED" and names what is absent.
   Never substitute a placeholder value, a sample size or an interval.
3. **Nothing is presented as more certain than it is.** Synthetic results are
   watermarked and labelled. The habitat contrast carries its "not a causal
   effect" statement. The sampling engine is named wherever results appear.

A pull request that weakens any of these needs an explicit argument in the
description, not just passing tests.

## Development setup

```bash
conda env create -f environment.yml
conda activate evonexus
PYTHONPATH=scripts python -m pytest tests -q
nextflow run . -profile test
```

## Coding conventions

- Python 3.11+, type hints on public functions, `ruff` clean.
- Every analysis step takes explicit paths and writes a `RunRecord`.
- Configuration belongs in `config/*.yaml`, never hard-coded in a script.
- New canonical columns go in `evonexus.core.schema` first; the data
  dictionary is generated from it (`evonexus.cli schema --markdown`).
- Nextflow processes stay thin: one CLI subcommand, no analysis logic in
  Groovy.

## Tests

- `tests/test_metadata.py` — configuration, mapping resolution, schema,
  harmonization invariants.
- `tests/test_models.py` — design construction, sampling, diagnostics, the
  estimand, parameter recovery on a known generating process.
- `tests/test_figures.py` — QC gates, connectivity, validation, and the
  "NOT GENERATED" behaviour.

New statistical code needs a test that would fail if the statistic were
wrong — a recovery check or an invariant, not only a smoke test.

## Data handling

- Never commit raw data, FASTQ, BAM, posterior traces or anything with a
  restrictive license.
- Never commit a mapping marked `approved: true` for the real dataset; that
  approval is a decision each analyst makes after reading the inventory.
- When adding a data source, declare it in `config/datasets.yaml` with its
  identifier, license and tier. Do not fetch anything not declared there.

## Review

Pull requests need: passing CI, a CHANGELOG entry, and — for anything that
changes a number in a figure or table — a note in the description saying
which outputs change and why.
