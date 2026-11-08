#!/usr/bin/env bash
#
# Run the whole EvoNexus pipeline on deterministic synthetic data, without
# Nextflow. Used by CI and useful for debugging a single step locally.
#
#   bin/run_synthetic_pipeline.sh [outdir] [libraries] [viruses]
#
# Outputs are watermarked and labelled synthetic. They validate software only.

set -euo pipefail

OUT="${1:-results/synthetic}"
LIBRARIES="${2:-90}"
VIRUSES="${3:-24}"
DRAWS="${EVONEXUS_DRAWS:-300}"
TUNE="${EVONEXUS_TUNE:-300}"
CHAINS="${EVONEXUS_CHAINS:-2}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="${ROOT}/scripts${PYTHONPATH:+:${PYTHONPATH}}"
export EVONEXUS_CONFIG_DIR="${ROOT}/config"
export MPLBACKEND=Agg

evonexus() { python -m evonexus.cli "$@"; }

echo "=== EvoNexus synthetic pipeline -> ${OUT} ==="
rm -rf "${OUT}"
mkdir -p "${OUT}"

echo "--- 01 acquire (synthetic, offline) ---"
evonexus acquire --out-dir "${OUT}/acquisition" --synthetic --offline \
    --libraries "${LIBRARIES}" --viruses "${VIRUSES}"

echo "--- 02 harmonize ---"
evonexus harmonize \
    --inventory "${OUT}/acquisition/inventory.json" \
    --raw-dir   "${OUT}/acquisition/raw" \
    --out-dir   "${OUT}/processed" \
    --skip-approval --allow-role-guess --synthetic

echo "--- 03 QC gates ---"
evonexus qc --processed-dir "${OUT}/processed" --out-dir "${OUT}/qc" --synthetic

PRIMARY_ALLOWED=$(python -c "import json;print(json.load(open('${OUT}/qc/qc_gates.json'))['primary_model_allowed'])")
echo "    primary model allowed: ${PRIMARY_ALLOWED}"

echo "--- 04 detection matrix ---"
evonexus matrix --processed-dir "${OUT}/processed" --out-dir "${OUT}/matrix" \
    --resolution lineage --synthetic

FRAME="${OUT}/matrix/analysis_frame_lineage"
TRUTH="${OUT}/acquisition/raw/synthetic/synthetic_truth.json"

echo "--- 05 Bayesian inference ---"
SKIP=""
if [ "${PRIMARY_ALLOWED}" != "True" ]; then
    echo "    a critical QC gate failed: running the descriptive branch only"
    SKIP="--skip-primary"
fi
evonexus bayesian \
    --frame "${FRAME}" --out-dir "${OUT}/bayesian" \
    --models M0 M3 --draws "${DRAWS}" --tune "${TUNE}" --chains "${CHAINS}" \
    --truth "${TRUTH}" --processed-dir "${OUT}/processed" \
    --no-sensitivity --synthetic ${SKIP}

if [ -z "${SKIP}" ]; then
    echo "--- 06 connectivity ---"
    evonexus connectivity \
        --posterior-dir "${OUT}/bayesian" --model M3 --frame "${FRAME}" \
        --out-dir "${OUT}/connectivity" --permutations 100 --synthetic

    echo "--- 07 validation ---"
    evonexus validate \
        --frame "${FRAME}" --out-dir "${OUT}/validation" \
        --posterior-dir "${OUT}/bayesian" --model M3 \
        --group-columns site_id --folds 3 \
        --draws 200 --tune 200 --chains 2 --synthetic
fi

echo "--- 08 evolution (no sequences supplied: expect a documented exclusion) ---"
evonexus evolution --processed-dir "${OUT}/processed" --out-dir "${OUT}/evolution" --synthetic

echo "--- 09 manuscript outputs ---"
evonexus report --results-dir "${OUT}" --out-dir "${OUT}"

echo
echo "=== done ==="
echo "QC report : ${OUT}/qc/qc_report.md"
echo "Run report: ${OUT}/manuscript/analysis_report.md"
echo "Figures   : ${OUT}/figures/"
echo
echo "SYNTHETIC DATA — these outputs validate software only."
