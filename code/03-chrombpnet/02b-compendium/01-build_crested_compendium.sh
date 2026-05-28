#!/bin/bash

# Purpose:
# Build a CREsted-style motif compendium directly from the per-dataset
# TF-MoDISco outputs produced by 01-train_models.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

load_optional_module() {
  local mod="$1"
  [[ -n "${mod}" ]] || return 0
  if module -t avail "${mod}" 2>&1 | grep -Fq "${mod}"; then
    module load "${mod}"
  else
    echo "[$(timestamp)] WARNING: module '${mod}' is unavailable; continuing without it."
  fi
}

load_requested_modules() {
  local modules_string="$1"
  local mod
  [[ -n "${modules_string}" ]] || return 0
  read -r -a modules <<< "${modules_string}"
  for mod in "${modules[@]}"; do
    load_optional_module "${mod}"
  done
}

activate_python() {
  if [[ -n "${CRESTED_PYTHON_BIN:-}" ]]; then
    PYTHON_BIN="${CRESTED_PYTHON_BIN}"
    return 0
  fi

  local venv_path="${CRESTED_VENV:-}"
  if [[ -n "${venv_path}" && -x "${venv_path}/bin/python" ]]; then
    # shellcheck disable=SC1090
    source "${venv_path}/bin/activate"
    PYTHON_BIN="python"
    return 0
  fi

  eval "$(conda shell.bash hook)"
  conda activate "${CRESTED_ENV_NAME:-modiscolite}"
  PYTHON_BIN="python"
}

PRE_MODULES="${PRE_MODULES:-legacy/CentOS7 gcc/8.3.0}"
load_requested_modules "${PRE_MODULES}"

BIAS_PARAMS="${BIAS_PARAMS:-${chrombpnet_bias_params}}"
DATASET_FILTER_REGEX="${CHROMBPNET_DATASET_FILTER_REGEX:-}"
DEFAULT_OUTPUT_DIR="${base_dir}/02b-compendium/crested_patterns"
OUTPUT_DIR="${CRESTED_COMPENDIUM_OUT_DIR:-${crested_compendium_dir:-${DEFAULT_OUTPUT_DIR}}}"
SIM_THRESHOLD="${CRESTED_SIM_THRESHOLD:-6.0}"
TRIM_IC_THRESHOLD="${CRESTED_TRIM_IC_THRESHOLD:-0.025}"
DISCARD_IC_THRESHOLD="${CRESTED_DISCARD_IC_THRESHOLD:-0.1}"
PATTERN_PARAMETER="${CRESTED_PATTERN_PARAMETER:-seqlet_count}"
NORMALIZE_PATTERN_MATRIX="${CRESTED_NORMALIZE_PATTERN_MATRIX:-0}"
WRITE_SIMILARITY_MATRIX="${CRESTED_WRITE_SIMILARITY_MATRIX:-0}"

mkdir -p "${OUTPUT_DIR}"

activate_python

echo "[$(timestamp)] ===== CREsted-style Compendium ====="
echo "[$(timestamp)] modisco_root=${modisco_dir%/}/bias_${BIAS_PARAMS}"
echo "[$(timestamp)] output_dir=${OUTPUT_DIR}"
echo "[$(timestamp)] dataset_filter=${DATASET_FILTER_REGEX:-<all datasets>}"
echo "[$(timestamp)] sim_threshold=${SIM_THRESHOLD}"
echo "[$(timestamp)] trim_ic_threshold=${TRIM_IC_THRESHOLD}"
echo "[$(timestamp)] discard_ic_threshold=${DISCARD_IC_THRESHOLD}"
echo "[$(timestamp)] pattern_parameter=${PATTERN_PARAMETER}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/01-build_crested_compendium.py" \
  --modisco-root "${modisco_dir%/}/bias_${BIAS_PARAMS}" \
  --keep-file "${chrombpnet_models_keep2}" \
  --keep-file-fallback "${chrombpnet_models_keep}" \
  --output-dir "${OUTPUT_DIR}" \
  --dataset-filter-regex "${DATASET_FILTER_REGEX}" \
  --sim-threshold "${SIM_THRESHOLD}" \
  --trim-ic-threshold "${TRIM_IC_THRESHOLD}" \
  --discard-ic-threshold "${DISCARD_IC_THRESHOLD}" \
  --pattern-parameter "${PATTERN_PARAMETER}" \
  --normalize-pattern-matrix "${NORMALIZE_PATTERN_MATRIX}" \
  --write-similarity-matrix "${WRITE_SIMILARITY_MATRIX}"
