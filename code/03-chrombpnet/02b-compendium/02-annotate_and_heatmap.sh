#!/bin/bash

# Purpose:
# Post-process the CREsted-style compendium outputs to add motif annotations
# and render a CREsted-style clustermap from the merged pattern matrix.

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

resolve_crested_repo() {
  local candidate
  if [[ -n "${CRESTED_REPO:-}" && -d "${CRESTED_REPO}" ]]; then
    export CRESTED_REPO
    return 0
  fi

  for candidate in \
    "/scratch1/${USER}/code/CREsted" \
    "${HOME}/code/CREsted" \
    "${HOME}/projects/CREsted"
  do
    if [[ -d "${candidate}/src/crested" ]]; then
      export CRESTED_REPO="${candidate}"
      return 0
    fi
  done

  return 1
}

PRE_MODULES="${PRE_MODULES:-legacy/CentOS7 gcc/8.3.0}"
load_requested_modules "${PRE_MODULES}"

BIAS_PARAMS="${BIAS_PARAMS:-${chrombpnet_bias_params}}"
DEFAULT_OUTPUT_DIR="${base_dir}/02b-compendium/crested_patterns"
OUTPUT_DIR="${CRESTED_COMPENDIUM_OUT_DIR:-${crested_compendium_dir:-${DEFAULT_OUTPUT_DIR}}}"
ANNOTATION_DIR="${CRESTED_COMPENDIUM_ANNOTATION_DIR:-${OUTPUT_DIR}/annotation}"
PLOTS_DIR="${CRESTED_COMPENDIUM_PLOTS_DIR:-${OUTPUT_DIR}/plots}"
IMPORTANCE_THRESHOLD="${CRESTED_HEATMAP_IMPORTANCE_THRESHOLD:-0}"
HEATMAP_WIDTH="${CRESTED_HEATMAP_WIDTH:-25}"
HEATMAP_HEIGHT="${CRESTED_HEATMAP_HEIGHT:-8}"
ANNOTATION_PVAL_THRESHOLD="${CRESTED_ANNOTATION_PVAL_THRESHOLD:-0.05}"

mkdir -p "${ANNOTATION_DIR}" "${PLOTS_DIR}"

activate_python
resolve_crested_repo || true

echo "[$(timestamp)] ===== CREsted-style Annotation + Heatmap ====="
echo "[$(timestamp)] modisco_root=${modisco_dir%/}/bias_${BIAS_PARAMS}"
echo "[$(timestamp)] compendium_dir=${OUTPUT_DIR}"
echo "[$(timestamp)] annotation_dir=${ANNOTATION_DIR}"
echo "[$(timestamp)] plots_dir=${PLOTS_DIR}"
echo "[$(timestamp)] crested_repo=${CRESTED_REPO:-<unset>}"

"${PYTHON_BIN}" "${SCRIPT_DIR}/02-annotate_and_heatmap.py" \
  --compendium-dir "${OUTPUT_DIR}" \
  --modisco-root "${modisco_dir%/}/bias_${BIAS_PARAMS}" \
  --annotation-dir "${ANNOTATION_DIR}" \
  --plots-dir "${PLOTS_DIR}" \
  --importance-threshold "${IMPORTANCE_THRESHOLD}" \
  --heatmap-width "${HEATMAP_WIDTH}" \
  --heatmap-height "${HEATMAP_HEIGHT}" \
  --annotation-pval-threshold "${ANNOTATION_PVAL_THRESHOLD}"
