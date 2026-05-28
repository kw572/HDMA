#!/bin/bash

# Purpose:
# Run Fi-NeMo hit calling per dataset using the original per-dataset TF-MoDISco
# motifs, then attach 02b merged-compendium IDs and annotations onto the hits.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

list_target_datasets() {
  if [[ -f "${chrombpnet_models_keep2}" ]]; then
    awk 'NF && $1 != "NA" {print $1}' "${chrombpnet_models_keep2}"
  elif [[ -f "${chrombpnet_models_keep}" ]]; then
    awk 'NF && $1 != "NA" {print $1}' "${chrombpnet_models_keep}"
  else
    find "${modisco_dir%/}/bias_${BIAS_PARAMS}" -mindepth 1 -maxdepth 1 -type d -print |
      while IFS= read -r path; do
        basename "${path}"
      done |
      awk 'NF && $1 != "NA" {print $1}'
  fi |
    sort -u |
    if [[ -n "${CHROMBPNET_DATASET_FILTER_REGEX:-}" ]]; then
      grep -E "${CHROMBPNET_DATASET_FILTER_REGEX}" || true
    else
      cat
    fi
}

BIAS_PARAMS="${BIAS_PARAMS:-${chrombpnet_bias_params}}"
COMPENDIUM_DIR="${CRESTED_COMPENDIUM_OUT_DIR:-}"
if [[ -z "${COMPENDIUM_DIR}" ]]; then
  COMPENDIUM_DIR="${base_dir%/}/02b-compendium/crested_patterns"
fi
DEFAULT_FINEMO_BASE="$(cd "${COMPENDIUM_DIR%/}/.." && pwd)"
OUTPUT_DIR="${CRESTED_FINEMO_OUT_DIR:-${DEFAULT_FINEMO_BASE}/finemo_hits}"
FINEMO_ALPHA="${CRESTED_FINEMO_ALPHA:-0.8}"
JOBSCRIPT="${SCRIPT_DIR}/03-call_finemo_hits_jobscript.sh"

mkdir -p "${OUTPUT_DIR}"

echo "[$(timestamp)] ===== 02b Fi-NeMo Hit Calling ====="
echo "[$(timestamp)] modisco_root=${modisco_dir%/}/bias_${BIAS_PARAMS}"
echo "[$(timestamp)] compendium_dir=${COMPENDIUM_DIR}"
echo "[$(timestamp)] output_dir=${OUTPUT_DIR}"
echo "[$(timestamp)] alpha=${FINEMO_ALPHA}"
echo "[$(timestamp)] dataset_filter=${CHROMBPNET_DATASET_FILTER_REGEX:-<all datasets>}"

datasets=$(list_target_datasets)

for dataset in ${datasets}; do
  peaks_bed="${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak"
  modisco_h5="${modisco_dir%/}/bias_${BIAS_PARAMS}/${dataset}/counts_modisco_output.h5"
  counts_shaps_h5="${contribs_dir%/}/bias_${BIAS_PARAMS}/${dataset}/average_shaps.counts.h5"
  dataset_out="${OUTPUT_DIR%/}/${dataset}"

  if [[ ! -f "${peaks_bed}" ]]; then
    echo "[$(timestamp)] WARNING: missing peaks for ${dataset}: ${peaks_bed}"
    continue
  fi
  if [[ ! -f "${modisco_h5}" ]]; then
    echo "[$(timestamp)] WARNING: missing modisco H5 for ${dataset}: ${modisco_h5}"
    continue
  fi
  if [[ ! -f "${counts_shaps_h5}" ]]; then
    echo "[$(timestamp)] WARNING: missing SHAP H5 for ${dataset}: ${counts_shaps_h5}"
    continue
  fi

  mkdir -p "${dataset_out}"

  if [[ -f "${dataset_out}/hits.compendium.tsv.gz" ]]; then
    echo "[$(timestamp)] ${dataset}: found reconciled hits, skipping"
    continue
  fi

  echo "[$(timestamp)] ${dataset}: submitting Fi-NeMo hit-calling job"
  sbatch -J "02b-finemo_${dataset}" "${JOBSCRIPT}" \
    "${dataset}" \
    "${peaks_bed}" \
    "${modisco_h5}" \
    "${counts_shaps_h5}" \
    "${dataset_out}" \
    "${COMPENDIUM_DIR}" \
    "${FINEMO_ALPHA}"

  sleep 2s
done
