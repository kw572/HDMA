#!/bin/bash

# Purpose:
# Convert MoDISco CWMs to PFMs for downstream gimme clustering.
# Treat each dataset independently so the wrapper can operate without keep.tsv files.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${SCRIPT_DIR}/../config.sh"

eval "$(conda shell.bash hook)"
conda activate modiscolite

compendium_dir="${compendium_dir:-${base_dir}/02-compendium}"
pfm_dir="${pfm_dir:-${compendium_dir}/pfm}"

mkdir -p "${pfm_dir}"

bias_params="${BIAS_PARAMS:-${chrombpnet_bias_params}}"
dataset_filter_regex="${CHROMBPNET_DATASET_FILTER_REGEX:-}"

search_root="${modisco_scratch%/}/bias_${bias_params}"

if [[ ! -d "${search_root}" ]]; then
  echo "@ ERROR: modisco search root missing:"
  echo "${search_root}"
  exit 1
fi

datasets=$(
  find "${search_root}" \
    -mindepth 1 \
    -maxdepth 1 \
    -type d |
  while IFS= read -r dataset_path; do
    dataset="$(basename "${dataset_path}")"
    case "${dataset}" in
      NA|Log|logs|tmp|scratch|qc)
        continue
        ;;
    esac
    if [[ -d "${dataset_path}/counts_modisco_report" ]]; then
      printf '%s\n' "${dataset}"
    fi
  done |
  sort -u
)

if [[ -n "${dataset_filter_regex}" ]]; then
  datasets=$(
    printf '%s\n' ${datasets} |
    grep -E "${dataset_filter_regex}" || true
  )
fi

for dataset in ${datasets}; do
  echo "================================================================="
  echo "@ DATASET: ${dataset}"
  echo "================================================================="

  key1="${dataset}.counts.pos_patterns"
  out_file1="${pfm_dir%/}/${key1}.pfm"

  key2="${dataset}.counts.neg_patterns"
  out_file2="${pfm_dir%/}/${key2}.pfm"

  config="${pfm_dir%/}/${dataset}.counts.config.tsv"
  [[ -f "${config}" ]] && rm -f "${config}"

  modisco_counts_h5="${modisco_scratch%/}/bias_${bias_params}/${dataset}/counts_modisco_output.h5"

  if [[ ! -f "${modisco_counts_h5}" ]]; then
    echo "@ missing COUNTS modisco output:"
    echo "${modisco_counts_h5}"
    continue
  fi

  echo -e "${dataset}\t${modisco_counts_h5}" >> "${config}"
  cat "${config}"

  echo "@ generating POS PFMs"
  python 01-modisco_to_pfm.py \
    -c "${config}" \
    -o "${out_file1}" \
    -p pos_patterns

  echo "@ generating NEG PFMs"
  python 01-modisco_to_pfm.py \
    -c "${config}" \
    -o "${out_file2}" \
    -p neg_patterns
done

echo "@ DONE"
