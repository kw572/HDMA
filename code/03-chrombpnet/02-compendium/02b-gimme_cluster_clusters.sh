#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/02/02/02b-gimme_cluster_all_%j.out
#SBATCH --partition=main
#SBATCH -t 02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G

# Purpose: run a second-pass gimme clustering on the per-dataset motif clusters.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${SCRIPT_DIR}/../config.sh"

eval "$(conda shell.bash hook)"
conda activate meme-suite

BATCH_SIZE=5
compendium_dir="${compendium_dir:-${base_dir}/02-compendium}"
pfm_dir="${pfm_dir:-${compendium_dir}/pfm}"
gimme_cluster_dir="${gimme_cluster_dir:-${compendium_dir}/gimme_cluster}"
UNIFIED_CLUSTER_KEY="${gimme_cluster_dir%/}/gimme_cluster_all_cluster_key.tsv"
dataset_filter_regex="${CHROMBPNET_DATASET_FILTER_REGEX:-}"
params="t0.8_n16"
AVAILABLE_PATTERN_CLASSES=()

cluster_key_path_for_dir() {
  local out_dir="$1"
  local candidate
  for candidate in "${out_dir}/cluster_key.txt" "${out_dir}/cluster_key.tsv"; do
    if [[ -f "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

append_unified_cluster_key() {
  local key_path="$1"
  local pattern_class="$2"
  local prefix="$3"
  local start_idx="$4"
  local idx="${start_idx}"
  local cluster_id
  local component_patterns
  local n_component_patterns
  local batch

  while IFS=$'\t' read -r cluster_id component_patterns _; do
    [[ -n "${cluster_id}" && -n "${component_patterns}" ]] || continue
    n_component_patterns=$(awk -F',' '{print NF+0}' <<< "${component_patterns}")
    batch=$(( (idx - 1) / BATCH_SIZE + 1 ))
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${prefix}.${cluster_id}" \
      "${pattern_class}" \
      "${component_patterns}" \
      "${batch}" \
      "${n_component_patterns}" \
      "${n_component_patterns}" >> "${UNIFIED_CLUSTER_KEY}"
    idx=$((idx + 1))
  done < "${key_path}"

  printf '%s\n' "${idx}"
}

list_dataset_keys() {
  find "${gimme_cluster_dir}" -maxdepth 1 -mindepth 1 -type d -name "*.counts.pos_patterns_${params}" -print |
    while IFS= read -r path; do
      basename "${path}" "_${params}"
    done |
    sed 's/\.counts\.pos_patterns$//' |
    sort -u |
    if [[ -n "${dataset_filter_regex}" ]]; then
      grep -E "${dataset_filter_regex}" || true
    else
      cat
    fi
}

build_all_pfm() {
  local pattern_class="$1"
  local all_pfm="$2"
  local dataset
  local raw_pfm

  [[ -f "${all_pfm}" ]] && rm -f "${all_pfm}"

  while IFS= read -r dataset; do
    [[ -n "${dataset}" ]] || continue
    raw_pfm="${pfm_dir%/}/${dataset}.counts.${pattern_class}.pfm"
    if [[ -f "${raw_pfm}" ]]; then
      # Use original per-dataset motif names here so downstream merge jobs can
      # map directly back to MoDISco patterns without synthetic Average_* ids.
      sed "s/>/>${dataset}__/g" "${raw_pfm}" >> "${all_pfm}"
    fi
  done < <(list_dataset_keys)
}

pfm_has_motifs() {
  local pfm_path="$1"
  [[ -s "${pfm_path}" ]] && grep -q '^>' "${pfm_path}"
}

count_pfm_motifs() {
  local pfm_path="$1"
  if [[ -f "${pfm_path}" ]]; then
    awk '/^>/{count++} END {print count+0}' "${pfm_path}"
  else
    printf '0\n'
  fi
}

run_second_pass_cluster() {
  local pattern_class="$1"
  local all_pfm="$2"
  local key="all_clustered_clusters.counts.${pattern_class}"
  local out_dir="${gimme_cluster_dir%/}/${key}_${params}"

  if ! pfm_has_motifs "${all_pfm}"; then
    echo "@ ${pattern_class} has no motifs after per-dataset clustering; skipping second-pass clustering."
    return 1
  fi

  [[ -d "${out_dir}" ]] || mkdir -p "${out_dir}"

  echo "@ input: ${all_pfm}"
  echo "@ output: ${out_dir}"
  echo "@ running: gimme cluster ${all_pfm} ${out_dir} 0.8"

  gimme cluster "${all_pfm}" "${out_dir}" -t 0.8 -N 16
}

all_pos_pfm="${pfm_dir%/}/all_clustered_clusters.counts.pos_patterns.pfm"
build_all_pfm "pos_patterns" "${all_pos_pfm}"
count_pfm_motifs "${all_pos_pfm}"
if run_second_pass_cluster "pos_patterns" "${all_pos_pfm}"; then
  AVAILABLE_PATTERN_CLASSES+=("pos_patterns")
fi

all_neg_pfm="${pfm_dir%/}/all_clustered_clusters.counts.neg_patterns.pfm"
build_all_pfm "neg_patterns" "${all_neg_pfm}"
count_pfm_motifs "${all_neg_pfm}"
if run_second_pass_cluster "neg_patterns" "${all_neg_pfm}"; then
  AVAILABLE_PATTERN_CLASSES+=("neg_patterns")
fi

rm -f "${UNIFIED_CLUSTER_KEY}"

next_idx=1
for pattern_class in "${AVAILABLE_PATTERN_CLASSES[@]}"; do
  key="all_clustered_clusters.counts.${pattern_class}"
  out_dir="${gimme_cluster_dir%/}/${key}_${params}"
  cluster_key_path="$(cluster_key_path_for_dir "${out_dir}")" || {
    echo "ERROR: no cluster key file found under ${out_dir}"
    exit 1
  }
  if [[ "${pattern_class}" == "pos_patterns" ]]; then
    prefix="pos"
  else
    prefix="neg"
  fi
  next_idx="$(append_unified_cluster_key "${cluster_key_path}" "${pattern_class}" "${prefix}" "${next_idx}")"
done

echo "@ wrote unified cluster key: ${UNIFIED_CLUSTER_KEY}"
