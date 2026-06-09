#!/bin/bash

# Purpose: cluster per-dataset PFMs with gimme motifs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${SCRIPT_DIR}/../config.sh"

slurm_job_active() {
  local job_name="$1"
  squeue -h -u "${USER}" -n "${job_name}" -t PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED 2>/dev/null | grep -q .
}

pfm_has_motifs() {
  local pfm_path="$1"
  [[ -s "${pfm_path}" ]] && grep -q '^>' "${pfm_path}"
}

compendium_dir="${compendium_dir:-${base_dir}/02-compendium}"
pfm_dir="${pfm_dir:-${compendium_dir}/pfm}"
gimme_cluster_dir="${gimme_cluster_dir:-${compendium_dir}/gimme_cluster}"

dataset_filter_regex="${CHROMBPNET_DATASET_FILTER_REGEX:-}"
pattern_types=(pos_patterns neg_patterns)

datasets=$(
  find "${pfm_dir}" -maxdepth 1 -type f -name '*.counts.pos_patterns.pfm' -print |
    while IFS= read -r path; do
      basename "${path}" '.counts.pos_patterns.pfm'
    done |
    sort -u
)

if [[ -n "${dataset_filter_regex}" ]]; then
  datasets=$(printf '%s\n' ${datasets} | grep -E "${dataset_filter_regex}" || true)
fi

for dataset in ${datasets}; do
  for pattern_type in "${pattern_types[@]}"; do
    key="${dataset}.counts.${pattern_type}"
    t=0.8
    ncpus="16"
    params="t${t}_n${ncpus}"
    input="${pfm_dir%/}/${key}.pfm"
    out_dir="${gimme_cluster_dir%/}/${key}_${params}"

    [[ -f "${input}" ]] || { echo "@ missing ${input}, skipping."; continue; }
    [[ -d "${out_dir}" ]] || mkdir -p "${out_dir}"

    job_name="02-gimme_cluster_${key}_${params}"
    if [[ -f "${out_dir}/clustered_motifs.pfm" ]]; then
      echo "@ found clustered motifs for ${key}, skipping."
    elif ! pfm_has_motifs "${input}"; then
      echo "@ ${key} has no motifs in ${input}; writing empty clustered output and skipping gimme."
      : > "${out_dir}/clustered_motifs.pfm"
    elif slurm_job_active "${job_name}"; then
      echo "@ ${job_name} already queued or running, skipping."
    else
      echo "@ submitting: gimme cluster ${input} ${out_dir} ${t}"
      sbatch -J "${job_name}" "${SCRIPT_DIR}/02a-jobscript.sh" "${input}" "${out_dir}" "${t}"
      sleep 2s
    fi
  done
done
