#!/bin/bash
#SBATCH --job-name=chrombpnet-compendium
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=Log/%x-%j.out
#SBATCH --error=Log/%x-%j.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
  set +u
  [[ -f /etc/profile ]] && source /etc/profile
  [[ -f /usr/share/lmod/lmod/init/bash ]] && source /usr/share/lmod/lmod/init/bash
  [[ -f /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash
  set -u
fi

source ~/.bashrc 2>/dev/null || true

activate_conda_env() {
  local env_name="$1"
  if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda is unavailable; cannot activate ${env_name}" >&2
    return 1
  fi
  eval "$(conda shell.bash hook)"
  conda activate "${env_name}"
}

RUN_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "${RUN_DIR}"

resolve_script_dir() {
  local candidate
  for candidate in \
    "${RUN_DIR}" \
    "${RUN_DIR}/code/03-chrombpnet/02-compendium" \
    "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" 2>/dev/null && pwd)"
  do
    if [[ -n "${candidate}" && -f "${candidate}/../config.sh" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

SCRIPT_DIR="$(resolve_script_dir)" || {
  echo "ERROR: could not locate 02-compendium directory from RUN_DIR=${RUN_DIR}"
  exit 1
}
cd "${SCRIPT_DIR}"
source "${SCRIPT_DIR}/../config.sh"

DATASET_FILTER_REGEX="${CHROMBPNET_DATASET_FILTER_REGEX:-}"
BIAS_PARAMS="${BIAS_PARAMS:-${chrombpnet_bias_params}}"
POLL_SECONDS="${POLL_SECONDS:-60}"
STOP_AFTER_STAGE="${STOP_AFTER_STAGE:-}"
compendium_dir="${compendium_dir:-${base_dir}/02-compendium}"
pfm_dir="${pfm_dir:-${compendium_dir}/pfm}"
gimme_cluster_dir="${gimme_cluster_dir:-${compendium_dir}/gimme_cluster}"
modisco_merged_dir="${modisco_merged_dir:-${compendium_dir}/modisco_merged}"
modisco_comp_dir="${modisco_comp_dir:-${compendium_dir}/modisco_compiled}"
modisco_comp_anno_dir="${modisco_comp_anno_dir:-${compendium_dir}/modisco_compiled_anno}"
hits_unified_dir="${hits_unified_dir:-${compendium_dir}/hits_unified}"
hits_reconciled_dir="${hits_reconciled_dir:-${compendium_dir}/hits_reconciled}"
nucleoatac_dir="${nucleoatac_dir:-${compendium_dir}/nucleoatac}"
vierstra_dir="${vierstra_dir:-${refs}/Vierstra_motifs}"
vierstra_meme_db="${vierstra_meme_db:-${vierstra_dir}/all.dbs.meme}"
compendium_annotation_tsv="${compendium_annotation_tsv:-${SCRIPT_DIR}/04d-ChromBPNet_de_novo_motifs.tsv}"
ANNOTATION_TSV="${COMPENDIUM_ANNOTATION_TSV:-${compendium_annotation_tsv}}"
VIERSTRA_MEME_DB="${VIERSTRA_MEME_DB:-${vierstra_meme_db}}"

export compendium_dir pfm_dir gimme_cluster_dir modisco_merged_dir modisco_comp_dir
export modisco_comp_anno_dir hits_unified_dir hits_reconciled_dir nucleoatac_dir
export vierstra_dir vierstra_meme_db compendium_annotation_tsv ANNOTATION_TSV VIERSTRA_MEME_DB

mkdir -p Log

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
    if [[ -n "${DATASET_FILTER_REGEX}" ]]; then
      grep -E "${DATASET_FILTER_REGEX}" || true
    else
      cat
    fi
}

wait_for_count() {
  local label="$1"
  local expected="$2"
  local cmd="$3"
  local current
  while true; do
    current=$(eval "$cmd")
    current=$(printf '%s' "${current}" | tail -n 1 | tr -d '[:space:]')
    echo "[$(timestamp)] ${label}: ${current}/${expected}"
    if [[ "${current}" == "${expected}" ]]; then
      break
    fi
    sleep "${POLL_SECONDS}"
  done
}

stop_if_requested() {
  local stage="$1"
  if [[ -n "${STOP_AFTER_STAGE}" && "${STOP_AFTER_STAGE}" == "${stage}" ]]; then
    echo "[$(timestamp)] STOP_AFTER_STAGE=${stage}; exiting after requested stage."
    exit 0
  fi
}

datasets_count() {
  list_target_datasets | awk 'NF {count++} END {print count+0}'
}

count_pfm_outputs() {
  local count=0
  local dataset
  local pattern_type
  while IFS= read -r dataset; do
    [[ -n "${dataset}" ]] || continue
    for pattern_type in pos_patterns neg_patterns; do
      if [[ -f "${pfm_dir%/}/${dataset}.counts.${pattern_type}.pfm" ]]; then
        count=$((count + 1))
      fi
    done
  done < <(list_target_datasets)
  printf '%s\n' "${count}"
}

count_gimme_outputs() {
  local count=0
  local dataset
  local pattern_type
  while IFS= read -r dataset; do
    [[ -n "${dataset}" ]] || continue
    for pattern_type in pos_patterns neg_patterns; do
      if [[ -f "${gimme_cluster_dir%/}/${dataset}.counts.${pattern_type}_t0.8_n16/clustered_motifs.pfm" ]]; then
        count=$((count + 1))
      fi
    done
  done < <(list_target_datasets)
  printf '%s\n' "${count}"
}

count_merge_done_files() {
  find "${modisco_merged_dir%/}" -maxdepth 1 -name '.*.done' -type f | awk 'END {print NR+0}'
}

count_hitcall_outputs() {
  local count=0
  local dataset
  while IFS= read -r dataset; do
    [[ -n "${dataset}" ]] || continue
    if [[ -f "${hits_unified_dir%/}/${dataset}/counts_v0.23_a0.8_all/hits.bed.gz" ]] && \
       [[ -f "${hits_unified_dir%/}/${dataset}/counts_v0.23_a0.8_nocompo/hits.bed.gz" ]]; then
      count=$((count + 1))
    fi
  done < <(list_target_datasets)
  printf '%s\n' "${count}"
}

count_nucleoatac_outputs() {
  local count=0
  local dataset
  while IFS= read -r dataset; do
    [[ -n "${dataset}" ]] || continue
    if [[ -f "${nucleoatac_dir%/}/${dataset}/${dataset}.nucmap_combined.bed.gz" ]]; then
      count=$((count + 1))
    fi
  done < <(list_target_datasets)
  printf '%s\n' "${count}"
}

echo "===== ChromBPNet Compendium Wrapper ====="
echo "Date: $(timestamp)"
echo "Run directory: ${RUN_DIR}"
echo "Script directory: ${SCRIPT_DIR}"
echo "Bias params: ${BIAS_PARAMS}"
echo "Dataset filter: ${DATASET_FILTER_REGEX:-<all datasets>}"
echo "Annotation TSV: ${ANNOTATION_TSV}"
echo "Vierstra MEME DB: ${VIERSTRA_MEME_DB}"

TARGET_DATASETS=$(datasets_count)
TARGET_DATASET_PATTERN_OUTPUTS=$(( TARGET_DATASETS * 2 ))

if [[ "${TARGET_DATASETS}" -eq 0 ]]; then
  echo "[$(timestamp)] ERROR: no datasets available from keep files or MoDISco directories"
  exit 1
fi

echo "[$(timestamp)] Step 01: MoDISco to PFM"
BIAS_PARAMS="${BIAS_PARAMS}" CHROMBPNET_DATASET_FILTER_REGEX="${DATASET_FILTER_REGEX}" bash 01-run_modisco_to_pfm_per_organ.sh
wait_for_count "PFM outputs" "${TARGET_DATASET_PATTERN_OUTPUTS}" "count_pfm_outputs"
stop_if_requested "01"

echo "[$(timestamp)] Step 02a: per-dataset gimme cluster"
CHROMBPNET_DATASET_FILTER_REGEX="${DATASET_FILTER_REGEX}" bash 02a-gimme_cluster_motifs.sh
wait_for_count "Per-dataset gimme outputs" "${TARGET_DATASET_PATTERN_OUTPUTS}" "count_gimme_outputs"
stop_if_requested "02a"

echo "[$(timestamp)] Step 02b: global gimme cluster"
bash 02b-gimme_cluster_clusters.sh
stop_if_requested "02b"

echo "[$(timestamp)] Step 03: merge MoDISco patterns"
BIAS_PARAMS="${BIAS_PARAMS}" bash 03-merge_modisco.sh
EXPECTED_BATCHES=$(cut -f 4 "${gimme_cluster_dir%/}/gimme_cluster_all_cluster_key.tsv" | sort -u | awk 'NF {count++} END {print count+0}')
wait_for_count "Merged MoDISco batches" "${EXPECTED_BATCHES}" "count_merge_done_files"
stop_if_requested "03"

echo "[$(timestamp)] Step 04a: compile merged MoDISco objects"
activate_conda_env modiscolite
python 04a-compile_modisco_obj.py \
  --modisco-merged-dir "${modisco_merged_dir}" \
  --cluster-key "${gimme_cluster_dir%/}/gimme_cluster_all_cluster_key.tsv" \
  --out-h5 "${modisco_comp_dir%/}/modisco_compiled.h5" \
  --out-tsv "${modisco_comp_dir%/}/modisco_compiled.tsv"
stop_if_requested "04a"

echo "[$(timestamp)] Step 04b: TOMTOM matching"
VIERSTRA_MEME_DB="${VIERSTRA_MEME_DB}" bash 04b-run.sh
stop_if_requested "04b"

echo "[$(timestamp)] Step 04f: export unified MEME"
bash 04f-CWM_to_MEME.sh
stop_if_requested "04f"

if [[ ! -f "${ANNOTATION_TSV}" ]]; then
  echo "[$(timestamp)] Manual motif annotation file not found: ${ANNOTATION_TSV}"
  echo "[$(timestamp)] Stopping after automated compendium build. Add the annotation TSV, then rerun from this wrapper."
  exit 0
fi

echo "[$(timestamp)] Step 04d/04g: finalize motif annotation products"
activate_conda_env modiscolite
python 04d-index_motifs.py --input "${ANNOTATION_TSV}"
MOTIF_ANNOTATION_FILE="${ANNOTATION_TSV}" MODISCO_COMP_ANNO_DIR="${modisco_comp_anno_dir}" Rscript 04g-prep_motif_names.R
stop_if_requested "04g"

echo "[$(timestamp)] Step 05: call hits"
BIAS_PARAMS="${BIAS_PARAMS}" CHROMBPNET_DATASET_FILTER_REGEX="${DATASET_FILTER_REGEX}" bash 05-call_hits.sh
wait_for_count "Unified hit calls" "${TARGET_DATASETS}" "count_hitcall_outputs"
stop_if_requested "05"

echo "[$(timestamp)] Step 06a: generate hitcall reports"
CHROMBPNET_DATASET_FILTER_REGEX="${DATASET_FILTER_REGEX}" bash 06a-generate_hitcall_reports.sh
stop_if_requested "06a"

echo "[$(timestamp)] Step 06b: reconcile hits"
CHROMBPNET_DATASET_FILTER_REGEX="${DATASET_FILTER_REGEX}" bash 06b-reconcile_hits.sh
stop_if_requested "06b"

echo "[$(timestamp)] Step 07: run NucleoATAC"
CHROMBPNET_DATASET_FILTER_REGEX="${DATASET_FILTER_REGEX}" bash 07-run_nucleoatac.sh
wait_for_count "NucleoATAC outputs" "${TARGET_DATASETS}" "count_nucleoatac_outputs"

echo "[$(timestamp)] Compendium wrapper completed."
