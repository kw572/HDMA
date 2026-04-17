#!/bin/bash
#SBATCH --job-name=chrombpnet-all
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=Log/%x-%j.out
#SBATCH --error=Log/%x-%j.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
  set +u
  if [[ -f /etc/profile ]]; then
    source /etc/profile
  fi
  if [[ -f /usr/share/lmod/lmod/init/bash ]]; then
    source /usr/share/lmod/lmod/init/bash
  elif [[ -f /usr/share/Modules/init/bash ]]; then
    source /usr/share/Modules/init/bash
  fi
  set -u
fi

source ~/.bashrc 2>/dev/null || true

RUN_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$RUN_DIR"

resolve_script_dir() {
  local candidate
  for candidate in \
    "$RUN_DIR" \
    "$RUN_DIR/code/03-chrombpnet/01-train_models" \
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
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] ERROR: could not locate 01-train_models directory from RUN_DIR=${RUN_DIR}"
  exit 1
}
cd "$SCRIPT_DIR"

source "$SCRIPT_DIR/../config.sh"

BIAS_CLUSTER="${BIAS_CLUSTER:-1-col_aspn_ogna}"
BIAS_MIN_THRESH="${BIAS_MIN_THRESH:-0.4}"
BIAS_FOLD="${BIAS_FOLD:-fold_0}"
BIAS_PARAMS="${BIAS_PARAMS:-${BIAS_CLUSTER}_thresh${BIAS_MIN_THRESH}}"
DATA_TYPE="${DATA_TYPE:-ATAC}"
MEME_DB="${MEME_DB:-${refs}/JASPAR_CORE_2016_vertebrates.meme}"
R_MODULE="${R_MODULE:-rstats/4.5.1}"
POLL_SECONDS="${POLL_SECONDS:-60}"
DATASET_FILTER_REGEX="${CHROMBPNET_DATASET_FILTER_REGEX:-}"
BIAS_SBATCH_ARGS="${BIAS_SBATCH_ARGS:---partition=gpu}"
FORCE_QC="${CHROMBPNET_FORCE_QC:-0}"

mkdir -p Log

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

activate_chrombpnet() {
  eval "$(conda shell.bash hook)"
  conda activate chrombpnet
}

activate_modiscolite() {
  eval "$(conda shell.bash hook)"
  conda activate modiscolite
}

activate_r() {
  module purge >/dev/null 2>&1 || true
  module load "${R_MODULE}"
  if [[ -n "${R_LIBS_USER:-}" ]]; then
    :
  elif [[ -n "${RENV_SJ:-}" ]]; then
    export R_LIBS_USER="${RENV_SJ}"
  else
    local renv_file
    for renv_file in \
      "${SCRIPT_DIR}/../../../SJ_RENV_DIR.txt" \
      "${RUN_DIR}/SJ_RENV_DIR.txt" \
      "${RUN_DIR}/../../../SJ_RENV_DIR.txt"
    do
      if [[ -f "${renv_file}" ]]; then
        export R_LIBS_USER
        IFS= read -r R_LIBS_USER < "${renv_file}"
        break
      fi
    done
  fi
  export OPENBLAS_NUM_THREADS=1
  export OMP_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  echo "[$(timestamp)] R_LIBS_USER=${R_LIBS_USER:-<unset>}"
}

assert_file_exists() {
  local label="$1"
  local path="$2"
  if [[ ! -f "$path" ]]; then
    echo "[$(timestamp)] ERROR: ${label} not found: ${path}"
    exit 1
  fi
}

dataset_count_all() {
  find "${chrombpnet_peaks_dir}" -maxdepth 1 -type f -name '*__peaks_bpnet.narrowPeak' ! -name '._*' -print |
    while IFS= read -r path; do
      basename "${path}" __peaks_bpnet.narrowPeak
    done |
    if [[ -n "${DATASET_FILTER_REGEX}" ]]; then
      grep -E "${DATASET_FILTER_REGEX}" || true
    else
      cat
    fi |
    awk 'NF {count++} END {print count+0}'
}

dataset_count_keep() {
  if [[ -f "${chrombpnet_models_keep}" ]]; then
    if [[ -n "${DATASET_FILTER_REGEX}" ]]; then
      awk -v pattern="${DATASET_FILTER_REGEX}" 'NF && $1 ~ pattern {count++} END {print count+0}' "${chrombpnet_models_keep}"
    else
      awk 'NF {count++} END {print count+0}' "${chrombpnet_models_keep}"
    fi
  else
    dataset_count_all
  fi
}

dataset_count_keep2() {
  if [[ -f "${chrombpnet_models_keep2}" ]]; then
    if [[ -n "${DATASET_FILTER_REGEX}" ]]; then
      awk -v pattern="${DATASET_FILTER_REGEX}" 'NF && $1 ~ pattern {count++} END {print count+0}' "${chrombpnet_models_keep2}"
    else
      awk 'NF {count++} END {print count+0}' "${chrombpnet_models_keep2}"
    fi
  else
    dataset_count_keep
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

echo "===== ChromBPNet Full Training Wrapper ====="
echo "Date: $(timestamp)"
echo "Run directory: ${RUN_DIR}"
echo "Script directory: ${SCRIPT_DIR}"
echo "Bias cluster: ${BIAS_CLUSTER}"
echo "Bias params: ${BIAS_PARAMS}"
echo "Data type: ${DATA_TYPE}"
echo "Dataset filter: ${DATASET_FILTER_REGEX:-<all datasets>}"
echo "MEME DB: ${MEME_DB}"
echo "R module: ${R_MODULE}"

assert_file_exists "Reference FASTA" "${ref_fasta}"
assert_file_exists "Chrom sizes" "${chromsizes}"
assert_file_exists "MEME database" "${MEME_DB}"

TOTAL_DATASETS=$(dataset_count_all)
EXPECTED_FOLDS=$(( TOTAL_DATASETS * 5 ))

if [[ "${TOTAL_DATASETS}" -eq 0 ]]; then
  echo "[$(timestamp)] ERROR: no ChromBPNet peak files found in ${chrombpnet_peaks_dir}"
  exit 1
fi

echo "[$(timestamp)] Found ${TOTAL_DATASETS} datasets and ${EXPECTED_FOLDS} expected fold jobs."

echo "[$(timestamp)] Step 01: bias model"
if [[ -f "${bias_dir%/}/${BIAS_PARAMS}/models/bias.h5" ]]; then
  echo "[$(timestamp)] Bias model already present, skipping submit."
else
  activate_chrombpnet
  BIAS_CLUSTER="${BIAS_CLUSTER}" \
  BIAS_MIN_THRESH="${BIAS_MIN_THRESH}" \
  BIAS_FOLD="${BIAS_FOLD}" \
  DATA_TYPE="${DATA_TYPE}" \
  sbatch ${BIAS_SBATCH_ARGS} 01-train_bias_model.sh
fi

wait_for_count \
  "Bias model" \
  "1" \
  "if [[ -f '${bias_dir%/}/${BIAS_PARAMS}/models/bias.h5' ]]; then echo 1; else echo 0; fi"

echo "[$(timestamp)] Step 02: ChromBPNet model training submit"
activate_chrombpnet
BIAS_PARAMS="${BIAS_PARAMS}" DATA_TYPE="${DATA_TYPE}" bash 02-train_chrombpnet_models.sh

wait_for_count \
  "Completed model evaluations" \
  "${EXPECTED_FOLDS}" \
  "find '${models_dir%/}/bias_${BIAS_PARAMS}' -path '*/evaluation/overall_report.html' -type f | wc -l"

echo "[$(timestamp)] Step 03: model QC"
QC_HTML="${SCRIPT_DIR}/03-model_QC.html"
if [[ "${FORCE_QC}" != "1" && -f "${qc_dir%/}/chrombpnet_metrics.tsv" && -f "${chrombpnet_models_keep}" && -f "${QC_HTML}" ]]; then
  echo "[$(timestamp)] QC outputs already present, skipping rerender of 03-model_QC.Rmd."
else
  activate_r
  BIAS_PARAMS="${BIAS_PARAMS}" CHROMBPNET_WORK_DIR="${base_dir}" Rscript -e "rmarkdown::render('03-model_QC.Rmd')"
fi

echo "[$(timestamp)] Step 04: contribution scores submit"
activate_chrombpnet
BIAS_PARAMS="${BIAS_PARAMS}" bash 04-get_contrib_scores.sh

wait_for_count \
  "Completed fold contribution exports" \
  "${EXPECTED_FOLDS}" \
  "find '${contribs_dir%/}/bias_${BIAS_PARAMS}' -path '*/peaks_shap.counts_scores.bw' -type f | wc -l"

wait_for_count \
  "Completed fold profile exports" \
  "${EXPECTED_FOLDS}" \
  "find '${contribs_dir%/}/bias_${BIAS_PARAMS}' -path '*/peaks_shap.profile_scores.bw' -type f | wc -l"

echo "[$(timestamp)] Step 05: fold-average DeepSHAP submit"
activate_chrombpnet
BIAS_PARAMS="${BIAS_PARAMS}" bash 05-average_deepshaps.sh

QC_DATASETS=$(dataset_count_keep)
wait_for_count \
  "Average counts shaps" \
  "${QC_DATASETS}" \
  "find '${contribs_dir%/}/bias_${BIAS_PARAMS}' -mindepth 2 -maxdepth 2 -name 'average_shaps.counts.h5' -type f | wc -l"

wait_for_count \
  "Average profile shaps" \
  "${QC_DATASETS}" \
  "find '${contribs_dir%/}/bias_${BIAS_PARAMS}' -mindepth 2 -maxdepth 2 -name 'average_shaps.profile.h5' -type f | wc -l"

echo "[$(timestamp)] Step 06: MoDISco submit"
activate_modiscolite
BIAS_PARAMS="${BIAS_PARAMS}" MEME_DB="${MEME_DB}" bash 06-modisco.sh

wait_for_count \
  "MoDISco outputs" \
  "${QC_DATASETS}" \
  "find '${modisco_dir%/}/bias_${BIAS_PARAMS}' -mindepth 2 -maxdepth 2 -name 'counts_modisco_output.h5' -type f | wc -l"

wait_for_count \
  "MoDISco reports" \
  "${QC_DATASETS}" \
  "find '${modisco_dir%/}/bias_${BIAS_PARAMS}' -path '*/counts_modisco_report/motifs.html' -type f | wc -l"

echo "[$(timestamp)] Step 06b/06c: post-MoDISco R summaries"
activate_r
BIAS_PARAMS="${BIAS_PARAMS}" Rscript 06b-modisco_report_to_tsv.R
BIAS_PARAMS="${BIAS_PARAMS}" Rscript 06c-post_modisco_qc.R

QC2_DATASETS=$(dataset_count_keep2)

echo "[$(timestamp)] Step 07: contribution bigwigs submit"
activate_chrombpnet
BIAS_PARAMS="${BIAS_PARAMS}" bash 07-make_contrib_bw.sh

wait_for_count \
  "Average counts contribution bigwigs" \
  "${QC2_DATASETS}" \
  "find '${contribs_dir%/}/bias_${BIAS_PARAMS}' -mindepth 2 -maxdepth 2 -name 'average_shaps.counts.bw' -type f | wc -l"

wait_for_count \
  "Average profile contribution bigwigs" \
  "${QC2_DATASETS}" \
  "find '${contribs_dir%/}/bias_${BIAS_PARAMS}' -mindepth 2 -maxdepth 2 -name 'average_shaps.profile.bw' -type f | wc -l"

echo "[$(timestamp)] Step 08: predictions submit"
activate_chrombpnet
BIAS_PARAMS="${BIAS_PARAMS}" bash 08-generate_predictions.sh

wait_for_count \
  "Bias-corrected prediction bigwigs" \
  "${QC2_DATASETS}" \
  "find '${preds_dir%/}/bias_corrected' -maxdepth 1 -name '*_chrombpnet_nobias.bw' -type f | wc -l"

wait_for_count \
  "Uncorrected prediction bigwigs" \
  "${QC2_DATASETS}" \
  "find '${preds_dir%/}/uncorrected' -maxdepth 1 -name '*_chrombpnet_uncorrected.bw' -type f | wc -l"

echo "[$(timestamp)] Step 09: final QC table"
activate_r
Rscript 09-prep_chrombpnet_QC_table.R

echo "[$(timestamp)] All stages completed."
