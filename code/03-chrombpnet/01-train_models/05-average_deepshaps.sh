#!/bin/bash
set -euo pipefail

# Purpose: average the Shap scores across folds, for counts and profile
# contribution scores separately.

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set bias model
bias_params="${BIAS_PARAMS:-${chrombpnet_bias_params}}"
out_dir="${contribs_scratch%/}/bias_${bias_params}/"

echo "@ contribs dir: "
echo ${out_dir}

slurm_job_active() {
  local job_name="$1"
  squeue -h -u "${USER}" -n "${job_name}" -t PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED 2>/dev/null | grep -q .
}

read -r -a avg_sbatch_extra_args <<< "${CHROMBPNET_AVG_SBATCH_ARGS:---partition=gpu}"
dataset_filter_regex="${CHROMBPNET_DATASET_FILTER_REGEX:-}"
num_folds="${chrombpnet_train_num_folds}"
effective_folds=$(( num_folds < 5 ? num_folds : 5 ))

completed_folds_for_dataset() {
  local dataset="$1"
  local folds_keep=""
  if [[ -f "${chrombpnet_models_keep}" ]]; then
    folds_keep=$(awk -v dataset="${dataset}" '$1 == dataset {print $2; exit}' "${chrombpnet_models_keep}")
  fi

  if [[ -n "${folds_keep}" ]]; then
    printf '%s\n' "${folds_keep}"
    return 0
  fi

  local fold
  local found=()
  for ((fold = 0; fold < effective_folds; fold++)); do
    if [[ -f "${models_dir%/}/bias_${bias_params}/${dataset}/fold_${fold}/evaluation/overall_report.html" ]]; then
      found+=("fold_${fold}")
    fi
  done

  if [[ "${#found[@]}" -gt 0 ]]; then
    IFS=,
    printf '%s\n' "${found[*]}"
    unset IFS
  fi
}

# CONSTRUCT COMMANDS -----------------------------------------------------------

if [[ -f "${chrombpnet_models_keep}" ]]; then
  datasets=$(awk '{print $1}' "${chrombpnet_models_keep}")
else
  datasets=$(
    find "${out_dir}" -mindepth 1 -maxdepth 1 -type d -print |
      while IFS= read -r path; do
        basename "${path}"
      done
  )
fi

if [[ -n "${dataset_filter_regex}" ]]; then
  datasets=$(printf '%s\n' ${datasets} | grep -E "${dataset_filter_regex}" || true)
fi

if [[ -z "${datasets}" ]]; then
  echo "No datasets matched CHROMBPNET_DATASET_FILTER_REGEX='${dataset_filter_regex}'"
  exit 1
fi

# DEBUG
# datasets=( Muscle_c7 Thymus_c16 Thyroid_c10 )

# organ="Heart"
# datasets_organ=( $(echo ${datasets[@]} | tr ' ' '\n' | grep ${organ}) )

for dataset in ${datasets}; do
  
  echo ${dataset}
  dataset_dir="${out_dir%/}/${dataset}"
  echo ${dataset_dir}

  folds_keep="$(completed_folds_for_dataset "${dataset}")"
  if [[ -z "${folds_keep}" ]]; then
    echo "@ no completed folds found for ${dataset}; skipping averaging."
    continue
  fi

  echo ${folds_keep}
  
  job_name="05-avg_${dataset}"
  JOBSCRIPT=05-jobscript.sh
 
  if [[ -f "${dataset_dir}/average_shaps.counts.h5" && -f "${dataset_dir}/average_shaps.profile.h5" ]]; then
    echo "@ done ${dataset} averaging."
  elif slurm_job_active "${job_name}"; then
    echo "@ ${dataset} averaging already queued or running; skipping resubmission."
  else 
    echo "Submitting python 05-average_deepshaps.py ${dataset_dir} ${folds_keep} {counts,profile}"
    sbatch "${avg_sbatch_extra_args[@]}" -J "${job_name}" "./${JOBSCRIPT}" "${dataset_dir}" "${folds_keep}"
    sleep 5s
  fi

done
