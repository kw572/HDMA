#!/bin/bash
set -euo pipefail

# Purpose: average the Shap scores across folds, for counts and profile
# contribution scores separately.

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set bias model
bias_params="${BIAS_PARAMS:-1-col_aspn_ogna_thresh0.4}"
out_dir="${contribs_scratch%/}/bias_${bias_params}/"

echo "@ contribs dir: "
echo ${out_dir}

slurm_job_active() {
  local job_name="$1"
  squeue -h -u "${USER}" -n "${job_name}" -t PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED 2>/dev/null | grep -q .
}

read -r -a avg_sbatch_extra_args <<< "${CHROMBPNET_AVG_SBATCH_ARGS:---partition=main}"
dataset_filter_regex="${CHROMBPNET_DATASET_FILTER_REGEX:-}"

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

  # grep the models_keep file for folds to keep
  if [[ -f "${chrombpnet_models_keep}" ]]; then
    folds_keep=$(grep -E "^${dataset}[[:space:]]" "${chrombpnet_models_keep}" | awk '{print $2}')
  else
    folds_keep="fold_0,fold_1,fold_2,fold_3,fold_4"
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
