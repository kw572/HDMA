#!/bin/bash
set -euo pipefail

# Purpose: this script generates bigwigs of the predicted bias-corrected or uncorrected
# chromatin accessibility profiles.

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set bias model
bias_params="${BIAS_PARAMS:-${chrombpnet_bias_params}}"
model_dir="${models_dir%/}/bias_${bias_params}"
ref_fasta="${ref_fasta}"

echo $model_dir
echo $ref_fasta
echo $chromsizes

JOBSCRIPT=08-jobscript.sh

slurm_job_active() {
  local job_name="$1"
  squeue -h -u "${USER}" -n "${job_name}" -t PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED 2>/dev/null | grep -q .
}

read -r -a predict_sbatch_extra_args <<< "${CHROMBPNET_PREDICT_SBATCH_ARGS:---partition=gpu}"
dataset_filter_regex="${CHROMBPNET_DATASET_FILTER_REGEX:-}"
num_folds="${chrombpnet_train_num_folds}"
effective_folds=$(( num_folds < 5 ? num_folds : 5 ))

completed_folds_for_dataset() {
  local dataset="$1"
  local folds_keep=""
  if [[ -f "${chrombpnet_models_keep2}" ]]; then
    folds_keep=$(awk -v dataset="${dataset}" '$1 == dataset {print $2; exit}' "${chrombpnet_models_keep2}")
  fi

  if [[ -z "${folds_keep}" && -f "${chrombpnet_models_keep}" ]]; then
    folds_keep=$(awk -v dataset="${dataset}" '$1 == dataset {print $2; exit}' "${chrombpnet_models_keep}")
  fi

  if [[ -n "${folds_keep}" ]]; then
    printf '%s\n' "${folds_keep}"
    return 0
  fi

  local fold
  local found=()
  for ((fold = 0; fold < effective_folds; fold++)); do
    if [[ -f "${model_dir}/${dataset}/fold_${fold}/evaluation/overall_report.html" ]]; then
      found+=("fold_${fold}")
    fi
  done

  if [[ "${#found[@]}" -gt 0 ]]; then
    IFS=,
    printf '%s\n' "${found[*]}"
    unset IFS
  fi
}

if [[ -f "${chrombpnet_models_keep2}" ]]; then
  datasets=$(awk '{print $1}' "${chrombpnet_models_keep2}")
else
  datasets=$(
    find "${model_dir}" -mindepth 1 -maxdepth 1 -type d -print |
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

# types of predictions
modes=("bias_corrected" "uncorrected")

for dataset in ${datasets}; do
  echo "Processing dataset: ${dataset}"

  for mode in "${modes[@]}"; do
    echo "  Mode: ${mode}"
    job_name="08-predict_${dataset}_${mode}"

    # construct paths based on the current mode
    if [[ "$mode" == "bias_corrected" ]]; then
      out_dir="${preds_scratch%/}/bias_${bias_params}/bias_corrected"
      final_out_dir="${preds_dir%/}/bias_${bias_params}/bias_corrected"
      out_key="nobias"
    else
      out_dir="${preds_scratch%/}/bias_${bias_params}/uncorrected"
      final_out_dir="${preds_dir%/}/bias_${bias_params}/uncorrected"
      out_key="uncorrected"
    fi
    mkdir -p "${out_dir}" "${final_out_dir}"

    folds_csv="$(completed_folds_for_dataset "${dataset}")"
    if [[ -z "${folds_csv}" ]]; then
      echo -e "\t\tno completed folds found for ${dataset}, skipping ${mode} predictions..."
      continue
    fi

    IFS=',' read -r -a folds <<< "${folds_csv}"
    model_paths=()
    for fold_name in "${folds[@]}"; do
      if [[ "$mode" == "bias_corrected" ]]; then
        model_file="${model_dir}/${dataset}/${fold_name}/models/chrombpnet_nobias.h5"
      else
        model_file="${model_dir}/${dataset}/${fold_name}/models/chrombpnet.h5"
      fi
      [[ -f "${model_file}" ]] || continue
      model_paths+=("${model_file}")
    done

    if [[ "${#model_paths[@]}" -eq 0 ]]; then
      echo -e "\t\tno model files found for ${dataset} (${mode}), skipping..."
      continue
    fi
  
    in_dir="${contribs_scratch%/}/bias_${bias_params}"
    peaks_file=${in_dir%/}/${dataset}/fold_0/peaks_shap.interpreted_regions.bed
    out_prefix="${out_dir%/}/${dataset}_avg"
    final_out_prefix="${final_out_dir%/}/${dataset}_avg"
    final_out_file1="${final_out_prefix}_chrombpnet_${out_key}_preds_w_logcounts.bed"
    final_out_file2="${final_out_prefix}_chrombpnet_${out_key}.bw"

    if [[ -f "${final_out_file1}" && -f "${final_out_file2}" ]]; then
      echo -e "\t\tfound completed ${mode} predictions for ${dataset}, skipping..."
    elif slurm_job_active "${job_name}"; then
      echo -e "\t\t${mode} prediction job already queued or running for ${dataset}, skipping resubmission..."
    else
      echo "Generating predictions for ${dataset} (${mode})"
      echo "Final Output File 1: ${final_out_file1}"
      echo "Final Output File 2: ${final_out_file2}"

      echo "@ generating ${out_prefix} ${out_key} predictions"

      echo "sbatch -J ${job_name} ${JOBSCRIPT} ${dataset}"
      echo " ${peaks_file}"
      echo " ${ref_fasta}"
      echo " ${chromsizes}"
      echo " ${out_prefix}"
      echo " ${out_key}"
      for model_file in "${model_paths[@]}"; do
        echo " ${model_file}"
      done

      sbatch "${predict_sbatch_extra_args[@]}" -J "${job_name}" "./${JOBSCRIPT}" "${dataset}" \
          "${peaks_file}" \
          "${ref_fasta}" \
          "${chromsizes}" \
          "${out_prefix}" \
          "${out_key}" \
          "${model_paths[@]}"

      sleep 3s
    fi
  done
done
