#!/bin/bash
set -euo pipefail

# Purpose: this script generates bigwigs based on the average contribution scores
# across folds, for each cell type.

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set bias model
bias_params="${BIAS_PARAMS:-${chrombpnet_bias_params}}"

# set outdir on scratch for faster IO
in_dir="${contribs_scratch%/}/bias_${bias_params}"

JOBSCRIPT=07-jobscript.sh

slurm_job_active() {
  local job_name="$1"
  squeue -h -u "${USER}" -n "${job_name}" -t PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED 2>/dev/null | grep -q .
}

read -r -a contrib_bw_sbatch_extra_args <<< "${CHROMBPNET_CONTRIB_BW_SBATCH_ARGS:---partition=main}"
dataset_filter_regex="${CHROMBPNET_DATASET_FILTER_REGEX:-}"

# CONSTRUCT COMMANDS -----------------------------------------------------------

if [[ -f "${chrombpnet_models_keep2}" ]]; then
  datasets=$(awk '{print $1}' "${chrombpnet_models_keep2}")
else
  datasets=$(
    find "${in_dir}" -mindepth 1 -maxdepth 1 -type d -print |
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

for dataset in ${datasets}; do

	echo ${dataset}
	job_name="07-contrib_bw_${dataset}"
	
	peaks_file=${in_dir%/}/${dataset}/fold_0/peaks_shap.interpreted_regions.bed
 
	counts_peak_shaps=${in_dir%/}/${dataset}/average_shaps.counts.h5
	counts_out_prefix=${in_dir%/}/${dataset}/average_shaps.counts
	
	profile_peak_shaps=${in_dir%/}/${dataset}/average_shaps.profile.h5
	profile_out_prefix=${in_dir%/}/${dataset}/average_shaps.profile
	
	if [[ -f "${counts_out_prefix}.bw" && -f "${profile_out_prefix}.bw" ]]; then
	  echo -e "@ Found bigwigs for ${dataset}; skipping."
  elif slurm_job_active "${job_name}"; then
    echo -e "@ Bigwig job for ${dataset} is already queued or running; skipping resubmission."
  else 
  	echo "@ Generating bigwigs for: "
  	echo -e "\t counts: ${counts_peak_shaps} to ${counts_out_prefix}.bw"
  	echo -e "\t profile: ${profile_peak_shaps} to ${profile_out_prefix}.bw"
  	echo -e "\t using ${peaks_file}"
  	sbatch "${contrib_bw_sbatch_extra_args[@]}" -J "${job_name}" "./${JOBSCRIPT}" "${dataset}" \
                                       "${counts_peak_shaps}" \
                                       "${counts_out_prefix}" \
                                       "${profile_peak_shaps}" \
                                       "${profile_out_prefix}" \
                                       "${peaks_file}" \
                                       "${chromsizes}"
  fi
  
  sleep 5s

done
