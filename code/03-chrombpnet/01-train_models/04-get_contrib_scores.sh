#!/bin/bash
set -euo pipefail

# Purpose: this script calls the chrombpnet helper to get contribution 
# score bigwigs, as computed with DeepLIFT: https://github.com/kundajelab/chrombpnet/wiki/Generate-contribution-score-bigwigs
# Based off code from Ryan Zhao / Salil Deshpande.

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set bias model
bias_params="${BIAS_PARAMS:-1-col_aspn_ogna_thresh0.4}"

# set parameters
ref_fasta="${ref_fasta}"
out_dir="${contribs_scratch%/}/bias_${bias_params}/"
[[ -d "${out_dir}" ]] || mkdir -p "${out_dir}"

slurm_job_active() {
  local job_name="$1"
  squeue -h -u "${USER}" -n "${job_name}" -t PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED 2>/dev/null | grep -q .
}

read -r -a contrib_sbatch_extra_args <<< "${CHROMBPNET_CONTRIB_SBATCH_ARGS:-}"
dataset_filter_regex="${CHROMBPNET_DATASET_FILTER_REGEX:-}"



# CONSTRUCT COMMANDS -----------------------------------------------------------

if [[ -f "${chrombpnet_models_keep}" ]]; then
  datasets=$(awk '{print $1}' "${chrombpnet_models_keep}")
else
  datasets=$(
    find "${models_dir%/}/bias_${bias_params}" -mindepth 1 -maxdepth 1 -type d -print |
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

	peaks_file="${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak"
	dataset_dir="${out_dir%/}/${dataset}"
	[[ -d "${dataset_dir}" ]] || mkdir -p "${dataset_dir}"

	for fold in {0..4}; do
		fold_name=fold_${fold}
		job_name="04-get_contrib_${dataset}_${fold_name}"
		echo -e "\t${fold_name}"

		model_file="${models_dir%/}/bias_${bias_params}/${dataset}/${fold_name}/models/chrombpnet_nobias.h5"
		
		[[ -f "${model_file}" ]] || { echo -e "\t\tModel file missing"; continue; }
		
		fold_dir="${dataset_dir}/${fold_name}"

		[[ -d "${fold_dir}" ]] || mkdir -p "${fold_dir}"

		out_prefix="${fold_dir}/peaks_shap"

		if [[ -f "${out_prefix}.counts_scores.bw" && -f "${out_prefix}.profile_scores.bw" ]]; then
			echo -e "\t\tfound completed shap scores, skipping..."
		elif slurm_job_active "${job_name}"; then
			echo -e "\t\tjob already queued or running, skipping resubmission..."
		else 
			JOBSCRIPT=04-jobscript.sh
			
			echo "Running chrombpnet contribs_bw"
			echo "--genome ${ref_fasta}"
			echo "--regions ${peaks_file}"
			echo "--model-h5 ${model_file}"
			echo "--output-prefix ${out_prefix}"
			echo "--chrom-sizes ${chromsizes}"
			
			sbatch "${contrib_sbatch_extra_args[@]}" -J "${job_name}" "./${JOBSCRIPT}" "${ref_fasta}" \
				"${peaks_file}" \
				"${model_file}" \
				"${out_prefix}" \
				"${chromsizes}"
			
			sleep 5s
		fi
	done
	sleep 20s
done
