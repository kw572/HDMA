#!/bin/bash
set -euo pipefail

# Purpose: train the bias-factorized ChromBPNet models.
# This script loops through all clusters, and for each cluster and each
# fold, constructs an sbatch command using the jobscript in 02-jobscript.sh,
# to train a model on that cluster/chromosome fold, via the chrombpnet pipeline command.
# https://github.com/kundajelab/chrombpnet/wiki/ChromBPNet-training
# The bias model used is set in the bias_params parameter.

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set bias model
bias_params="${chrombpnet_bias_params}"
bias_model="${bias_dir%/}/${bias_params}/models/bias.h5"

# set parameters
out_dir="${models_dir%/}/bias_${bias_params}/"
num_folds="${chrombpnet_train_num_folds}"

# make outdir
echo ${out_dir}
mkdir -p "${out_dir}"

slurm_job_active() {
	local job_name="$1"
	squeue -h -u "${USER}" -n "${job_name}" -t PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED 2>/dev/null | grep -q .
}

read -r -a train_sbatch_extra_args <<< "${CHROMBPNET_TRAIN_SBATCH_ARGS:---partition=gpu}"
dataset_filter_regex="${CHROMBPNET_DATASET_FILTER_REGEX:-}"



# CONSTRUCT COMMANDS -----------------------------------------------------------

all_pseudorep_dir="${cluster_frags_dir%/}/fragments"

# for every completed peaks file, extract [cluster]
datasets=$(
  find "${chrombpnet_peaks_dir}" -maxdepth 1 -type f -name '*__peaks_bpnet.narrowPeak' ! -name '._*' -print |
    while IFS= read -r path; do
      basename "${path}" __peaks_bpnet.narrowPeak
    done
)

if [[ -n "${dataset_filter_regex}" ]]; then
	datasets=$(printf '%s\n' ${datasets} | grep -E "${dataset_filter_regex}" || true)
fi

if [[ -z "${datasets}" ]]; then
	echo "No datasets matched CHROMBPNET_DATASET_FILTER_REGEX='${dataset_filter_regex}'"
	exit 1
fi

# datasets=( Eye_c21 )

for dataset in ${datasets}; do
	echo ${dataset}

	dataset_dir="${out_dir%/}/${dataset}"
	[[ -d "${dataset_dir}" ]] || mkdir -p "${dataset_dir}"

	for ((fold = 0; fold < num_folds; fold++)); do
		fold_name=fold_${fold}
		job_name="02-train_${dataset}_${fold_name}"
		echo -e "\t${fold_name}"

		# this is the output directory for that fold, for that cluster
		fold_dir="${dataset_dir}/${fold_name}/"
		[[ -d "${fold_dir}" ]] || mkdir -p "${fold_dir}"

		if [[ -f "${fold_dir}/evaluation/overall_report.html" ]]; then
			echo -e "\t\ttraining completed, skipping..."
		elif slurm_job_active "${job_name}"; then
			echo -e "\t\tjob already queued or running, skipping resubmission..."
		else
			JOBSCRIPT="./02-jobscript.sh"

	  		sleep 5s

			# If a prior run left a partial fold directory behind, wipe it before
			# resubmitting so chrombpnet can recreate its output tree cleanly.
			if [[ -d "${fold_dir}" ]] && [[ -n "$(find "${fold_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
				echo "existing partial outputs found under ${fold_dir}; removing fold directory before resubmission"
				rm -rf "${fold_dir}"
				mkdir -p "${fold_dir}"
			fi

			echo "Running chrombpnet pipeline for ${dataset} ${fold_name}"
			echo "config=${PWD}/../config.sh"
			echo "bias_model=${bias_model}"
			echo "output_dir=${fold_dir}"

			sbatch "${train_sbatch_extra_args[@]}" \
				--export=ALL,CHROMBPNET_TRAIN_DATASET="${dataset}",CHROMBPNET_TRAIN_FOLD="${fold_name}",CHROMBPNET_TRAIN_BIAS_PARAMS="${bias_params}",CHROMBPNET_TRAIN_OUT_DIR="${fold_dir}" \
				-J "${job_name}" \
				"${JOBSCRIPT}"
		fi
	done
	sleep 20s
done
