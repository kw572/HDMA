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
bias_params="${BIAS_PARAMS:-1-col_aspn_ogna_thresh0.4}"
bias_model="${bias_dir%/}/${bias_params}/models/bias.h5"

# set parameters
data_type="${DATA_TYPE:-ATAC}"
ref_fasta="${ref_fasta}"
out_dir="${models_dir%/}/bias_${bias_params}/"

# make outdir
echo ${out_dir}
mkdir -p "${out_dir}"

slurm_job_active() {
	local job_name="$1"
	squeue -h -u "${USER}" -n "${job_name}" -t PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED 2>/dev/null | grep -q .
}

read -r -a train_sbatch_extra_args <<< "${CHROMBPNET_TRAIN_SBATCH_ARGS:-}"
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

	frag_file="${cluster_frags_dir%/}/fragments/${dataset}__sorted.tsv"
	peaks_file="${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak"

	dataset_dir="${out_dir%/}/${dataset}"
	[[ -d "${dataset_dir}" ]] || mkdir -p "${dataset_dir}"

	for fold in {0..4}; do
		fold_name=fold_${fold}
		job_name="02-train_${dataset}_${fold_name}"
		echo -e "\t${fold_name}"

		negatives_file="${negatives_dir%/}/${dataset}/${fold_name}/output_negatives.bed"
		split_file="${split_dir%/}/${fold_name}.json"

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
	  
	  		# clean old logs if it exists - echo the command instead of running it.
	  		if [[ -d "${fold_dir}/logs" ]]; then
	  			echo "existing logs found under ${fold_dir}; jobscript will clean interpret_subsample if needed"
			fi

			echo "Running chrombpnet pipeline"
			echo "--input-fragment-file ${frag_file}"
			echo "--genome ${ref_fasta}"
			echo "--chrom-sizes ${chromsizes}"
			echo "--peaks ${peaks_file}"
			echo "--nonpeaks ${negatives_file}"
			echo "--chr-fold-path ${split_file}"
			echo "--bias-model-path ${bias_model}"
			echo "--output-dir ${fold_dir}"
			echo "--data-type ${data_type}"

			sbatch "${train_sbatch_extra_args[@]}" -J "${job_name}" "${JOBSCRIPT}" "${frag_file}" \
						"${ref_fasta}" \
						"${chromsizes}" \
						"${peaks_file}" \
						"${negatives_file}" \
						"${split_file}" \
						"${bias_model}" \
						"${fold_dir}" \
						"${data_type}"
		fi
	done
	sleep 20s
done
