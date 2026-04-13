#!/bin/bash
set -euo pipefail

# Purpose: this script does a deep MoDisco run with 1M seqlets for each of the
# positive and negative sets, on averaged contribution scores for counts head,
# for each cell type.

# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set bias model
bias_params="${BIAS_PARAMS:-1-col_aspn_ogna_thresh0.4}"

# set MODISCO params
max_seqlets=1000000
num_leiden=2
meme_db="${MEME_DB:-}"
num_matches=10

# set outdir on scratch for faster IO
in_dir="${contribs_scratch%/}/bias_${bias_params}"
out_dir="${modisco_scratch%/}/bias_${bias_params}"
[[ -d "${out_dir}" ]] || mkdir -p "${out_dir}"

JOBSCRIPT=06-jobscript.sh

slurm_job_active() {
  local job_name="$1"
  squeue -h -u "${USER}" -n "${job_name}" -t PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED 2>/dev/null | grep -q .
}

read -r -a modisco_sbatch_extra_args <<< "${CHROMBPNET_MODISCO_SBATCH_ARGS:-}"
dataset_filter_regex="${CHROMBPNET_DATASET_FILTER_REGEX:-}"

if [[ -z "${meme_db}" ]]; then
  echo "MEME_DB is not set. Export MEME_DB=/path/to/motifs.meme before running step 06."
  exit 1
fi

if [[ ! -f "${meme_db}" ]]; then
  echo "MEME_DB does not exist: ${meme_db}"
  exit 1
fi

# CONSTRUCT COMMANDS -----------------------------------------------------------

if [[ -f "${chrombpnet_models_keep}" ]]; then
  datasets=$(awk '{print $1}' "${chrombpnet_models_keep}")
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

  # set up common inputs/outputs
	img_suffix_dir="./"
	
	# run for COUNTS -----------------------------------------
	job_name="06-modisco_${dataset}_counts"
	counts_output_memedb="${out_dir}/${dataset}/${dataset}_memedb.counts.txt"
	counts_peak_shaps="${in_dir}/${dataset}/average_shaps.counts.h5"
	counts_modisco_output="${out_dir}/${dataset}/counts_modisco_output.h5"
	counts_report_outdir="${out_dir}/${dataset}/counts_modisco_report"
	[[ -d "${counts_report_outdir}" ]] || mkdir -p "${counts_report_outdir}"

  # check if outs exist already
  if [[ -f $counts_modisco_output && -f ${counts_report_outdir}/motifs.html ]]; then 
    echo "@ done ${dataset} modisco on counts."
  elif slurm_job_active "${job_name}"; then
    echo "@ ${dataset} modisco already queued or running; skipping resubmission."
  else
    echo "@ Running modisco with: "
  	echo -e "\t job name: ${job_name}"
  	echo -e "\t counts shaps: ${counts_peak_shaps}"
  	echo -e "\t counts modisco output: ${counts_modisco_output}"
  	echo -e "\t counts report out: ${counts_report_outdir}"
  
  	sbatch "${modisco_sbatch_extra_args[@]}" -J "${job_name}" "./${JOBSCRIPT}" "${counts_peak_shaps}" \
                                       "${max_seqlets}" \
                                       "${num_leiden}" \
                                       "${counts_modisco_output}" \
                                       "${counts_report_outdir}" \
                                       "${img_suffix_dir}" \
                                       "${meme_db}" \
                                       "${num_matches}" \
                                       "${counts_output_memedb}"

  	sleep 2s
  fi
	
done
