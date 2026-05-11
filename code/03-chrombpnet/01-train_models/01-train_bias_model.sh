#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/%x-%j.out
#SBATCH --partition=gpu
#SBATCH -t 2-0
#SBATCH -c 4
#SBATCH --mem=60G
#SBATCH --gres=gpu:1

# Purpose: trains an enzymatic bias model for a specific cluster, with a specific
# bias-threshold-factor, using the chrombpnet bias pipeline command, following
# https://github.com/kundajelab/chrombpnet/wiki/Bias-model-training.
# The specific cluster / threshold to use are set in the
# bias_cluster and min_thresh params. Multiple models are trained and the best
# one (which captures Tn5 sequence preferences but not the motifs of real TFs)
# is used for training all bias-factorized models.



# ENVIRONMENT ------------------------------------------------------------------

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

load_optional_module() {
  local mod="$1"
  [[ -n "${mod}" ]] || return 0
  if module -t avail "${mod}" 2>&1 | grep -Fq "${mod}"; then
    module load "${mod}"
  else
    echo "WARNING: module '${mod}' is unavailable; continuing without it."
  fi
}

load_requested_modules() {
  local modules_string="$1"
  local mod
  [[ -n "${modules_string}" ]] || return 0
  read -r -a modules <<< "${modules_string}"
  for mod in "${modules[@]}"; do
    load_optional_module "${mod}"
  done
}

# load modules
PRE_MODULES="${PRE_MODULES:-legacy/CentOS7 gcc/8.3.0}"
CUDA_MODULE="${CUDA_MODULE:-cuda/11.2.0}"
CUDNN_MODULE="${CUDNN_MODULE:-cudnn/8.1.0.77-11.2}"
load_requested_modules "${PRE_MODULES}"
load_optional_module "${CUDA_MODULE}"
load_optional_module "${CUDNN_MODULE}"
for mod in system libxml2 libxslt perl zlib ghostscript cairo pango; do
  load_optional_module "${mod}"
done



# PARAMETERS -------------------------------------------------------------------

# source configuration variables
source ../config.sh

# set parameters
min_thresh="${BIAS_MIN_THRESH:-0.4}"
bias_cluster="${BIAS_CLUSTER:-1_Jaw_Hyoid__peaks_bpnet}"
bias_fold="${BIAS_FOLD:-fold_0}"
data_type="${DATA_TYPE:-ATAC}"

ref_fasta="${ref_fasta}"
frag_file="${cluster_frags_dir}/fragments/${bias_cluster}__sorted.tsv"
peaks_file="${chrombpnet_peaks_dir}/${bias_cluster}__peaks_bpnet.narrowPeak"
negatives_file="${negatives_dir}/${bias_cluster}/${bias_fold}/output_negatives.bed"
split_file="${split_dir}/${bias_fold}.json"
out_dir="${bias_dir%/}/${bias_cluster}_thresh${min_thresh}/"

mkdir -p "${out_dir}"



# RUN --------------------------------------------------------------------------
                                   
echo "Running: chrombpnet bias pipeline --genome ${ref_fasta}"
echo "--input-fragment-file ${frag_file}"
echo "--peaks ${peaks_file}"
echo "--nonpeaks ${negatives_file}"
echo "--chr-fold-path ${split_file}"
echo "--bias-threshold-factor ${min_thresh}"
echo "--output-dir ${out_dir}"
echo "--chrom-sizes ${chromsizes}"
echo "--data-type ${data_type}"
                                 
chrombpnet bias pipeline --genome "${ref_fasta}" \
                         --input-fragment-file "${frag_file}" \
                         --peaks "${peaks_file}" \
                         --nonpeaks "${negatives_file}" \
                         --chr-fold-path "${split_file}" \
                         --bias-threshold-factor "${min_thresh}" \
                         --output-dir "${out_dir}" \
                         --chrom-sizes "${chromsizes}" \
                         --data-type "${data_type}"
                                   
