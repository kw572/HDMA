#!/bin/bash
set -euo pipefail
#SBATCH --output=../../logs/03-chrombpnet/01/08/%x-%j.out
#SBATCH -p akundaje,owners,gpu,wjg
#SBATCH -t 2-0
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH --gres=gpu:a40:1
#SBATCH --requeue
#SBATCH --open-mode=append

celltype="${1}"
peaks_file="${2}"
ref_fasta="${3}"
chrom_sizes="${4}"
out_prefix="${5}"
out_key="${6}"
model_fold_0="${7}"
model_fold_1="${8}"
model_fold_2="${9}"
model_fold_3="${10}"
model_fold_4="${11}"
predict_batch_size="${PREDICT_BATCH_SIZE:-16}"

echo "celltype=${1}"
echo "peaks_file=${2}"
echo "ref_fasta=${3}"
echo "chrom_sizes=${4}"
echo "out_prefix=${5}"
echo "out_key=${6}"
echo "model_fold_0=${7}"
echo "model_fold_1=${8}"
echo "model_fold_2=${9}"
echo "model_fold_3=${10}"
echo "model_fold_4=${11}"


# source configuration variables
source ../config.sh

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

CUDA_MODULE="${CUDA_MODULE:-cuda/11.2}"
CUDNN_MODULE="${CUDNN_MODULE:-cudnn/8.1}"
load_optional_module "${CUDA_MODULE}"
load_optional_module "${CUDNN_MODULE}"
for mod in system cairo pango; do
  load_optional_module "${mod}"
done

echo "[$(date +"%m/%d/%Y (%r)")] starting ${celltype}"

python ./08-predict_and_avg.py --regions "${peaks_file}" \
  --genome "${ref_fasta}" \
  --chrom-sizes "${chrom_sizes}" \
  --output-prefix "${out_prefix}" \
  --output-key "${out_key}" \
  --output-bed True \
  --batch-size "${predict_batch_size}" \
  --chrombpnet-model "${model_fold_0}" \
  --chrombpnet-model "${model_fold_1}" \
  --chrombpnet-model "${model_fold_2}" \
  --chrombpnet-model "${model_fold_3}" \
  --chrombpnet-model "${model_fold_4}"

echo "[$(date +"%m/%d/%Y (%r)")] done!"
