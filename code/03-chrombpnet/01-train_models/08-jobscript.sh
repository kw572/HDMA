#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/08/%x-%j.out
#SBATCH --partition=gpu
#SBATCH -t 2-0
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH --gres=gpu:1
#SBATCH --requeue
#SBATCH --open-mode=append

set -euo pipefail
set +u
source /etc/profile

module purge
module load legacy/CentOS7
module load gcc/8.3.0
module load cuda/11.2.0
module load cudnn/8.1.0.77-11.2-cuda

# optional graphics deps
module load cairo || true
module load pango || true

eval "$(conda shell.bash hook)"
conda activate chrombpnet
set -u

echo "=== NVIDIA ==="
nvidia-smi
celltype="${1}"
peaks_file="${2}"
ref_fasta="${3}"
chrom_sizes="${4}"
out_prefix="${5}"
out_key="${6}"
shift 6
model_paths=("$@")
predict_batch_size="${PREDICT_BATCH_SIZE:-16}"

echo "celltype=${celltype}"
echo "peaks_file=${peaks_file}"
echo "ref_fasta=${ref_fasta}"
echo "chrom_sizes=${chrom_sizes}"
echo "out_prefix=${out_prefix}"
echo "out_key=${out_key}"
printf 'model_paths=%s\n' "${model_paths[*]}"

if [[ "${#model_paths[@]}" -eq 0 ]]; then
  echo "No model paths were provided."
  exit 1
fi


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

load_requested_modules() {
  local modules_string="$1"
  local mod
  [[ -n "${modules_string}" ]] || return 0
  read -r -a modules <<< "${modules_string}"
  for mod in "${modules[@]}"; do
    load_optional_module "${mod}"
  done
}

PRE_MODULES="${PRE_MODULES:-legacy/CentOS7 gcc/8.3.0}"
CUDA_MODULE="${CUDA_MODULE:-cuda/11.2.0}"
CUDNN_MODULE="${CUDNN_MODULE:-cudnn/8.1.0.77-11.2}"
load_requested_modules "${PRE_MODULES}"
load_optional_module "${CUDA_MODULE}"
load_optional_module "${CUDNN_MODULE}"
for mod in system cairo pango; do
  load_optional_module "${mod}"
done

echo "[$(date +"%m/%d/%Y (%r)")] starting ${celltype}"

predict_cmd=(
  python ./08-predict_and_avg.py
  --regions "${peaks_file}"
  --genome "${ref_fasta}"
  --chrom-sizes "${chrom_sizes}"
  --output-prefix "${out_prefix}"
  --output-key "${out_key}"
  --output-bed True
  --batch-size "${predict_batch_size}"
)

for model_path in "${model_paths[@]}"; do
  predict_cmd+=(--chrombpnet-model "${model_path}")
done

"${predict_cmd[@]}"

echo "[$(date +"%m/%d/%Y (%r)")] done!"
