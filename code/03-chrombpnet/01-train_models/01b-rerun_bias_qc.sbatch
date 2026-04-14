#!/bin/bash
#SBATCH --job-name=bias-qc
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=2-00:00:00
#SBATCH --output=../../logs/03-chrombpnet/01/%x-%j.out
#SBATCH --error=../../logs/03-chrombpnet/01/%x-%j.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
  set +u
  if [[ -f /etc/profile ]]; then
    source /etc/profile
  fi
  if [[ -f /usr/share/lmod/lmod/init/bash ]]; then
    source /usr/share/lmod/lmod/init/bash
  elif [[ -f /usr/share/Modules/init/bash ]]; then
    source /usr/share/Modules/init/bash
  fi
  set -u
fi

source ~/.bashrc 2>/dev/null || true

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
source ../config.sh

load_optional_module() {
  local mod="$1"
  [[ -n "${mod}" ]] || return 0
  if module -t avail "${mod}" 2>&1 | grep -Fq "${mod}"; then
    module load "${mod}"
  else
    echo "WARNING: module '${mod}' is unavailable; continuing without it."
  fi
}

eval "$(conda shell.bash hook)"
conda activate chrombpnet

CUDA_MODULE="${CUDA_MODULE:-cuda/11.2.0}"
CUDNN_MODULE="${CUDNN_MODULE:-}"
load_optional_module "${CUDA_MODULE}"
load_optional_module "${CUDNN_MODULE}"
for mod in system libxml2 libxslt perl zlib ghostscript cairo pango; do
  load_optional_module "${mod}"
done

BIAS_CLUSTER="${BIAS_CLUSTER:-1-col_aspn_ogna}"
BIAS_MIN_THRESH="${BIAS_MIN_THRESH:-0.4}"
BIAS_FOLD="${BIAS_FOLD:-fold_0}"
DATA_TYPE="${DATA_TYPE:-ATAC}"
BIAS_PARAMS="${BIAS_PARAMS:-${BIAS_CLUSTER}_thresh${BIAS_MIN_THRESH}}"
BIAS_OUT_DIR="${bias_dir%/}/${BIAS_PARAMS}"

BIAS_MODEL="${BIAS_OUT_DIR}/models/bias.h5"
BIAS_BIGWIG="${BIAS_OUT_DIR}/auxiliary/data_unstranded.bw"
PEAKS_FILE="${chrombpnet_peaks_dir}/${BIAS_CLUSTER}__peaks_bpnet.narrowPeak"
NEGATIVES_FILE="${negatives_dir}/${BIAS_CLUSTER}/${BIAS_FOLD}/output_negatives.bed"
SPLIT_FILE="${split_dir}/${BIAS_FOLD}.json"

for req in "${BIAS_MODEL}" "${BIAS_BIGWIG}" "${PEAKS_FILE}" "${NEGATIVES_FILE}" "${SPLIT_FILE}"; do
  if [[ ! -e "${req}" ]]; then
    echo "ERROR: required file missing: ${req}"
    exit 1
  fi
done

# Bias QC reruns can fail on stale partial outputs because chrombpnet creates
# these directories with exist_ok=False.
rm -rf \
  "${BIAS_OUT_DIR}/auxiliary/interpret_subsample" \
  "${BIAS_OUT_DIR}/evaluation/modisco_profile" \
  "${BIAS_OUT_DIR}/evaluation/modisco_counts"

echo "Running chrombpnet bias qc"
echo "--bigwig ${BIAS_BIGWIG}"
echo "--bias-model ${BIAS_MODEL}"
echo "--genome ${ref_fasta}"
echo "--chrom-sizes ${chromsizes}"
echo "--output-dir ${BIAS_OUT_DIR}"
echo "--data-type ${DATA_TYPE}"
echo "--peaks ${PEAKS_FILE}"
echo "--nonpeaks ${NEGATIVES_FILE}"
echo "--chr-fold-path ${SPLIT_FILE}"

chrombpnet bias qc \
  --bigwig "${BIAS_BIGWIG}" \
  --bias-model "${BIAS_MODEL}" \
  --genome "${ref_fasta}" \
  --chrom-sizes "${chromsizes}" \
  --output-dir "${BIAS_OUT_DIR}" \
  --data-type "${DATA_TYPE}" \
  --peaks "${PEAKS_FILE}" \
  --nonpeaks "${NEGATIVES_FILE}" \
  --chr-fold-path "${SPLIT_FILE}"
