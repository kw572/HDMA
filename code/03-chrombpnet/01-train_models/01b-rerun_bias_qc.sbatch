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

RUN_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "${RUN_DIR}"

resolve_script_dir() {
  local candidate
  for candidate in \
    "${RUN_DIR}" \
    "${RUN_DIR}/code/03-chrombpnet/01-train_models" \
    "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" 2>/dev/null && pwd)"
  do
    if [[ -n "${candidate}" && -f "${candidate}/../config.sh" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

SCRIPT_DIR="$(resolve_script_dir)" || {
  echo "ERROR: could not locate 01-train_models directory from RUN_DIR=${RUN_DIR}"
  exit 1
}
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
FRAG_FILE="${cluster_frags_dir}/fragments/${BIAS_CLUSTER}__sorted.tsv"
PEAKS_FILE="${chrombpnet_peaks_dir}/${BIAS_CLUSTER}__peaks_bpnet.narrowPeak"
NEGATIVES_FILE="${negatives_dir}/${BIAS_CLUSTER}/${BIAS_FOLD}/output_negatives.bed"
SPLIT_FILE="${split_dir}/${BIAS_FOLD}.json"

for req in "${BIAS_MODEL}" "${FRAG_FILE}" "${PEAKS_FILE}" "${NEGATIVES_FILE}" "${SPLIT_FILE}"; do
  if [[ ! -e "${req}" ]]; then
    echo "ERROR: required file missing: ${req}"
    exit 1
  fi
done

if [[ ! -f "${BIAS_BIGWIG}" ]]; then
  echo "Bias QC input bigwig is missing; regenerating ${BIAS_BIGWIG} from ${FRAG_FILE}"
  mkdir -p "${BIAS_OUT_DIR}/auxiliary"
  READS_TO_BIGWIG=$(python - <<'PY'
import chrombpnet.helpers.preprocessing.reads_to_bigwig as m
print(m.__file__)
PY
)
  python "${READS_TO_BIGWIG}" \
    --genome "${ref_fasta}" \
    --input-fragment-file "${FRAG_FILE}" \
    --chrom-sizes "${chromsizes}" \
    --output-prefix "${BIAS_OUT_DIR}/auxiliary/data" \
    --data-type "${DATA_TYPE}"
fi

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
