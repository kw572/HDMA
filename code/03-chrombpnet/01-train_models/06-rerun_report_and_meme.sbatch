#!/bin/bash
#SBATCH --job-name=chrombpnet-modisco-repair
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=../../logs/03-chrombpnet/01/06/%x-%j.out
#SBATCH --error=../../logs/03-chrombpnet/01/06/%x-%j.err

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
cd "$RUN_DIR"

resolve_script_dir() {
  local candidate
  for candidate in \
    "$RUN_DIR" \
    "$RUN_DIR/code/03-chrombpnet/01-train_models" \
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
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] ERROR: could not locate 01-train_models directory from RUN_DIR=${RUN_DIR}"
  exit 1
}
cd "$SCRIPT_DIR"

source "$SCRIPT_DIR/../config.sh"

BIAS_PARAMS="${BIAS_PARAMS:-${chrombpnet_bias_params}}"
DATASET_FILTER_REGEX="${CHROMBPNET_DATASET_FILTER_REGEX:-}"
REGENERATE_REPORT="${REGENERATE_REPORT:-0}"
REGENERATE_MEME="${REGENERATE_MEME:-0}"

MODISCO_BASE="${modisco_dir%/}/bias_${BIAS_PARAMS}"

if [[ ! -d "${MODISCO_BASE}" ]]; then
  echo "MoDISco output directory not found: ${MODISCO_BASE}"
  exit 1
fi

if [[ ! -f "${MEME_DB}" ]]; then
  echo "MEME database not found: ${MEME_DB}"
  exit 1
fi

eval "$(conda shell.bash hook)"
conda activate modiscolite

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export PATH="$PATH:$HOME/.conda/envs/chrombpnet/bin:$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.5"

echo "MoDISco base: ${MODISCO_BASE}"
echo "MEME DB: ${MEME_DB}"
echo "Dataset filter: ${DATASET_FILTER_REGEX:-<all datasets>}"
echo "Regenerate report: ${REGENERATE_REPORT}"
echo "Regenerate meme: ${REGENERATE_MEME}"
echo "modisco: $(command -v modisco)"
echo "tomtom: $(command -v tomtom || true)"

for dataset_dir in "${MODISCO_BASE}"/*; do
  [[ -d "${dataset_dir}" ]] || continue

  dataset="$(basename "${dataset_dir}")"
  if [[ -n "${DATASET_FILTER_REGEX}" && ! "${dataset}" =~ ${DATASET_FILTER_REGEX} ]]; then
    continue
  fi

  h5="${dataset_dir}/counts_modisco_output.h5"
  report_dir="${dataset_dir}/counts_modisco_report"
  meme_out="${dataset_dir}/${dataset}_memedb.counts.txt"

  if [[ ! -f "${h5}" ]]; then
    echo "Skipping ${dataset}: missing ${h5}"
    continue
  fi

  mkdir -p "${report_dir}"

  if [[ "${REGENERATE_REPORT}" == "1" || ! -f "${report_dir}/motifs.html" ]]; then
    echo "Running modisco report for ${dataset}"
    modisco report \
      -i "${h5}" \
      -o "${report_dir}" \
      -s ./ \
      -m "${MEME_DB}" \
      -n 10
  else
    echo "Skipping modisco report for ${dataset}: motifs.html already exists"
  fi

  if [[ "${REGENERATE_MEME}" == "1" || ! -f "${meme_out}" ]]; then
    echo "Running modisco meme for ${dataset}"
    modisco meme \
      -i "${h5}" \
      -t PFM \
      -o "${meme_out}"
  else
    echo "Skipping modisco meme for ${dataset}: MEME export already exists"
  fi
done
