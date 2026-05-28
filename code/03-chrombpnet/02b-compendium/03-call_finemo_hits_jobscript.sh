#!/bin/bash
#SBATCH --output=Log/%x-%j.out
#SBATCH --partition=gpu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --open-mode=append

set -eo pipefail

dataset="${1}"
peaks_bed="${2}"
modisco_h5="${3}"
shaps_h5="${4}"
out_dir="${5}"
compendium_dir="${6}"
alpha="${7}"

RUN_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
SCRIPT_DIR="${RUN_DIR}"
mkdir -p "${RUN_DIR}/Log" "${out_dir}"

FINEMO_VENV="${FINEMO_VENV:-/scratch1/${USER}/venvs/finemo310}"
if [[ ! -x "${FINEMO_VENV}/bin/activate" && -x "/scratch1/${USER}/venvs/finemo/bin/activate" ]]; then
  FINEMO_VENV="/scratch1/${USER}/venvs/finemo"
fi
if [[ -x "${FINEMO_VENV}/bin/activate" ]]; then
  # shellcheck disable=SC1090
  source "${FINEMO_VENV}/bin/activate"
else
  eval "$(conda shell.bash hook)"
  conda activate finemo
fi

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

echo "@ dataset ${dataset}"
echo "@ peaks ${peaks_bed}"
echo "@ modisco_h5 ${modisco_h5}"
echo "@ shaps_h5 ${shaps_h5}"
echo "@ out_dir ${out_dir}"
echo "@ compendium_dir ${compendium_dir}"
echo "@ alpha ${alpha}"

finemo_npz="${out_dir%/}/intermediate_inputs.npz"

if [[ ! -f "${out_dir%/}/hits.tsv" ]]; then
  finemo extract-regions-chrombpnet-h5 \
    --h5s "${shaps_h5}" \
    --out-path "${finemo_npz}" \
    --region-width 1000

  finemo call-hits \
    -r "${finemo_npz}" \
    -m "${modisco_h5}" \
    -p "${peaks_bed}" \
    --a "${alpha}" \
    -o "${out_dir}" \
    -b 200

  finemo report \
    -r "${finemo_npz}" \
    -H "${out_dir}/hits.tsv" \
    -p "${peaks_bed}" \
    -m "${modisco_h5}" \
    -o "${out_dir}" \
    --no-recall \
    --no-seqlets
fi

if [[ -f "${out_dir}/hits.bed" && ! -f "${out_dir}/hits.bed.gz" ]]; then
  load_optional_module biology
  load_optional_module samtools
  bgzip -c "${out_dir}/hits.bed" > "${out_dir}/hits.bed.gz"
  tabix -p bed "${out_dir}/hits.bed.gz"
fi

python "${SCRIPT_DIR}/03-map_finemo_to_compendium.py" \
  --compendium-dir "${compendium_dir}" \
  --dataset "${dataset}" \
  --hits-tsv "${out_dir}/hits.tsv" \
  --output-tsv-gz "${out_dir}/hits.compendium.tsv.gz" \
  --summary-json "${out_dir}/hits.compendium.summary.json"

echo "[$(date +"%m/%d/%Y (%r)")] done"
