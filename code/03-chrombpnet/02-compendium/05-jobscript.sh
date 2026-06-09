#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/02/05/%x-%j.out
#SBATCH --partition=gpu
#SBATCH -t 24:00:00
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH --gres=gpu:1
#SBATCH --requeue
#SBATCH --open-mode=append

set -eo pipefail

# get args
celltype=${1}
peaks_bed=${2}
modisco_h5=${3}
shaps_h5=${4}
finemo_out=${5}
finemo_out_nocompo=${6}
shap_type=${7}
alpha=${8}
alpha_nocompo=${9}
motifs_nocompo=${10}

echo -e "\t@ ${celltype}"
echo -e "\t@ ${peaks_bed}"
echo -e "\t@ ${modisco_h5}"
echo -e "\t@ ${shaps_h5}"
echo -e "\t@ ${finemo_out}"
echo -e "\t@ ${finemo_out_nocompo}"
echo -e "\t@ ${alpha}"
echo -e "\t@ ${alpha_nocompo}"
echo -e "\t@ ${motifs_nocompo}"

# Prefer an explicit virtualenv on HPC, but keep the legacy conda path as a fallback.
FINEMO_VENV="${FINEMO_VENV:-/scratch1/${USER}/venvs/finemo310}"
if [[ ! -x "${FINEMO_VENV}/bin/activate" && -x "/scratch1/${USER}/venvs/finemo/bin/activate" ]]; then
  FINEMO_VENV="/scratch1/${USER}/venvs/finemo"
elif [[ ! -x "${FINEMO_VENV}/bin/activate" && -x "$HOME/venvs/finemo/bin/activate" ]]; then
  FINEMO_VENV="$HOME/venvs/finemo"
fi
if [[ -x "${FINEMO_VENV}/bin/activate" ]]; then
  # shellcheck disable=SC1090
  source "${FINEMO_VENV}/bin/activate"
else
  # NOTE: this call to conda triggers an unbound variable error, hence we 
  # use set -eo pipefail instead of -euo pipefail
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

# print version for debugging
pip freeze | grep finemo


# RUN FINEMO

i=1
outs=( "$finemo_out" "$finemo_out_nocompo" )
for out in ${outs[@]}; do

  echo ${out}

  finemo_npz="${out%/}/intermediate_inputs.npz"

  if [[ -f "${out%/}/hits.bed.gz" ]]; then

    echo "@ found hits."

  else

    echo "[$(date +"%m/%d/%Y (%r)")] starting ${celltype}"
    echo "@ extract-regions-h5"
    finemo extract-regions-chrombpnet-h5 --h5s ${shaps_h5} --out-path ${finemo_npz} --region-width 1000

    echo "[$(date +"%m/%d/%Y (%r)")]"
    echo "@ call-hits"

    # the first time, use all motifs
    if [ $i -eq 1 ]; then
      finemo call-hits -r ${finemo_npz} -m ${modisco_h5} -p ${peaks_bed} --a ${alpha} -o ${out} -b 200
    # the second time, use non-composite motifs only
    elif [ $i -eq 2 ]; then
      finemo call-hits -r ${finemo_npz} -m ${modisco_h5} -p ${peaks_bed} --a ${alpha_nocompo} -o ${out} -b 200 --motifs-include ${motifs_nocompo}
    fi

    # make report
    finemo report -r ${finemo_npz} -H ${out}/hits.tsv -p ${peaks_bed} -m ${modisco_h5} -o ${out} --no-recall --no-seqlets

    echo "[$(date +"%m/%d/%Y (%r)")]"

    # INDEX HITS -----------------------------------------------------------------
    if [[ -f "${out}/hits.bed" ]]; then

    	echo "@ found hits; renaming and indexing"

    	# bgzip and index hits
    	module load biology samtools

    	bgzip -c ${out}/hits.bed > ${out}/hits.bed.gz
    	tabix -p bed ${out}/hits.bed.gz

    fi

  fi

  ((i++))

done

echo "[$(date +"%m/%d/%Y (%r)")] done"
