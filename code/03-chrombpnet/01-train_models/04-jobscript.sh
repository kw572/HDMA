#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/04/%x-%j.out
#SBATCH -p akundaje,wjg,biochem,sfgf
#SBATCH --time=2-00:00:00
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH --gres=gpu:a40:1
#SBATCH --requeue
#SBATCH --open-mode=append

set -euo pipefail

# NOTE: for jobs that timeout with 2 days,
# switch, set a longer time limit and submit to akundaje partition only
# #SBATCH -p akundaje
# #SBATCH --time=3-00:00:00
# otherwise:
# #SBATCH -p akundaje,owners,gpu
# #SBATCH --time=2-00:00:00

ref_fasta="${1}"
peaks_file="${2}"
model_file="${3}"
out_prefix="${4}"
chrom_sizes="${5}"

function timestamp {
    # Function to get the current time with the new line character
    # removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

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

echo "--- $(timestamp): Beginning interpretation ---"

chrombpnet contribs_bw --genome "${ref_fasta}" \
                                  --regions "${peaks_file}" \
                                  --model-h5 "${model_file}" \
                                  --output-prefix "${out_prefix}" \
                                  --chrom-sizes "${chrom_sizes}"

echo "--- $(timestamp): Completed interpretation ---"
