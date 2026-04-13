#!/bin/bash
set -euo pipefail
#SBATCH --output=../../logs/03-chrombpnet/01/02/%x-%j.out
#SBATCH -p akundaje,gpu
#SBATCH -t 2-0
#SBATCH -c 4
#SBATCH --mem=60G
#SBATCH --gres=gpu:a40:1
#SBATCH --requeue
#SBATCH --open-mode=append


# PARSE ARGUMENTS --------------------------------------------------------------

frag_file="${1}"
ref_fasta="${2}"
chromsizes="${3}"
peaks_file="${4}"
negatives_file="${5}"
split_file="${6}"
bias_model="${7}"
out_dir="${8}"
data_type="${9}"

function timestamp {
    # Function to get the current time with the new line character removed 
    
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
for mod in system cairo pango; do
    load_optional_module "${mod}"
done

echo "--- $(timestamp): Beginning training ---"

# Sometimes the job starts then goes back to pending (e.g. if preempted), and then
# when it restarts, chrombpnet qc errors because the interpret_subsample dir is already
# created. have the jobscript clear the interpret_subsample dir before running
# chrombpnet qc again just in case
if [[ -d "${out_dir%/}/auxiliary/interpret_subsample" ]]; then
    echo "Found existing 'interpret_subsample' dir. Deleting previous subsample data..."
    rm -rf "${out_dir%/}/auxiliary/interpret_subsample"
fi


chrombpnet pipeline \
        --input-fragment-file "${frag_file}" \
        --genome "${ref_fasta}" \
        --chrom-sizes "${chromsizes}" \
        --peaks "${peaks_file}" \
        --nonpeaks "${negatives_file}" \
        --chr-fold-path "${split_file}" \
        --bias-model-path "${bias_model}" \
        --output-dir "${out_dir}" \
        --data-type "${data_type}"
        

echo "--- $(timestamp): Completed training ---"
# done.
