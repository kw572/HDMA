#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/06/%x-%j.out
#SBATCH --partition=main
#SBATCH -t 2-0
#SBATCH --mem=50G
#SBATCH -C NO_GPU
#SBATCH --cpus-per-task=12

set -euo pipefail


peak_shaps="${1}"
max_seqlets="${2}"
num_leiden="${3}"
modisco_output="${4}"
report_outdir="${5}"
img_suffix_dir="${6}"
meme_db="${7}"
num_matches="${8}"
output_memedb="${9}"

# load conda environment
eval "$(conda shell.bash hook)"
conda activate modiscolite

load_optional_module() {
    local mod="$1"
    [[ -n "${mod}" ]] || return 0
    if module -t avail "${mod}" 2>&1 | grep -Fq "${mod}"; then
        module load "${mod}"
    else
        echo "WARNING: module '${mod}' is unavailable; continuing without it."
    fi
}

for mod in system libxml2 libxslt perl zlib ghostscript cairo; do
    load_optional_module "${mod}"
done

export PATH="$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.5:$PATH"

echo "[$(date +"%m/%d/%Y (%r)")] running modisco motifs..."
modisco motifs -i "${peak_shaps}" -n "${max_seqlets}" -l "${num_leiden}" -o "${modisco_output}"

echo "[$(date +"%m/%d/%Y (%r)")] running modisco report..."
modisco report -i "${modisco_output}" -o "${report_outdir}" -s "${img_suffix_dir}" -m "${meme_db}" -n "${num_matches}"

echo "[$(date +"%m/%d/%Y (%r)")] exporting modisco meme..."
modisco meme -i "${modisco_output}" -t PFM -o "${output_memedb}"
