#!/usr/bin/bash
#SBATCH --job-name="09-idx_peaks"
#SBATCH --time=02:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=main
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G

# Compress and sort peaks for viewing in the WashU browser.

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

export LC_ALL=C
export LC_CTYPE=C
export LANG=C

peaksets=$(
    find "${chrombpnet_peaks_dir}" -maxdepth 1 -type f -name '*.narrowPeak' ! -name '._*' -print
)
echo "${peaksets}"

idx_pk () {

    peakset=$1
    tmp_sorted="${peakset}.sorted"
    bedtools sort -faidx "${chromsizes}" -i "${peakset}" > "${tmp_sorted}"
    mv "${tmp_sorted}" "${peakset}"
    bgzip -c "${peakset}" > "${peakset}.gz"
    tabix -p bed "${peakset}.gz"

}
export -f idx_pk

while IFS= read -r peakset; do
    [[ -n "${peakset}" ]] || continue
    idx_pk "${peakset}"
done <<< "${peaksets}"
