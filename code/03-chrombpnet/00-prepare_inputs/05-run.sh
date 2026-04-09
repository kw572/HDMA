#!/usr/bin/bash
#SBATCH --job-name="05-format_peaks"
#SBATCH --time=01:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=akundaje
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

out_dir=$chrombpnet_peaks_dir
export LC_ALL=C
export LC_CTYPE=C
export LANG=C

datasets=$(
  find "${peaks_dir}" -maxdepth 1 -type f -name '*__peaks_overlap_filtered.narrowPeak' ! -name '._*' -print |
    while IFS= read -r path; do
      basename "${path}" __peaks_overlap_filtered.narrowPeak
    done
)

echo "@ formatting peaks for clusters: ${datasets}"

# DEBUG:
# python 05-format_peaks.py Adrenal_c0 ${peaks_dir} ${out_dir}

for dataset in ${datasets}; do
  python 05-format_peaks.py "${dataset}" "${peaks_dir}" "${out_dir}" "${chromsizes}"
done

echo "@ done"
