#!/usr/bin/bash
#SBATCH --job-name="10-subset_background"
#SBATCH --time=03:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=main
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

export LC_ALL=C
export LC_CTYPE=C
export LANG=C
num_folds="${chrombpnet_train_num_folds}"

datasets=$(
  find "${peaks_dir}" -maxdepth 1 -type f -name '*__peaks_overlap_filtered.narrowPeak' ! -name '._*' -print |
    while IFS= read -r path; do
      basename "${path}" __peaks_overlap_filtered.narrowPeak
    done
)

echo "@ formatting peaks for clusters: ${datasets}"

# DEBUG:
for dataset in ${datasets}; do
  peak_file="${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak"
  if [[ ! -s "${peak_file}" ]]; then
    echo "@ ${dataset} has no chrombpnet peaks; skipping"
    continue
  fi

  python 10-subset_backgrounds.py \
    --cluster "${dataset}" \
    --negatives-dir "${negatives_dir}" \
    --output-dir "${negatives_subset_dir}" \
    --num-folds "${num_folds}"
  
done

echo "@ done"
