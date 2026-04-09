#!/bin/bash
#SBATCH --job-name=12-annotate_hits
#SBATCH --output=../../logs/03-chrombpnet/00/12/%x-%j.out
#SBATCH --partition=akundaje,wjg,biochem,sfgf
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=6
#SBATCH --time=01:00:00

# Purpose: runner script to annotate peaks.
#
# Usage:
# $ organs=( Adrenal Brain Eye Heart Liver Lung Muscle Skin Spleen Stomach Thymus Thyroid )
# $ for i in ${organs[@]}; do sbatch 12-annotate_peaks.sh $i; sleep 1s; done

set -euo pipefail

source ../config.sh

export LC_ALL=C
export LC_CTYPE=C
export LANG=C

genome_id="${1:-GRCz11}"

datasets=$(
  find "${chrombpnet_peaks_dir}" -maxdepth 1 -type f -name '*__peaks_bpnet.narrowPeak' ! -name '._*' -print |
    while IFS= read -r path; do
      basename "${path}" __peaks_bpnet.narrowPeak
    done
)

datasets_to_do=$(
  for dataset in ${datasets}; do
    out_file="${anno_peaks_dir%/}/${dataset}__peaks_bpnet.annotated.tsv"
    if [[ ! -f "${out_file}" ]]; then
      echo "${dataset}"
    fi
  done
)

echo "@ Annotating cell types: ${datasets_to_do}"

for dataset in ${datasets_to_do}; do
  peaks_file="${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak"
  out_tsv="${anno_peaks_dir%/}/${dataset}__peaks_bpnet.annotated.tsv"
  echo "@ ${dataset}"
  Rscript 12-annotate_peaks.R "${peaks_file}" "${out_tsv}" "${genome_id}"
done
