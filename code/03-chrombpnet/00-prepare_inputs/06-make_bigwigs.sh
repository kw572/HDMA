#!/usr/bin/env bash
#SBATCH --job-name="06-make_bw"
#SBATCH --time=12:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=main
#SBATCH --cpus-per-task=6
#SBATCH --mem=30G

set -euo pipefail

# Create unnormalized cut-site bigWigs from each cluster's sorted fragment file.

set +u
eval "$(conda shell.bash hook)"
conda activate chrombpnet
set -u

source ../config.sh

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
python_script="${script_dir}/06-make_bigwigs.py"
fragments_dir="${cluster_frags_dir%/}/fragments"
out_dir="${bigwigs_dir%/}"

datasets=$(
  find "${fragments_dir}" -maxdepth 1 -type f -name '*__sorted.tsv' ! -name '._*' -print |
    while IFS= read -r path; do
      dataset="$(basename "${path}" __sorted.tsv)"
      done_file="${out_dir}/${dataset}_unstranded.bw"
      if [[ ! -f "${done_file}" ]]; then
        printf '%s\n' "${dataset}"
      fi
    done
)

if [[ -z "${datasets}" ]]; then
  echo "@ all cluster bigwigs already exist in ${out_dir}"
  exit 0
fi

echo "@ Producing bigwigs for: ${datasets}"

for dataset in ${datasets}; do
  frag_file="${fragments_dir}/${dataset}__sorted.tsv"
  out_prefix="${out_dir}/${dataset}"
  echo "@ ${dataset}: ${frag_file} -> ${out_prefix}_unstranded.bw"
  python "${python_script}" \
    --fragments "${frag_file}" \
    --chromsizes "${chromsizes}" \
    --output-prefix "${out_prefix}"
done

echo "@ done"
