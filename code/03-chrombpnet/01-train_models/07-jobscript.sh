#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/07/%x-%j.out
#SBATCH -p akundaje,owners,wjg,sfgf,biochem
#SBATCH -t 01:00:00
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -C NO_GPU

set -euo pipefail

celltype="${1}"
counts_shaps="${2}"
counts_out="${3}"
profile_shaps="${4}"
profile_out="${5}"
peaks_file="${6}"
chrom_sizes="${7}"

echo -e "@ Using counts shaps: ${counts_shaps} --> ${counts_out}.bw"
echo -e "@ Using profile shaps: ${profile_shaps} --> ${profile_out}.bw"

# source configuration variables
source ../config.sh

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# locate the helper script from the installed chrombpnet package
script_loc=$(python - <<'PY'
import chrombpnet.evaluation.make_bigwigs.importance_hdf5_to_bigwig as m
print(m.__file__)
PY
)

# generate bigwigs if they don't exist yet
if [[ ! -f "${counts_out}.bw" ]]; then
  echo "[$(date +"%m/%d/%Y (%r)")] starting ${celltype} - generating bigwig for counts"
  echo ${counts_shaps}
  echo ${peaks_file}
  echo ${counts_out}
  python "${script_loc}" --hdf5 "${counts_shaps}" \
                         --regions "${peaks_file}" \
                         --chrom-sizes "${chrom_sizes}" \
                         --output-prefix "${counts_out}"
fi

if [[ ! -f "${profile_out}.bw" ]]; then
  echo "[$(date +"%m/%d/%Y (%r)")] starting ${celltype} - generating bigwig for profile"
  python "${script_loc}" --hdf5 "${profile_shaps}" \
                         --regions "${peaks_file}" \
                         --chrom-sizes "${chrom_sizes}" \
                         --output-prefix "${profile_out}"
fi

echo "[$(date +"%m/%d/%Y (%r)")] done!"
