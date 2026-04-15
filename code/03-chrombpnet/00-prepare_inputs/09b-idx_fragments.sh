#!/usr/bin/bash
#SBATCH --job-name="10b-idx_frags"
#SBATCH --time=06:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=main
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G

# Compress and sort peaks for space saving.

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

export LC_ALL=C
export LC_CTYPE=C
export LANG=C

fragsets=$(
    find "${cluster_frags_dir}/fragments" -maxdepth 1 -type f -name '*__sorted.tsv' ! -name '._*' -print
)
echo "${fragsets}"

idx_frag () {

    fragset=$1
    bgzip -c "${fragset}" > "${fragset}.gz"
    tabix -p bed "${fragset}.gz"

}
export -f idx_frag

while IFS= read -r fragset; do
    [[ -n "${fragset}" ]] || continue
    idx_frag "${fragset}"
done <<< "${fragsets}"
