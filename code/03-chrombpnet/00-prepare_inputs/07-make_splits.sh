#!/bin/bash
#SBATCH --job-name="make_splits"
#SBATCH --time=01:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
# This file creates the ChromBPNet chromosomal folds files specifying train/test/validation
# folds, using the splits in https://zenodo.org/records/7445373

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# get configuration variables
source ../config.sh

chrombpnet_bin="${CHROMBPNET_BIN:-/home1/kuangtse/.conda/envs/chrombpnet/bin/chrombpnet}"

chroms=()
while IFS= read -r chrom; do
  [[ -n "${chrom}" ]] || continue
  chroms+=("${chrom}")
done < <(awk 'NF {print $1}' "${chromsizes}")
num_folds=5

if (( ${#chroms[@]} < num_folds + 2 )); then
  echo "Not enough chromosomes in ${chromsizes} to build ${num_folds} folds" >&2
  exit 1
fi

declare -a buckets
for ((i=0; i<${#chroms[@]}; i++)); do
  bucket_idx=$((i % num_folds))
  buckets[$bucket_idx]="${buckets[$bucket_idx]:-} ${chroms[$i]}"
done

for ((fold=0; fold<num_folds; fold++)); do
  test_idx=$fold
  valid_idx=$(((fold + 1) % num_folds))

  read -r -a test_chroms <<< "${buckets[$test_idx]}"
  read -r -a valid_chroms <<< "${buckets[$valid_idx]}"

  "${chrombpnet_bin}" prep splits \
    -c "${chromsizes}" \
    --test-chroms "${test_chroms[@]}" \
    --valid-chroms "${valid_chroms[@]}" \
    -op "${split_dir}/fold_${fold}"
done
