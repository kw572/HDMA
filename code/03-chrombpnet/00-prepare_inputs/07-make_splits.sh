#!/usr/bin/env bash
#SBATCH --job-name="make_splits"
#SBATCH --time=01:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G

set -euo pipefail

# Create ChromBPNet chromosome fold JSONs from the chromosomes present in config.sh.

set +u
eval "$(conda shell.bash hook)"
conda activate chrombpnet
set -u

source ../config.sh

chroms=()
while IFS= read -r chrom; do
  [[ -n "${chrom}" ]] || continue
  chroms+=("${chrom}")
done < <(awk 'NF {print $1}' "${chromsizes}")

num_folds="${chrombpnet_train_num_folds}"
if (( ${#chroms[@]} < num_folds + 2 )); then
  echo "Not enough chromosomes in ${chromsizes} to build ${num_folds} folds" >&2
  exit 1
fi

declare -a buckets
for ((i = 0; i < ${#chroms[@]}; i++)); do
  bucket_idx=$((i % num_folds))
  buckets[$bucket_idx]="${buckets[$bucket_idx]:-} ${chroms[$i]}"
done

for ((fold = 0; fold < num_folds; fold++)); do
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
