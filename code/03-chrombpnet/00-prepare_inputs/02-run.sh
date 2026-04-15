#!/usr/bin/bash
#SBATCH --job-name="02-process_frags"
#SBATCH --time=01:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=main
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

# source configuration variables
source ../config.sh

export LC_ALL=C
export LC_CTYPE=C
export LANG=C

input_parallel=8
datasets1=$(
  find "${sample_frags_dir}/fragments" -maxdepth 1 -type f -name '*.tsv' ! -name '._*' -print |
    while IFS= read -r path; do
      basename "${path}" .tsv
    done
)

datasets2=$( for dataset in ${datasets1[@]}; do

    # check if the out file exists
    # if not, keep the dataset & get list of all unique ones
    done_file="${sample_frags_dir}/pseudorepT/${dataset}.tsv"
    if [[ ! -f "$done_file" ]]; then
        echo "${dataset}.tsv"
    fi

done | uniq )


echo "@ processing fragments file: ${datasets2}"

for dataset in ${datasets2}; do
  python 02-process_frags.py "${dataset}" "${sample_frags_dir}" "${chromsizes}"
done

echo "@ done"
