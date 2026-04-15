#!/usr/bin/bash

# Purpose: extract CWMs from merged MoDISco patterns in MEME format for downstream
# use in R and for visualization.

#SBATCH --output=../../logs/03-chrombpnet/02/04/04f-cwm_to_meme-%j.out
#SBATCH --partition=main
#SBATCH -t 04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

# fail explicitly for any errors, nonexistent variables, etc
# https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425#set--u
set -euo pipefail

source ../config.sh

# load conda env
eval "$(conda shell.bash hook)"
conda activate modiscolite

load_optional_module() {
  local mod="$1"
  [[ -n "${mod}" ]] || return 0
  if module -t avail "${mod}" 2>&1 | grep -Fq "${mod}"; then
    module load "${mod}"
  else
    echo "WARNING: module '${mod}' is unavailable; continuing without it."
  fi
}

for mod in system cairo; do
  load_optional_module "${mod}"
done

compiled_h5="${modisco_comp_dir}/modisco_compiled.h5"
echo $compiled_h5

python -u 04f-CWM_to_MEME.py --filename $compiled_h5 \
  --datatype CWM \
  --output "${modisco_comp_dir}/modisco_compiled.memedb.txt"

echo "done."
