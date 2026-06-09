#!/usr/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/02/04/04b-get_tomtom-%j.out
#SBATCH --partition=main
#SBATCH -t 01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

# fail explicitly for any errors, nonexistent variables, etc
# https://gist.github.com/mohanpedala/1e2ff5661761d3abd0385e8223e16425#set--u
# set -euo pipefail

source ../config.sh

# load conda env
eval "$(conda shell.bash hook)"
conda activate modiscolite

if [[ -z "${TOMTOM_EXEC_PATH:-}" && -x "${HOME}/.conda/envs/meme-suite/bin/tomtom" ]]; then
  export TOMTOM_EXEC_PATH="${HOME}/.conda/envs/meme-suite/bin/tomtom"
fi

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

compiled_h5=$modisco_comp_dir/modisco_compiled.h5
out_dir=$modisco_comp_dir
meme_db="${VIERSTRA_MEME_DB:-${vierstra_meme_db}}"

echo $compiled_h5
echo $out_dir

python -u 04b-get_tomtom_matches.py --modisco-h5 $compiled_h5 \
    --out-dir ${out_dir} \
    --meme-db ${meme_db} \
    --tomtom-exec "${TOMTOM_EXEC_PATH:-tomtom}" \
    --verbose True
                    
echo "done."
