#!/bin/bash
#SBATCH --job-name=chrombpnet-crested-compendium
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=Log/%x-%j.out
#SBATCH --error=Log/%x-%j.err

set -euo pipefail

if ! command -v module >/dev/null 2>&1; then
  set +u
  [[ -f /etc/profile ]] && source /etc/profile
  [[ -f /usr/share/lmod/lmod/init/bash ]] && source /usr/share/lmod/lmod/init/bash
  [[ -f /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash
  set -u
fi

source ~/.bashrc 2>/dev/null || true

RUN_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "${RUN_DIR}"

resolve_script_dir() {
  local candidate
  for candidate in \
    "${RUN_DIR}" \
    "${RUN_DIR}/code/03-chrombpnet/02b-compendium" \
    "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" 2>/dev/null && pwd)"
  do
    if [[ -n "${candidate}" && -f "${candidate}/../config.sh" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

SCRIPT_DIR="$(resolve_script_dir)" || {
  echo "ERROR: could not locate 02b-compendium directory from RUN_DIR=${RUN_DIR}" >&2
  exit 1
}
cd "${SCRIPT_DIR}"

mkdir -p Log

bash "${SCRIPT_DIR}/01-build_crested_compendium.sh"
bash "${SCRIPT_DIR}/02-annotate_and_heatmap.sh"
