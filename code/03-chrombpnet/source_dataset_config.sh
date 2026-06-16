#!/usr/bin/env bash

set -euo pipefail

script_dir="${1:-$(pwd)}"
repo_root="$(cd "${script_dir}/.." && pwd)"
root_config="${repo_root}/config.sh"
candidate=""

if [[ -n "${CHROMBPNET_CONFIG:-}" ]]; then
  candidate="${CHROMBPNET_CONFIG}"
  if [[ "${candidate}" == "${root_config}" ]]; then
    echo "ERROR: CHROMBPNET_CONFIG points to repo root config.sh; use a dataset-specific config instead." >&2
    exit 1
  fi
  if [[ ! -f "${candidate}" ]]; then
    echo "ERROR: CHROMBPNET_CONFIG does not exist: ${candidate}" >&2
    exit 1
  fi
elif [[ -n "${CHROMBPNET_DATA_ROOT:-}" ]]; then
  candidate="${CHROMBPNET_DATA_ROOT%/}/config.sh"
  if [[ ! -f "${candidate}" ]]; then
    echo "ERROR: dataset-specific config missing under CHROMBPNET_DATA_ROOT: ${candidate}" >&2
    exit 1
  fi
else
  echo "ERROR: set CHROMBPNET_DATA_ROOT or CHROMBPNET_CONFIG to a dataset-specific config before running training scripts." >&2
  exit 1
fi

export CHROMBPNET_CONFIG="${candidate}"
source "${candidate}"
