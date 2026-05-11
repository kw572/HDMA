#!/usr/bin/env bash
set -euo pipefail

export COPYFILE_DISABLE=1

cleanup_macos_sidecars() {
  local data_root=""
  if [[ -n "${CHROMBPNET_DATA_ROOT:-}" ]]; then
    data_root="${CHROMBPNET_DATA_ROOT}"
  else
    data_root="$(cd ../data && pwd)"
  fi

  if [[ -d "${data_root}" ]]; then
    find "${data_root}" -name '._*' -type f -delete 2>/dev/null || true
  fi
}

scripts=(
  "01-run_seurat.sh"
  "02-run.sh"
  "03a-cat_frags.sh"
  "03c-sort_frags.sh"
  "04-call_peaks.sh"
  "04b-make_signal_bw.sh"
  "05-run.sh"
  "06-make_bigwigs.sh"
  "07-make_splits.sh"
  "08-write_cmds.sh"
  "08-run.sh"
  "09-idx_peaks.sh"
  "09b-idx_fragments.sh"
  "10-run.sh"
  "11-count_peaks.sh"
  "12-annotate_peaks.sh"
)

for script in "${scripts[@]}"; do
    [[ -f "${script}" ]] || {
        echo "Missing preprocess script: ${script}" >&2
        exit 1
    }

    cleanup_macos_sidecars

    echo "====================================="
    echo "Running $script"
    echo "====================================="

    bash "$script"
done

cleanup_macos_sidecars
