#!/usr/bin/env bash
set -euo pipefail

mapfile -t scripts < <(
    find . -maxdepth 1 -type f -name "*.sh" \
        ! -name "00_run_preprocess.sh" \
        ! -name "01-run_seurat.sh" \
        ! -name "01-run.sh" \
        ! -name "._*" \
    | sed 's|^\./||' \
    | sort
)

for script in "${scripts[@]}"; do
    echo "====================================="
    echo "Running $script"
    echo "====================================="

    bash "$script"
done