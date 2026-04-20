#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/02/%x-%j.out
#SBATCH --partition=gpu
#SBATCH -t 2-0
#SBATCH -c 4
#SBATCH --mem=60G
#SBATCH --gres=gpu:1
#SBATCH --requeue
#SBATCH --open-mode=append

set -euo pipefail


# PARSE ARGUMENTS --------------------------------------------------------------

frag_file="${1}"
ref_fasta="${2}"
chromsizes="${3}"
peaks_file="${4}"
negatives_file="${5}"
split_file="${6}"
bias_model="${7}"
out_dir="${8}"
data_type="${9}"

function timestamp {
    # Function to get the current time with the new line character removed 
    
    # current time
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

# load conda environment
eval "$(conda shell.bash hook)"
conda activate chrombpnet

load_optional_module() {
    local mod="$1"
    [[ -n "${mod}" ]] || return 0
    if module -t avail "${mod}" 2>&1 | grep -Fq "${mod}"; then
        module load "${mod}"
    else
        echo "WARNING: module '${mod}' is unavailable; continuing without it."
    fi
}

load_requested_modules() {
    local modules_string="$1"
    local mod
    [[ -n "${modules_string}" ]] || return 0
    read -r -a modules <<< "${modules_string}"
    for mod in "${modules[@]}"; do
        load_optional_module "${mod}"
    done
}

PRE_MODULES="${PRE_MODULES:-legacy/CentOS7 gcc/8.3.0}"
CUDA_MODULE="${CUDA_MODULE:-cuda/11.2.0}"
CUDNN_MODULE="${CUDNN_MODULE:-cudnn/8.1.0.77-11.2}"
load_requested_modules "${PRE_MODULES}"
load_optional_module "${CUDA_MODULE}"
load_optional_module "${CUDNN_MODULE}"
for mod in system cairo pango; do
    load_optional_module "${mod}"
done

echo "--- $(timestamp): Beginning training ---"

prepare_negatives_file() {
    local original_negatives="$1"
    local peaks_bed="$2"
    local split_json="$3"
    local ratio="$4"
    local tmp_root="${5:-${TMPDIR:-/tmp}}"

    python - "$original_negatives" "$peaks_bed" "$split_json" "$ratio" "$tmp_root" <<'PY'
import json
import math
import os
import random
import shutil
import sys
from collections import defaultdict
from pathlib import Path

negatives_path = Path(sys.argv[1])
peaks_path = Path(sys.argv[2])
split_path = Path(sys.argv[3])
ratio = float(sys.argv[4])
tmp_root = Path(sys.argv[5])

split_to_chroms = json.loads(split_path.read_text())
chrom_to_split = {
    chrom: split_name
    for split_name, chroms in split_to_chroms.items()
    for chrom in chroms
}

peak_counts = defaultdict(int)
with peaks_path.open() as handle:
    for line in handle:
        chrom = line.split("\t", 1)[0]
        split_name = chrom_to_split.get(chrom)
        if split_name is not None:
            peak_counts[split_name] += 1

negative_rows = defaultdict(list)
with negatives_path.open() as handle:
    for raw_line in handle:
        chrom = raw_line.split("\t", 1)[0]
        split_name = chrom_to_split.get(chrom)
        if split_name is not None:
            negative_rows[split_name].append(raw_line)

extras = []
for split_name, chroms in split_to_chroms.items():
    required = math.floor(ratio * peak_counts.get(split_name, 0))
    available = len(negative_rows.get(split_name, []))
    if available >= required:
        continue
    if available == 0 and required > 0:
        print(
            f"ERROR: split '{split_name}' has 0 negatives but requires {required} "
            f"for ratio {ratio}",
            file=sys.stderr,
        )
        sys.exit(1)
    needed = required - available
    print(
        f"Padding negatives for split '{split_name}': have {available}, need {required} "
        f"(adding {needed} duplicate rows)",
        file=sys.stderr,
    )
    rng = random.Random(0)
    rows = negative_rows[split_name]
    extras.extend(rng.choice(rows) for _ in range(needed))

if not extras:
    print(str(negatives_path))
    sys.exit(0)

tmp_root.mkdir(parents=True, exist_ok=True)
prepared_path = tmp_root / f"{negatives_path.stem}.padded_{os.getpid()}.bed"
shutil.copyfile(negatives_path, prepared_path)
with prepared_path.open("a") as handle:
    for row in extras:
        handle.write(row)

print(str(prepared_path))
PY
}

# ChromBPNet creates these directories with exist_ok=False. If a prior run died
# after partially populating the fold directory, remove the stale stage outputs
# so reruns can start cleanly without tripping on FileExistsError.
if [[ ! -f "${out_dir%/}/evaluation/overall_report.html" ]]; then
    for stale_dir in logs auxiliary models evaluation; do
        if [[ -d "${out_dir%/}/${stale_dir}" ]]; then
            echo "Found stale '${stale_dir}' dir from an incomplete run. Removing it before restart..."
            rm -rf "${out_dir%/}/${stale_dir}"
        fi
    done
fi

negative_sampling_ratio="${CHROMBPNET_NEGATIVE_SAMPLING_RATIO:-0.1}"
prepared_negatives_file="$(
    prepare_negatives_file \
        "${negatives_file}" \
        "${peaks_file}" \
        "${split_file}" \
        "${negative_sampling_ratio}" \
        "${TMPDIR:-/tmp}"
)"

if [[ "${prepared_negatives_file}" != "${negatives_file}" ]]; then
    echo "Using padded negatives file: ${prepared_negatives_file}"
fi


chrombpnet pipeline \
        --input-fragment-file "${frag_file}" \
        --genome "${ref_fasta}" \
        --chrom-sizes "${chromsizes}" \
        --peaks "${peaks_file}" \
        --nonpeaks "${prepared_negatives_file}" \
        --chr-fold-path "${split_file}" \
        --bias-model-path "${bias_model}" \
        --output-dir "${out_dir}" \
        --data-type "${data_type}"
        

echo "--- $(timestamp): Completed training ---"
# done.
