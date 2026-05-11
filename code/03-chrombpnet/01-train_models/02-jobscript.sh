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

source ../config.sh

dataset="${CHROMBPNET_TRAIN_DATASET:?CHROMBPNET_TRAIN_DATASET is required}"
fold_name="${CHROMBPNET_TRAIN_FOLD:?CHROMBPNET_TRAIN_FOLD is required}"
bias_params="${CHROMBPNET_TRAIN_BIAS_PARAMS:-${chrombpnet_bias_params}}"
out_dir="${CHROMBPNET_TRAIN_OUT_DIR:-${models_dir%/}/bias_${bias_params}/${dataset}/${fold_name}}"

frag_file="${cluster_frags_dir%/}/fragments/${dataset}__sorted.tsv"
peaks_file="${chrombpnet_peaks_dir%/}/${dataset}__peaks_bpnet.narrowPeak"
negatives_file="${negatives_dir%/}/${dataset}/${fold_name}/output_negatives.bed"
split_file="${split_dir%/}/${fold_name}.json"
bias_model="${bias_dir%/}/${bias_params}/models/bias.h5"
data_type="${chrombpnet_data_type}"

timestamp() {
    date +"%Y-%m-%d_%H-%M-%S" | tr -d '\n'
}

export LANG=C
export LC_ALL=C

set +u
source /etc/profile

module purge
module load legacy/CentOS7
module load gcc/8.3.0
module load cuda/11.2.0
module load cudnn/8.1.0.77-11.2-cuda

# optional graphics deps
module load cairo || true
module load pango || true

eval "$(conda shell.bash hook)"
conda activate chrombpnet
set -u

echo "=== NVIDIA ==="
nvidia-smi

echo "=== TENSORFLOW GPU CHECK ==="
python - <<'EOF'
import tensorflow as tf
print("TF version:", tf.__version__)
print("GPUs:", tf.config.list_physical_devices('GPU'))
EOF

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

negative_sampling_ratio="${chrombpnet_negative_sampling_ratio}"
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

echo "Running chrombpnet pipeline"
echo "--input-fragment-file ${frag_file}"
echo "--genome ${ref_fasta}"
echo "--chrom-sizes ${chromsizes}"
echo "--peaks ${peaks_file}"
echo "--nonpeaks ${prepared_negatives_file}"
echo "--chr-fold-path ${split_file}"
echo "--bias-model-path ${bias_model}"
echo "--output-dir ${out_dir}"
echo "--data-type ${data_type}"


python - "${frag_file}" "${ref_fasta}" "${chromsizes}" "${peaks_file}" "${prepared_negatives_file}" "${split_file}" "${bias_model}" "${out_dir}" "${data_type}" <<'PY'
import sys

from chrombpnet import CHROMBPNET
from chrombpnet.training.data_generators import initializers

orig_fetch = initializers.fetch_data_and_model_params_based_on_mode


def patched_fetch(mode, args, parameters, nonpeak_regions, peak_regions):
    if mode != "valid" or nonpeak_regions is None or peak_regions is None:
        return orig_fetch(mode, args, parameters, nonpeak_regions, peak_regions)

    inputlen = int(parameters["inputlen"])
    outputlen = int(parameters["outputlen"])
    requested = int(float(parameters["negative_sampling_ratio"]) * peak_regions.shape[0])
    available = nonpeak_regions.shape[0]
    sample_n = min(requested, available)

    if sample_n < requested:
        print(
            f"Validation negatives capped from {requested} to {sample_n} "
            f"for split '{mode}' because only {available} are available after filtering."
        )

    nonpeak_regions = nonpeak_regions.sample(
        n=sample_n,
        replace=False,
        random_state=args.seed,
    )

    negative_sampling_ratio = 1.0
    max_jitter = 0
    add_revcomp = False
    shuffle_at_epoch_start = False

    return (
        inputlen,
        outputlen,
        nonpeak_regions,
        negative_sampling_ratio,
        max_jitter,
        add_revcomp,
        shuffle_at_epoch_start,
    )


initializers.fetch_data_and_model_params_based_on_mode = patched_fetch

(
    frag_file,
    ref_fasta,
    chromsizes,
    peaks_file,
    prepared_negatives_file,
    split_file,
    bias_model,
    out_dir,
    data_type,
) = sys.argv[1:]

sys.argv = [
    "chrombpnet",
    "pipeline",
    "--input-fragment-file",
    frag_file,
    "--genome",
    ref_fasta,
    "--chrom-sizes",
    chromsizes,
    "--peaks",
    peaks_file,
    "--nonpeaks",
    prepared_negatives_file,
    "--chr-fold-path",
    split_file,
    "--bias-model-path",
    bias_model,
    "--output-dir",
    out_dir,
    "--data-type",
    data_type,
]

CHROMBPNET.main()
PY
        

echo "--- $(timestamp): Completed training ---"
# done.
