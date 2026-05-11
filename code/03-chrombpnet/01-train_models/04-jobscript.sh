#!/bin/bash
#SBATCH --output=../../logs/03-chrombpnet/01/04/%x-%j.out
#SBATCH --partition=gpu
#SBATCH --time=2-00:00:00
#SBATCH -c 8
#SBATCH --mem=100G
#SBATCH --gres=gpu:1
#SBATCH --requeue
#SBATCH --open-mode=append

set -eo pipefail

ref_fasta="${1}"
peaks_file="${2}"
model_file="${3}"
out_prefix="${4}"
chrom_sizes="${5}"

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

echo "--- $(timestamp): Beginning interpretation ---"

chrombpnet contribs_bw \
    --genome "${ref_fasta}" \
    --regions "${peaks_file}" \
    --model-h5 "${model_file}" \
    --output-prefix "${out_prefix}" \
    --chrom-sizes "${chrom_sizes}"

echo "--- $(timestamp): Completed interpretation ---"