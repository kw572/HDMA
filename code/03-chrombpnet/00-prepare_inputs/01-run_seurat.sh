#!/usr/bin/bash
#SBATCH --job-name="01-frags-rds"
#SBATCH --time=04:00:00
#SBATCH --output=../../logs/03-chrombpnet/00/%x-%j.out
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

set -euo pipefail

if [[ -f "../../SJ_RENV_DIR.txt" ]]; then
  renv_dir=""
  IFS= read -r renv_dir < "../../SJ_RENV_DIR.txt"
  export R_LIBS_USER=$renv_dir
fi

source ../config.sh

MAPPING_TSV="${sample_frags_dir}/barcode_to_output_stem.tsv"
MANIFEST_TSV="${sample_frags_dir}/fragments_manifest.tsv"
STATS_TSV="${sample_frags_dir}/fragments_stats.tsv"

R_CMD=(
  Rscript 01a-export_rds_metadata.R
  --rds "$rds_file"
  --mapping-out "$MAPPING_TSV"
  --manifest-out "$MANIFEST_TSV"
  --rds "$rds_file"
  --cluster-col "$cluster_col"
)

if [[ -n "${fragment_file:-}" ]]; then
  R_CMD+=(--fragments "$fragment_file")
fi

if [[ -n "${sample_col:-}" ]]; then
  R_CMD+=(--sample-col "$sample_col")
fi

if [[ -n "${sample_name:-}" ]]; then
  R_CMD+=(--sample-name "$sample_name")
fi

if [[ -n "${assay_name:-}" ]]; then
  :
fi

"${R_CMD[@]}"

PY_CMD=(
  python 01b-route_fragments.py
  --fragments "$fragment_file"
  --mapping "$MAPPING_TSV"
  --chromsizes "$chromsizes"
  --outdir "${sample_frags_dir}/fragments"
  --stats-out "$STATS_TSV"
)

"${PY_CMD[@]}"
