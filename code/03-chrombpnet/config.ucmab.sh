#!/usr/bin/env bash

# Local config for running code/03-chrombpnet/00-prepare_inputs on this Mac.
# Edit the input file paths and metadata column names below to match your dataset.

set -euo pipefail

script_source="${BASH_SOURCE[0]:-$0}"
chrombpnet_root="$(cd "$(dirname "${script_source}")" && pwd)"
data_root="${chrombpnet_root}/data/local_run_chr"

# Core input files
rds_file="${data_root}/inputs/ucmab_14dpf.rds"
fragment_file="${data_root}/inputs/fragments.chr.tsv.gz"

# Metadata columns in the Seurat object
cluster_col="predicted.id"
sample_col=""
sample_name=""
assay_name="peaks"

# Reference files
refs="${data_root}/../local_run_chr/refs"
ref_fasta="${refs}/Danio_rerio.GRCz11.chr.dna.toplevel.fa"
chromsizes="${refs}/zf.chr.major.chrom.sizes"
blacklist="${refs}/blacklist.GRCz11.chr.bed"
genome_size=""
chrombpnet_bias_params="1-col_aspn_ogna_thresh0.4"
chrombpnet_train_num_folds=5
# Working directories
base_dir="${data_root}/work"
sample_frags_dir="${base_dir}/00-inputs/sample_fragments"
cluster_frags_dir="${base_dir}/00-inputs/cluster_fragments"
peaks_dir="${base_dir}/00-inputs/peaks"
chrombpnet_peaks_dir="${base_dir}/00-inputs/chrombpnet_peaks"
split_dir="${base_dir}/00-inputs/splits"
negatives_dir="${base_dir}/00-inputs/negatives"
negatives_subset_dir="${base_dir}/00-inputs/negatives_subset"
anno_peaks_dir="${base_dir}/00-inputs/annotated_peaks"

# Optional directories used by later steps.
models_dir="${base_dir}/01-models"
bias_dir="${base_dir}/01-models/bias"
qc_dir="${base_dir}/01-models/qc"
contribs_dir="${base_dir}/01-models/contribs"
contribs_scratch="${contribs_dir}"
modisco_dir="${base_dir}/01-models/modisco"
bigwigs_dir="${base_dir}/01-models/bigwigs"
bigwigs_signal_dir="${base_dir}/01-models/bigwigs_signal"

modisco_scratch="${modisco_dir}"
preds_dir="${base_dir}/01-models/predictions"
preds_scratch="${preds_dir}"
chrombpnet_models_keep="${qc_dir}/chrombpnet_models_keep.tsv"
chrombpnet_models_keep2="${qc_dir}/chrombpnet_models_keep2.tsv"
# Tool and training defaults
chrombpnet_bin="${CHROMBPNET_BIN:-chrombpnet}"
chrombpnet_train_num_folds="${CHROMBPNET_TRAIN_NUM_FOLDS:-5}"
chrombpnet_data_type="${DATA_TYPE:-ATAC}"
chrombpnet_negative_sampling_ratio="${CHROMBPNET_NEGATIVE_SAMPLING_RATIO:-0.1}"
chrombpnet_bias_cluster="${BIAS_CLUSTER:-1-col_aspn_ogna}"
chrombpnet_bias_min_thresh="${BIAS_MIN_THRESH:-0.4}"
chrombpnet_bias_fold="${BIAS_FOLD:-fold_0}"
chrombpnet_bias_params="${BIAS_PARAMS:-${chrombpnet_bias_cluster}_thresh${chrombpnet_bias_min_thresh}}"
# Temp space for local Mac runs. Keep large intermediates off the nearly-full system volume.
labcluster_scratch="${base_dir}/tmp"

mkdir -p \
  "${data_root}/inputs" \
  "${data_root}/refs" \
  "${base_dir}" \
  "${sample_frags_dir}/fragments" \
  "${sample_frags_dir}/pseudorep1" \
  "${sample_frags_dir}/pseudorep2" \
  "${sample_frags_dir}/pseudorepT" \
  "${cluster_frags_dir}/fragments" \
  "${cluster_frags_dir}/pseudorep1" \
  "${cluster_frags_dir}/pseudorep2" \
  "${cluster_frags_dir}/pseudorepT" \
  "${peaks_dir}" \
  "${chrombpnet_peaks_dir}" \
  "${bigwigs_dir}" \
  "${bigwigs_signal_dir}" \
  "${split_dir}" \
  "${negatives_dir}" \
  "${negatives_subset_dir}" \
  "${anno_peaks_dir}" \
  "${models_dir}" \
  "${bias_dir}" \
  "${qc_dir}" \
  "${contribs_dir}" \
  "${modisco_dir}" \
  "${preds_dir}" \

