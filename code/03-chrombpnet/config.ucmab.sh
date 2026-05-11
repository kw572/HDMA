#!/usr/bin/env bash
set -euo pipefail

chrombpnet_root="/scratch1/kuangtse/HDMA/code/03-chrombpnet"
data_root="${chrombpnet_root}/data/local_run_chr"

rds_file="${data_root}/inputs/ucmab_14dpf.rds"
fragment_file="${data_root}/inputs/fragments.chr.tsv.gz"

cluster_col="cell_type"
sample_col="orig.ident"
sample_name=""
assay_name="peaks"

refs="${data_root}/refs"
ref_fasta="${refs}/Danio_rerio.GRCz11.chr.dna.toplevel.fa"
chromsizes="${refs}/zf.chr.major.chrom.sizes"
blacklist="${refs}/blacklist.GRCz11.chr.bed"
genome_size=""

base_dir="${data_root}/work"
sample_frags_dir="${base_dir}/00-inputs/sample_fragments"
cluster_frags_dir="${base_dir}/00-inputs/cluster_fragments"
peaks_dir="${base_dir}/00-inputs/peaks"
chrombpnet_peaks_dir="${base_dir}/00-inputs/chrombpnet_peaks"
split_dir="${base_dir}/00-inputs/splits"
negatives_dir="${base_dir}/00-inputs/negatives"
negatives_subset_dir="${base_dir}/00-inputs/negatives_subset"
anno_peaks_dir="${base_dir}/00-inputs/annotated_peaks"

models_dir="${base_dir}/01-models"
bias_dir="${base_dir}/01-models/bias"
qc_dir="${base_dir}/01-models/qc"
contribs_dir="${base_dir}/01-models/contribs"
contribs_scratch="${contribs_dir}"
modisco_dir="${base_dir}/01-models/modisco"
modisco_scratch="${modisco_dir}"
preds_dir="${base_dir}/01-models/predictions"
preds_scratch="${preds_dir}"
chrombpnet_models_keep="${qc_dir}/chrombpnet_models_keep.tsv"
chrombpnet_models_keep2="${qc_dir}/chrombpnet_models_keep2.tsv"

labcluster_scratch="${TMPDIR:-/tmp}"

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
  "${split_dir}" \
  "${negatives_dir}" \
  "${negatives_subset_dir}" \
  "${anno_peaks_dir}" \
  "${models_dir}" \
  "${bias_dir}" \
  "${qc_dir}" \
  "${contribs_dir}" \
  "${modisco_dir}" \
  "${preds_dir}"
