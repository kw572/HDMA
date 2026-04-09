#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

parse_args <- function(args) {
  parsed <- list()
  idx <- 1L
  while (idx <= length(args)) {
    key <- args[[idx]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    if (idx == length(args)) {
      stop("Missing value for argument: ", key)
    }
    parsed[[sub("^--", "", key)]] <- args[[idx + 1L]]
    idx <- idx + 2L
  }
  parsed
}

sanitize_token <- function(values) {
  cleaned <- gsub("[^A-Za-z0-9._-]+", "_", as.character(values))
  cleaned <- gsub("^_+|_+$", "", cleaned)
  cleaned[cleaned == ""] <- "unknown"
  cleaned
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("rds", "mapping-out", "manifest-out", "cluster-col")
missing_args <- required_args[!required_args %in% names(args)]
if (length(missing_args) > 0L) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

rds_path <- normalizePath(args[["rds"]], mustWork = TRUE)
mapping_out <- normalizePath(args[["mapping-out"]], mustWork = FALSE)
manifest_out <- normalizePath(args[["manifest-out"]], mustWork = FALSE)
cluster_col <- args[["cluster-col"]]
sample_col <- args[["sample-col"]] %||% NULL
sample_name <- args[["sample-name"]] %||% NULL

dir.create(dirname(mapping_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(manifest_out), recursive = TRUE, showWarnings = FALSE)

message("@ loading Seurat object: ", rds_path)
obj <- readRDS(rds_path)
meta <- obj@meta.data
meta$cell_name <- rownames(meta)

if (!cluster_col %in% colnames(meta)) {
  stop("Cluster column not found in Seurat metadata: ", cluster_col)
}

if (is.null(sample_col)) {
  if ("orig.ident" %in% colnames(meta)) {
    sample_col <- "orig.ident"
    message("@ using sample column from metadata: ", sample_col)
  } else if (!is.null(sample_name) && nzchar(sample_name)) {
    sample_col <- ".sample_name"
    meta[[sample_col]] <- sample_name
    message("@ using fixed sample name: ", sample_name)
  } else if (!is.null(obj@project.name) && nzchar(obj@project.name)) {
    sample_col <- ".sample_name"
    meta[[sample_col]] <- obj@project.name
    message("@ using Seurat project name as sample: ", obj@project.name)
  } else {
    sample_col <- ".sample_name"
    meta[[sample_col]] <- "sample1"
    message("@ sample column not provided; using fallback sample name: sample1")
  }
}

if (!sample_col %in% colnames(meta)) {
  stop("Sample column not found in Seurat metadata: ", sample_col)
}

meta <- meta[!is.na(meta[[cluster_col]]) & !is.na(meta[[sample_col]]), , drop = FALSE]
if (nrow(meta) == 0L) {
  stop("No cells remain after filtering missing cluster/sample metadata")
}

meta$cluster_label <- as.character(meta[[cluster_col]])
meta$sample_label <- as.character(meta[[sample_col]])
meta$cluster_token <- sanitize_token(meta$cluster_label)
meta$sample_token <- sanitize_token(meta$sample_label)
meta$output_stem <- paste(meta$cluster_token, meta$sample_token, sep = "__")

mapping_df <- meta[, c("cell_name", "output_stem")]
write.table(mapping_df, mapping_out, sep = "\t", quote = FALSE, row.names = FALSE)

manifest_df <- unique(
  meta[, c("cluster_label", "sample_label", "cluster_token", "sample_token", "output_stem")]
)
manifest_df <- manifest_df[order(manifest_df$output_stem), , drop = FALSE]
write.table(manifest_df, manifest_out, sep = "\t", quote = FALSE, row.names = FALSE)

message("@ wrote mapping: ", mapping_out)
message("@ wrote manifest: ", manifest_out)
message("@ exported cells: ", nrow(mapping_df))
