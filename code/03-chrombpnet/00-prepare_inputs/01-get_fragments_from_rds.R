#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
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

barcode_variants <- function(barcode) {
  values <- c(
    barcode,
    gsub("#", "_", barcode, fixed = TRUE),
    gsub("_", "#", barcode, fixed = TRUE)
  )

  if (grepl("#", barcode, fixed = TRUE)) {
    values <- c(values, sub("^.*#", "", barcode))
  }
  if (grepl("_", barcode, fixed = TRUE)) {
    values <- c(values, sub("^.*_", "", barcode))
  }

  unique(values[nzchar(values)])
}

build_named_map <- function(keys, values) {
  stats::setNames(as.list(as.character(values)), keys)
}

collect_unique_bare_matches <- function(cells, indices) {
  bare <- sub("^[^#_]+[#_]", "", cells)
  bare_matches <- split(indices, bare)
  bare_matches[vapply(bare_matches, length, integer(1)) == 1L]
}

resolve_fragment_path <- function(obj, assay_name, explicit_path) {
  if (!is.null(explicit_path) && nzchar(explicit_path)) {
    return(normalizePath(explicit_path, mustWork = TRUE))
  }

  if (!assay_name %in% names(obj@assays)) {
    stop("Assay not found in Seurat object: ", assay_name)
  }

  assay_obj <- obj[[assay_name]]
  fragment_objects <- tryCatch(Signac::Fragments(assay_obj), error = function(e) NULL)
  if (is.null(fragment_objects) || length(fragment_objects) == 0L) {
    stop(
      "No fragment object found in assay '", assay_name,
      "'. Pass --fragments explicitly."
    )
  }

  fragment_path <- tryCatch({
    slot(fragment_objects[[1]], "path")
  }, error = function(e) NULL)

  if (is.null(fragment_path) || !nzchar(fragment_path)) {
    stop(
      "Unable to resolve fragment path from assay '", assay_name,
      "'. Pass --fragments explicitly."
    )
  }

  normalizePath(fragment_path, mustWork = TRUE)
}

load_reference_seqlevels <- function(chromsizes_path) {
  if (is.null(chromsizes_path) || !nzchar(chromsizes_path)) {
    return(character(0))
  }

  chromsizes_path <- normalizePath(chromsizes_path, mustWork = TRUE)
  chrom_df <- read.table(
    chromsizes_path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
  )
  unique(as.character(chrom_df[[1]]))
}

normalize_chrom_name <- function(chrom, reference_seqlevels) {
  if (length(reference_seqlevels) == 0L) {
    return(chrom)
  }

  if (chrom %in% reference_seqlevels) {
    return(chrom)
  }

  alt <- if (startsWith(chrom, "chr")) {
    sub("^chr", "", chrom)
  } else {
    paste0("chr", chrom)
  }

  if (alt %in% reference_seqlevels) {
    return(alt)
  }

  NA_character_
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

required_args <- c("workdir", "outdir", "rds", "cluster-col")
missing_args <- required_args[!required_args %in% names(args)]
if (length(missing_args) > 0L) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

workdir <- normalizePath(args[["workdir"]], mustWork = FALSE)
outdir <- normalizePath(args[["outdir"]], mustWork = FALSE)
rds_path <- normalizePath(args[["rds"]], mustWork = TRUE)
cluster_col <- args[["cluster-col"]]
sample_col <- args[["sample-col"]] %||% NULL
sample_name <- args[["sample-name"]] %||% NULL
assay_name <- args[["assay"]] %||% "peaks"
chromsizes_path <- args[["chromsizes"]] %||% NULL

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "fragments"), recursive = TRUE, showWarnings = FALSE)

setwd(workdir)

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

manifest <- unique(
  meta[, c("cluster_label", "sample_label", "cluster_token", "sample_token", "output_stem")]
)
manifest <- manifest[order(manifest$output_stem), , drop = FALSE]
manifest_path <- file.path(outdir, "fragments_manifest.tsv")
write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("@ wrote manifest: ", manifest_path)

fragment_path <- resolve_fragment_path(obj, assay_name, args[["fragments"]] %||% NULL)
message("@ using fragment file: ", fragment_path)
reference_seqlevels <- load_reference_seqlevels(chromsizes_path)
if (length(reference_seqlevels) > 0L) {
  message("@ loaded ", length(reference_seqlevels), " reference seqlevels from: ", normalizePath(chromsizes_path, mustWork = TRUE))
}

cell_indices <- seq_len(nrow(meta))
exact_map <- build_named_map(meta$cell_name, cell_indices)
hash_map <- build_named_map(gsub("#", "_", meta$cell_name, fixed = TRUE), cell_indices)
underscore_map <- build_named_map(gsub("_", "#", meta$cell_name, fixed = TRUE), cell_indices)
bare_map <- collect_unique_bare_matches(meta$cell_name, cell_indices)

resolve_index <- function(barcode) {
  candidates <- integer(0)

  for (lookup in list(exact_map, hash_map, underscore_map)) {
    hit <- lookup[[barcode]]
    if (!is.null(hit)) {
      candidates <- c(candidates, as.integer(hit))
    }
  }

  bare_hit <- bare_map[[barcode]]
  if (!is.null(bare_hit)) {
    candidates <- c(candidates, as.integer(bare_hit))
  }

  candidates <- unique(candidates)
  if (length(candidates) == 1L) {
    candidates[[1]]
  } else {
    NA_integer_
  }
}

open_input <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

writers <- new.env(parent = emptyenv())
get_writer <- function(stem) {
  if (exists(stem, envir = writers, inherits = FALSE)) {
    return(get(stem, envir = writers, inherits = FALSE))
  }
  path <- file.path(outdir, "fragments", paste0(stem, ".tsv"))
  con <- file(path, open = "wt")
  assign(stem, con, envir = writers)
  con
}

on.exit({
  open_names <- ls(envir = writers, all.names = TRUE)
  for (name in open_names) {
    close(get(name, envir = writers, inherits = FALSE))
  }
}, add = TRUE)

counts <- c(total = 0L, written = 0L, unmatched = 0L, malformed = 0L, dropped_chrom = 0L, normalized_chrom = 0L)
input_con <- open_input(fragment_path)
on.exit(close(input_con), add = TRUE)

repeat {
  line <- readLines(input_con, n = 1L)
  if (length(line) == 0L) {
    break
  }

  counts[["total"]] <- counts[["total"]] + 1L
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(fields) < 4L) {
    counts[["malformed"]] <- counts[["malformed"]] + 1L
    next
  }

  normalized_chrom <- normalize_chrom_name(fields[[1]], reference_seqlevels)
  if (is.na(normalized_chrom)) {
    counts[["dropped_chrom"]] <- counts[["dropped_chrom"]] + 1L
    next
  }
  if (!identical(normalized_chrom, fields[[1]])) {
    counts[["normalized_chrom"]] <- counts[["normalized_chrom"]] + 1L
  }
  fields[[1]] <- normalized_chrom

  barcode <- fields[[4]]
  idx <- resolve_index(barcode)
  if (is.na(idx)) {
    counts[["unmatched"]] <- counts[["unmatched"]] + 1L
    next
  }

  stem <- meta$output_stem[[idx]]
  writer <- get_writer(stem)
  writeLines(paste(fields[1:4], collapse = "\t"), writer)
  counts[["written"]] <- counts[["written"]] + 1L
}

stats_path <- file.path(outdir, "fragments_stats.tsv")
stats_df <- data.frame(
  metric = names(counts),
  value = as.integer(counts),
  stringsAsFactors = FALSE
)
write.table(stats_df, stats_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("@ finished fragment extraction")
message("@ total lines: ", counts[["total"]])
message("@ written lines: ", counts[["written"]])
message("@ unmatched lines: ", counts[["unmatched"]])
message("@ malformed lines: ", counts[["malformed"]])
