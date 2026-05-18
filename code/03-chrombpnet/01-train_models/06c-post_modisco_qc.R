# Purpose: after QC-ing modisco reports, we exclude models where we
# seemed to be strongly affected by GC bias. These generally came from low-coverage
# cell types with <1M fragments. Therefore, we remove those from our models-to-keep file.

library(here)

script_path <- normalizePath(
  commandArgs(trailingOnly = FALSE),
  winslash = "/",
  mustWork = FALSE
)

script_file <- sub(
  "^--file=",
  "",
  script_path[grepl("^--file=", script_path)]
)

if (length(script_file) == 0) {
  script_file <- file.path(
    getwd(),
    "code/03-chrombpnet/01-train_models/06c-post_modisco_qc.R"
  )
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(scales)
library(glue)
library(purrr)
library(stringr)
library(rvest)
library(universalmotif)

bias_params <- Sys.getenv("BIAS_PARAMS")
work_dir <- Sys.getenv("CHROMBPNET_WORK_DIR")

stopifnot(nzchar(bias_params))
stopifnot(nzchar(work_dir))

message("bias_params = ", bias_params)
message("work_dir = ", work_dir)

out <- file.path(work_dir, "01-models/qc")

dir.create(out, showWarnings = FALSE, recursive = TRUE)

chrombpnet_models_keep <- read_tsv(
  file.path(out, "chrombpnet_models_keep.tsv"),
  col_names = c("Cluster", "Folds_keep", "Cluster_ID"),
  show_col_types = FALSE
)

message("Initial models: ",
        length(unique(chrombpnet_models_keep$Cluster)))

message("Post-MoDISco filtering disabled.")
message("Copying keep.tsv -> keep2.tsv unchanged.")

chrombpnet_models_keep %>%
  write_tsv(file.path(out, "chrombpnet_models_keep2.tsv"))

message("@ done.")