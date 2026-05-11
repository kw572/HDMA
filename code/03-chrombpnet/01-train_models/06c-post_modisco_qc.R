# Purpose: after QC-ing modisco reports, we exclude models where we
# seemed to be strongly affected by GC bias. These generally came from low-coverage
# cell types with <1M fragments. Therefore, we remove those from our models-to-keep file.

library(here)
script_path <- normalizePath(commandArgs(trailingOnly = FALSE), winslash = "/", mustWork = FALSE)
script_file <- sub("^--file=", "", script_path[grepl("^--file=", script_path)])
if (length(script_file) == 0) {
  script_file <- file.path(getwd(), "code/03-chrombpnet/01-train_models/06c-post_modisco_qc.R")
}
source(file.path(dirname(normalizePath(script_file[1], winslash = "/", mustWork = FALSE)), "../config.sh"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(scales)
library(glue)
library(purrr)
library(stringr)
library(rvest) # for parsing modisco HTML reports
library(universalmotif) # for working with motifs

resolve_bias_params <- function() {
  bias_params <- Sys.getenv("BIAS_PARAMS", unset = "")
  if (nzchar(bias_params)) {
    return(bias_params)
  }

  config_path <- file.path(dirname(normalizePath(script_file[1], winslash = "/", mustWork = FALSE)), "../config.sh")
  cmd <- sprintf("source %s >/dev/null 2>&1 && printf '%%s' \"$chrombpnet_bias_params\"", shQuote(config_path))
  trimws(system2("bash", c("-lc", cmd), stdout = TRUE))
}

bias_params <- resolve_bias_params()
out         <- file.path(base_dir, "01-models/qc")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

chrombpnet_models_keep <- read_tsv(file.path(out, "chrombpnet_models_keep.tsv"),
                                   col_names = c("Cluster", "Folds_keep", "Cluster_ID"))
length(unique(chrombpnet_models_keep$Cluster))

models_rm_post_modisco <- Sys.getenv("MODELS_RM_POST_MODISCO", unset = "") %>%
  str_split(",") %>%
  .[[1]] %>%
  trimws() %>%
  discard(~ .x == "")

length(setdiff(chrombpnet_models_keep$Cluster, models_rm_post_modisco))

chrombpnet_models_keep %>% 
  filter(!(Cluster %in% models_rm_post_modisco)) %>% 
  write_tsv(file.path(out, "chrombpnet_models_keep2.tsv"))
