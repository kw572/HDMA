# Purpose: get the models metrics and specify which are kept. Prepare the sup
# table for the manuscript.

library(here)
script_path <- normalizePath(commandArgs(trailingOnly = FALSE), winslash = "/", mustWork = FALSE)
script_file <- sub("^--file=", "", script_path[grepl("^--file=", script_path)])
if (length(script_file) == 0) {
  script_file <- file.path(getwd(), "code/03-chrombpnet/01-train_models/09-prep_chrombpnet_QC_table.R")
}
source(file.path(dirname(normalizePath(script_file[1], winslash = "/", mustWork = FALSE)), "../config.sh"))
library(dplyr)
library(tidyr)
library(readr)
library(glue)
library(purrr)

out         <- file.path(base_dir, "01-models/qc")

dir.create(out, showWarnings = FALSE, recursive = TRUE)

# load outs
all_metrics <- read_tsv(file.path(out, "chrombpnet_metrics.tsv"))
models_keep1 <- read_tsv(file.path(out, "chrombpnet_models_keep.tsv"), col_names = c("Cluster_ChromBPNet", "Folds", "Cluster"))
models_keep2 <- read_tsv(file.path(out, "chrombpnet_models_keep2.tsv"), col_names = c("Cluster_ChromBPNet", "Folds", "Cluster"))

all_metrics <- all_metrics %>% 
  mutate(Pass_QC = case_when(
    Cluster_ChromBPNet %in% models_keep2$Cluster_ChromBPNet ~ TRUE,
    TRUE ~ FALSE
  ), Reason = case_when(
    Pass_QC ~ "N/A",
    !(Cluster_ChromBPNet %in% models_keep1$Cluster_ChromBPNet) ~ "Failed Spearman correlation threshold",
    TRUE ~ "Failed post-MoDISco QC"
  )) %>%
  dplyr::select(Organ = organ, Cluster_ChromBPNet, L1_annot, Total_n_fragments = total_frags, Model, Fold, Pass_QC, Reason, matches("metrics"), matches("tn5"))

all_metrics %>%
  write_tsv(glue("{out}/TABLE_chrombpnet_qc_metrics.tsv"))
