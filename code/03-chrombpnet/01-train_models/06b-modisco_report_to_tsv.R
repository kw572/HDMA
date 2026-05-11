# Purpose: this script loops through all modisco reports and converts the HTMLs
# to TSV for easy parsing and re-use later.

library(here)
script_path <- normalizePath(commandArgs(trailingOnly = FALSE), winslash = "/", mustWork = FALSE)
script_file <- sub("^--file=", "", script_path[grepl("^--file=", script_path)])
if (length(script_file) == 0) {
  script_file <- file.path(getwd(), "code/03-chrombpnet/01-train_models/06b-modisco_report_to_tsv.R")
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
work_dir    <- base_dir
out         <- file.path(work_dir, "01-models/modisco_tsv")
modisco_dir <- file.path(work_dir, "01-models/modisco", glue("bias_{bias_params}"))

dir.create(out, showWarnings = FALSE, recursive = TRUE)

#' Convert the motifs.html output from modisco to a TSV, adding the name of the
#' model ("Cluster") as well as the model head used for interpretation ("Model_head")
#' 
#' @param key character, Name of model/cluster
#' @param key characger, Either "counts" or "profile"
modisco_report_to_tsv <- function(key, model_head = "counts") {
  
  message("@ ", key)
  
  modisco_report <- rvest::read_html(glue(
    "{modisco_dir}/{key}/{model_head}_modisco_report/motifs.html"))
  
  # read in MEME format (only for pos patterns though)
  meme <- universalmotif::read_meme(glue(
    "{modisco_dir}/{key}/{key}_memedb.{model_head}.txt"))
  
  # get consensus sequences using universal motif
  meme_df <- data.frame(pattern = map_chr(meme, ~ glue('{key}__pos_patterns.{.x["name"]}')),
                        consensus_fwd = map_chr(meme, ~ .x["consensus"]),
                        consensus_rev = map_chr(meme, ~ universalmotif::motif_rc(.x)["consensus"]))
  
  
  col_names <- modisco_report %>% html_elements("thead") %>% html_table() %>% getElement(1) %>% colnames()
  modisco_df <- modisco_report %>% html_elements("tbody") %>% html_table() %>% 
    getElement(1) %>% 
    set_colnames(col_names) %>% 
    mutate(pattern = paste0(key, "__", pattern)) %>% 
    tibble::add_column("Cluster" = key, "Model_head" = model_head, .before = 1) %>% 
    left_join(meme_df, by = "pattern") %>% 
    dplyr::relocate(consensus_fwd, .after = 3) %>% 
    dplyr::relocate(consensus_rev, .after = consensus_fwd)
  
  data.table::fwrite(modisco_df, file = glue("{out}/{key}_modisco_report.counts.tsv"), sep = "\t")
  
  return(modisco_df)
  
}

models_keep_path <- file.path(work_dir, "01-models/qc/chrombpnet_models_keep.tsv")
if (file.exists(models_keep_path)) {
  chrombpnet_models_keep <- read_tsv(models_keep_path,
                                     col_names = c("Cluster", "Folds_keep", "Cluster_ID"))
  datasets <- unique(chrombpnet_models_keep$Cluster)
} else {
  datasets <- fs::dir_ls(modisco_dir, type = "directory", recurse = FALSE) %>%
    basename()
}

counts_modisco_reports <- map_dfr(datasets, ~ modisco_report_to_tsv(.x))

write_tsv(counts_modisco_reports, glue("{out}/all_modisco_report.motifs.counts.tsv"))

message("@ done.")
