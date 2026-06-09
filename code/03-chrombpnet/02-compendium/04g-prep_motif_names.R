# Purpose: extract out a set of compiled motif names for non-composites and composites,
# in preparation for running hit calling on the two sets. Motif names
# need to be in the form "pos_patterns.pattern_3", and these are the MoDISco
# motif names, not the custom labels.

library(dplyr)
library(tidyr)
library(readr)
library(glue)

out <- Sys.getenv("MODISCO_COMP_ANNO_DIR", unset = "")
if (!nzchar(out)) {
  out <- here::here("output/03-chrombpnet/02-compendium/modisco_compiled_anno/")
}
annotation_tsv <- Sys.getenv("MOTIF_ANNOTATION_FILE", unset = "04d-ChromBPNet_de_novo_motifs.tsv")

motifs_compiled <- read_tsv(annotation_tsv)

table(motifs_compiled$category)
# 
# base  base_with_flank composite_hetero   composite_homo          exclude          partial           repeat       unresolved 
# 147              128              261              137              107                4                7               43

motifs_compiled %>%
  filter(!grepl("composite", category)) %>% 
  mutate(pattern = case_when(
      pattern_class == "pos" ~ paste0("pos_patterns.", pattern),
      pattern_class == "neg" ~ paste0("neg_patterns.", pattern))) %>%
  pull(pattern) %>% 
  write_lines(glue("{out}/motif_names.non_composites.tsv"))

motifs_compiled %>%
  filter(grepl("composite", category)) %>% 
  mutate(pattern = case_when(
      pattern_class == "pos" ~ paste0("pos_patterns.", pattern),
      pattern_class == "neg" ~ paste0("neg_patterns.", pattern))) %>%
  pull(pattern) %>% 
  write_lines(glue("{out}/motif_names.composites.tsv"))
