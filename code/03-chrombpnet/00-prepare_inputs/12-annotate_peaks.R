# Purpose: annotate ChromBPNet training peaks using ArchR:::.fastAnnoPeaks, which annotates peaks
#   according to the genomic features they overlap (distal, intronic, exonic, promoter),
#   distance to TSS, GC content, etc
#   see https://github.com/GreenleafLab/ArchR/blob/d9e741c980c7c64e5348c97a74d146cc95f8ba76/R/ReproduciblePeakSet.R#L590

# libraries
library(here)
library(ArchR)
library(BSgenome)
library(dplyr)
library(tidyr)
library(data.table)
library(readr)
library(ggplot2)
library(cowplot)
library(GenomeInfoDb)


# 0. GET ARGS ------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
print(args)

peaks_bed     = args[1]
peaks_tsv_out = args[2]
genome_id     = ifelse(length(args) >= 3, args[3], "hg38")


get_annotation_resources <- function(genome_id) {
  if (genome_id %in% c("danRer11", "GRCz11", "drerio")) {
    if (!requireNamespace("BSgenome.Drerio.UCSC.danRer11", quietly = TRUE)) {
      stop("Missing package: BSgenome.Drerio.UCSC.danRer11")
    }
    if (!requireNamespace("TxDb.Drerio.UCSC.danRer11.refGene", quietly = TRUE)) {
      stop("Missing package: TxDb.Drerio.UCSC.danRer11.refGene")
    }
    if (!requireNamespace("org.Dr.eg.db", quietly = TRUE)) {
      stop("Missing package: org.Dr.eg.db")
    }

    list(
      archr_genome = "danRer11",
      bsgenome = BSgenome.Drerio.UCSC.danRer11::BSgenome.Drerio.UCSC.danRer11,
      txdb = TxDb.Drerio.UCSC.danRer11.refGene::TxDb.Drerio.UCSC.danRer11.refGene,
      orgdb = org.Dr.eg.db::org.Dr.eg.db
    )
  } else if (genome_id %in% c("hg38", "GRCh38", "human")) {
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      stop("Missing package: BSgenome.Hsapiens.UCSC.hg38")
    }
    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
      stop("Missing package: TxDb.Hsapiens.UCSC.hg38.knownGene")
    }
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Missing package: org.Hs.eg.db")
    }

    list(
      archr_genome = "hg38",
      bsgenome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
      txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
      orgdb = org.Hs.eg.db::org.Hs.eg.db
    )
  } else {
    stop(paste("Unsupported genome_id:", genome_id))
  }
}

resources <- get_annotation_resources(genome_id)

normalize_peak_seqlevels <- function(peaks, resources) {
  peak_levels <- as.character(seqlevels(peaks))
  genome_levels <- as.character(seqlevels(resources$bsgenome))

  if (length(peak_levels) == 0 || length(genome_levels) == 0) {
    return(peaks)
  }

  has_peak_chr <- any(grepl("^chr", peak_levels))
  has_genome_chr <- any(grepl("^chr", genome_levels))

  if (!has_peak_chr && has_genome_chr) {
    rename_map <- stats::setNames(paste0("chr", peak_levels), peak_levels)
    peaks <- renameSeqlevels(peaks, rename_map)
  } else if (has_peak_chr && !has_genome_chr) {
    rename_map <- stats::setNames(sub("^chr", "", peak_levels), peak_levels)
    peaks <- renameSeqlevels(peaks, rename_map)
  }

  common_levels <- intersect(seqlevels(peaks), genome_levels)
  peaks <- keepSeqlevels(peaks, common_levels, pruning.mode = "coarse")
  peaks
}



# 1. LOAD DATA -----------------------------------------------------------------

message("@ loading data...")

peaks <- rtracklayer::import.bed(peaks_bed, extraCols = c("a" = "character",
                                                          "b" = "character",
                                                          "c" = "character",
                                                          "d" = "character",
                                                          "e" = "character",
                                                          "f" = "character",
                                                          "summit" = "numeric"))
print(head(peaks))
print(length(peaks))

peaks <- normalize_peak_seqlevels(peaks, resources)
print(head(peaks))
print(length(peaks))



# 2. COMPUTE GENOMIC ANNO ------------------------------------------------------

message("@ computing genomic annotation...")
# First we want to get an annotation per hit that says what type of genomic
# location it occurs in (promoter, intronic, exonic, distal) as per the ArchR definitions.

# set default genome
addArchRGenome(resources$archr_genome)

# get gene annotation

# Avoid - this uses a very old v24 Gencode annotation packaged with ArchR
# gene_annotation <- getArchRGenome(geneAnnotation = TRUE, genomeAnnotation = FALSE) %>% as.list

# Get the latest anno
gene_annotation <- ArchR::createGeneAnnotation(OrgDb = resources$orgdb, TxDb = resources$txdb)

# > gene_annotation$TSS@metadata %>% unlist() %>% tibble::enframe()
# # A tibble: 18 × 2
# name                                                value                                       
# <chr>                                               <chr>                                       
#   1 genomeInfo.Db type                                  TxDb                                        
# 2 genomeInfo.Supporting package                       GenomicFeatures                             
# 3 genomeInfo.Data source                              UCSC                                        
# 4 genomeInfo.Genome                                   hg38                                        
# 5 genomeInfo.Organism                                 Homo sapiens                                
# 6 genomeInfo.Taxonomy ID                              9606                                        
# 7 genomeInfo.UCSC Table                               knownGene                                   
# 8 genomeInfo.UCSC Track                               GENCODE V38                                 
# 9 genomeInfo.Resource URL                             http://genome.ucsc.edu/                     
#   10 genomeInfo.Type of Gene ID                          Entrez Gene ID                              
# 11 genomeInfo.Full dataset                             yes                                         
# 12 genomeInfo.miRBase build ID                         NA                                          
# 13 genomeInfo.Nb of transcripts                        258145                                      
# 14 genomeInfo.Db created by                            GenomicFeatures package from Bioconductor   
# 15 genomeInfo.Creation time                            2021-10-19 10:58:00 -0700 (Tue, 19 Oct 2021)
# 16 genomeInfo.GenomicFeatures version at creation time 1.45.2                                      
# 17 genomeInfo.RSQLite version at creation time         2.2.8                                       
# 18 genomeInfo.DBSCHEMAVERSION                          1.2    


gene_annotation <- gene_annotation %>% as.list()

# annotate peaks
# https://github.com/GreenleafLab/ArchR/blob/d9e741c980c7c64e5348c97a74d146cc95f8ba76/R/ReproduciblePeakSet.R#L590

# about 10 min
Sys.time()
peaks_anno <- ArchR:::.fastAnnoPeaks(peaks, BSgenome = resources$bsgenome, geneAnnotation = gene_annotation, promoterRegion = c(2000, 2000))
Sys.time()
print(length(peaks_anno))



# 3. SAVE TSV OUT --------------------------------------------------------------

message("@ writing out...")

# convert from 1-based coords (GRanges) to 0-based coordinates (BED)
peaks_anno_df <- data.frame(peaks_anno)
peaks_anno_df <- peaks_anno_df %>%
  dplyr::select(chr = seqnames, start, end, width, strand, summit, distToGeneStart,
  nearestGene, peakType, distToTSS, nearestTSS, GC, N)
peaks_anno_df$start <- peaks_anno_df$start - 1
print(head(peaks_anno_df))

# write out to TSV
write_tsv(peaks_anno_df, file = peaks_tsv_out)

message("@ done.")
