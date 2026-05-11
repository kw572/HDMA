library(AnnotationDbi)
library(BSgenome)
library(Biostrings)
library(GenomicFeatures)
library(GenomeInfoDb)
library(GenomicRanges)
library(readr)
library(rtracklayer)


args <- commandArgs(trailingOnly = TRUE)
print(args)

peaks_bed <- args[1]
peaks_tsv_out <- args[2]
genome_id <- ifelse(length(args) >= 3, args[3], "hg38")


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
      bsgenome = BSgenome.Drerio.UCSC.danRer11::BSgenome.Drerio.UCSC.danRer11,
      txdb = TxDb.Drerio.UCSC.danRer11.refGene::TxDb.Drerio.UCSC.danRer11.refGene,
      orgdb = org.Dr.eg.db::org.Dr.eg.db,
      orgdb_keytype = "ENTREZID"
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
      bsgenome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
      txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
      orgdb = org.Hs.eg.db::org.Hs.eg.db,
      orgdb_keytype = "ENTREZID"
    )
  } else {
    stop(paste("Unsupported genome_id:", genome_id))
  }
}


normalize_peak_seqlevels <- function(peaks, bsgenome) {
  peak_levels <- as.character(seqlevels(peaks))
  genome_levels <- as.character(seqlevels(bsgenome))

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
  keepSeqlevels(peaks, common_levels, pruning.mode = "coarse")
}


make_gene_annotation <- function(txdb, orgdb, orgdb_keytype) {
  genes_gr <- suppressWarnings(GenomicFeatures::genes(txdb))
  exons_gr <- GenomicFeatures::exons(txdb)
  tx_gr <- GenomicFeatures::transcripts(txdb, columns = c("tx_name", "gene_id"))

  gene_ids <- as.character(mcols(genes_gr)$gene_id)
  gene_symbols <- AnnotationDbi::mapIds(
    orgdb,
    keys = gene_ids,
    column = "SYMBOL",
    keytype = orgdb_keytype,
    multiVals = "first"
  )
  mcols(genes_gr)$symbol <- unname(gene_symbols[gene_ids])

  tx_gene_ids <- as.character(mcols(tx_gr)$gene_id)
  tx_symbols <- AnnotationDbi::mapIds(
    orgdb,
    keys = tx_gene_ids,
    column = "SYMBOL",
    keytype = orgdb_keytype,
    multiVals = "first"
  )
  mcols(tx_gr)$symbol <- unname(tx_symbols[tx_gene_ids])

  list(
    genes = genes_gr,
    exons = exons_gr,
    tss = resize(tx_gr, width = 1, fix = "start")
  )
}


annotate_peaks <- function(peaks, bsgenome, gene_annotation, promoter_region = c(2000, 2000)) {
  peak_summits <- resize(peaks, width = 1, fix = "center")

  gene_starts <- resize(gene_annotation$genes, width = 1, fix = "start")
  gene_dist <- distanceToNearest(peak_summits, gene_starts, ignore.strand = TRUE)
  mcols(peaks)$distToGeneStart <- mcols(gene_dist)$distance
  mcols(peaks)$nearestGene <- mcols(gene_annotation$genes)$symbol[subjectHits(gene_dist)]

  promoters_gr <- promoters(gene_annotation$genes, upstream = promoter_region[1], downstream = promoter_region[2])
  overlap_promoters <- overlapsAny(peak_summits, promoters_gr, ignore.strand = TRUE)
  overlap_genes <- overlapsAny(peak_summits, gene_annotation$genes, ignore.strand = TRUE)
  overlap_exons <- overlapsAny(peak_summits, gene_annotation$exons, ignore.strand = TRUE)

  peak_type <- rep("Distal", length(peaks))
  peak_type[overlap_genes & overlap_exons] <- "Exonic"
  peak_type[overlap_genes & !overlap_exons] <- "Intronic"
  peak_type[overlap_promoters] <- "Promoter"
  mcols(peaks)$peakType <- peak_type

  tss_dist <- distanceToNearest(peak_summits, gene_annotation$tss, ignore.strand = TRUE)
  mcols(peaks)$distToTSS <- mcols(tss_dist)$distance
  nearest_tss_symbol <- mcols(gene_annotation$tss)$symbol[subjectHits(tss_dist)]
  nearest_tss_name <- mcols(gene_annotation$tss)$tx_name[subjectHits(tss_dist)]
  mcols(peaks)$nearestTSS <- ifelse(
    is.na(nearest_tss_symbol) | nearest_tss_symbol == "",
    nearest_tss_name,
    nearest_tss_symbol
  )

  nuc_freq <- Biostrings::alphabetFrequency(BSgenome::getSeq(bsgenome, peaks))
  denom_cols <- intersect(colnames(nuc_freq), c("A", "C", "G", "T", "N"))
  denom <- rowSums(nuc_freq[, denom_cols, drop = FALSE])
  gc_num <- rowSums(nuc_freq[, intersect(colnames(nuc_freq), c("G", "C")), drop = FALSE])
  n_num <- nuc_freq[, "N", drop = TRUE]
  denom[denom == 0] <- NA_real_
  mcols(peaks)$GC <- round(gc_num / denom, 4)
  mcols(peaks)$N <- round(n_num / denom, 4)

  peaks
}


message("@ loading resources...")
resources <- get_annotation_resources(genome_id)
gene_annotation <- make_gene_annotation(resources$txdb, resources$orgdb, resources$orgdb_keytype)

message("@ loading peaks...")
peaks <- rtracklayer::import.bed(
  peaks_bed,
  extraCols = c(
    "a" = "character",
    "b" = "character",
    "c" = "character",
    "d" = "character",
    "e" = "character",
    "f" = "character",
    "summit" = "numeric"
  )
)
peaks <- normalize_peak_seqlevels(peaks, resources$bsgenome)

message("@ annotating peaks...")
peaks_anno <- annotate_peaks(peaks, resources$bsgenome, gene_annotation, promoter_region = c(2000, 2000))

message("@ writing out...")
peaks_anno_df <- data.frame(
  chr = as.character(seqnames(peaks_anno)),
  start = start(peaks_anno) - 1L,
  end = end(peaks_anno),
  width = width(peaks_anno),
  strand = as.character(strand(peaks_anno)),
  summit = mcols(peaks_anno)$summit,
  distToGeneStart = mcols(peaks_anno)$distToGeneStart,
  nearestGene = mcols(peaks_anno)$nearestGene,
  peakType = mcols(peaks_anno)$peakType,
  distToTSS = mcols(peaks_anno)$distToTSS,
  nearestTSS = mcols(peaks_anno)$nearestTSS,
  GC = mcols(peaks_anno)$GC,
  N = mcols(peaks_anno)$N,
  check.names = FALSE
)
readr::write_tsv(peaks_anno_df, file = peaks_tsv_out)

message("@ done.")
