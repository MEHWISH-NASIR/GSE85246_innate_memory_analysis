# ============================================================
# 12d — GENOME-WIDE DAY-6 PROMOTER QUANTIFICATION
#
# GSE85245
#
# Purpose:
#   Build one hg19 promoter (TSS +/- 2 kb) per gene and
#   quantify mean BigWig signal for all 18 Day-6 tracks.
#
# Marks:
#   H3K27ac
#   H3K4me1
#   H3K4me3
#
# Conditions:
#   RPMI rep1 / rep2
#   LPS  rep1 / rep2
#   BG   rep1 / rep2
#
# IMPORTANT:
#   - NO scaling or re-normalisation is performed here.
#   - The three GEO NotNormalized RPMI_rep2 tracks remain
#     exactly as deposited.
#   - Scaling will be handled in a separate subsequent step.
# ============================================================

rm(list = ls())

required_packages <- c(
  "rtracklayer",
  "GenomicRanges",
  "GenomicFeatures",
  "GenomeInfoDb",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste(
      "Missing required package(s):",
      paste(missing_packages, collapse = ", ")
    )
  )
}

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(AnnotationDbi)

outdir <- "results/12_kinase_reevaluation"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n========================================\n")
cat("12d — GENOME-WIDE PROMOTER QUANTIFICATION\n")
cat("========================================\n\n")


# ============================================================
# 1. BIGWIG MANIFEST
# ============================================================

files <- list.files(
  "data/chip/day6_formal",
  pattern = "\\.bw$",
  recursive = TRUE,
  full.names = TRUE
)

files <- sort(files)

if (length(files) != 18) {
  stop(
    paste(
      "Expected 18 Day-6 BigWigs but found",
      length(files)
    )
  )
}

sample_manifest <- data.frame(
  file = files,
  mark = basename(dirname(files)),
  sample = basename(files),
  stringsAsFactors = FALSE
)

sample_manifest$condition <- ifelse(
  grepl("^RPMI_", sample_manifest$sample),
  "RPMI",
  ifelse(
    grepl("^LPS_", sample_manifest$sample),
    "LPS",
    ifelse(
      grepl("^BG_", sample_manifest$sample),
      "BG",
      NA_character_
    )
  )
)

sample_manifest$replicate <- ifelse(
  grepl("_rep1_", sample_manifest$sample),
  "rep1",
  ifelse(
    grepl("_rep2_", sample_manifest$sample),
    "rep2",
    NA_character_
  )
)

sample_manifest$normalization_label <- ifelse(
  grepl(
    "NotNormalized",
    sample_manifest$sample
  ),
  "NotNormalized",
  "Normalized"
)

sample_manifest$sample_id <- paste(
  sample_manifest$mark,
  sample_manifest$condition,
  sample_manifest$replicate,
  sep = "_"
)

if (anyNA(sample_manifest$condition) ||
    anyNA(sample_manifest$replicate)) {
  stop("Could not parse condition/replicate from one or more filenames.")
}

if (anyDuplicated(sample_manifest$sample_id)) {
  stop("Duplicate sample IDs detected.")
}

write.csv(
  sample_manifest,
  file.path(
    outdir,
    "12d_day6_sample_manifest.csv"
  ),
  row.names = FALSE
)

cat("BigWigs:", nrow(sample_manifest), "\n")

cat("\nTracks by mark / condition:\n")
print(
  with(
    sample_manifest,
    table(mark, condition)
  )
)

cat("\nNormalization labels:\n")
print(
  with(
    sample_manifest,
    table(mark, condition, normalization_label)
  )
)


# ============================================================
# 2. BUILD hg19 GENE-LEVEL PROMOTERS
# ============================================================

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

gene_gr <- suppressMessages(
  GenomicFeatures::genes(
    txdb,
    single.strand.genes.only = TRUE
  )
)

standard_chr <- c(
  paste0("chr", 1:22),
  "chrX",
  "chrY"
)

gene_gr <- gene_gr[
  as.character(seqnames(gene_gr)) %in%
    standard_chr
]

cat(
  "\nGene ranges on standard chromosomes:",
  length(gene_gr),
  "\n"
)


# ============================================================
# 3. ANNOTATE ENTREZ -> SYMBOL / ENSEMBL
# ============================================================

entrez <- names(gene_gr)

symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = entrez,
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

ensembl <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = entrez,
  keytype = "ENTREZID",
  column = "ENSEMBL",
  multiVals = "first"
)

mcols(gene_gr)$ENTREZID <- entrez
mcols(gene_gr)$SYMBOL <-
  unname(symbol[entrez])
mcols(gene_gr)$ENSEMBL <-
  unname(ensembl[entrez])


# ============================================================
# 4. PROMOTERS = TSS +/- 2 kb
# ============================================================

promoter_gr <- GenomicRanges::promoters(
  gene_gr,
  upstream = 2000,
  downstream = 2000
)

promoter_gr <- GenomicRanges::trim(
  promoter_gr
)

mcols(promoter_gr)$ENTREZID <-
  mcols(gene_gr)$ENTREZID

mcols(promoter_gr)$SYMBOL <-
  mcols(gene_gr)$SYMBOL

mcols(promoter_gr)$ENSEMBL <-
  mcols(gene_gr)$ENSEMBL


# ============================================================
# 5. REMOVE INVALID / DUPLICATED GENE RECORDS
# ============================================================

valid <- (
  width(promoter_gr) > 0 &
  !is.na(mcols(promoter_gr)$ENTREZID) &
  mcols(promoter_gr)$ENTREZID != ""
)

promoter_gr <- promoter_gr[valid]

if (anyDuplicated(mcols(promoter_gr)$ENTREZID)) {
  stop(
    "Unexpected duplicated ENTREZIDs after gene-level promoter construction."
  )
}

cat(
  "Final gene-level promoters:",
  length(promoter_gr),
  "\n"
)


# ============================================================
# 6. PROMOTER MANIFEST
# ============================================================

promoter_manifest <- data.frame(
  ENTREZID =
    mcols(promoter_gr)$ENTREZID,
  SYMBOL =
    mcols(promoter_gr)$SYMBOL,
  ENSEMBL =
    mcols(promoter_gr)$ENSEMBL,
  chr =
    as.character(seqnames(promoter_gr)),
  start =
    start(promoter_gr),
  end =
    end(promoter_gr),
  strand =
    as.character(strand(promoter_gr)),
  width =
    width(promoter_gr),
  stringsAsFactors = FALSE
)

write.csv(
  promoter_manifest,
  file.path(
    outdir,
    "12d_hg19_TSS_2kb_promoters.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 7. CHECK TARGET GENES
# ============================================================

targets <- c(
  "JAK3",
  "EPHB2",
  "MET",
  "MAP3K8",
  "BMPR1A"
)

cat("\n========================================\n")
cat("TARGET PROMOTER CHECK\n")
cat("========================================\n\n")

target_manifest <- promoter_manifest[
  promoter_manifest$SYMBOL %in% targets,
]

print(
  target_manifest,
  row.names = FALSE
)

missing_targets <- setdiff(
  targets,
  target_manifest$SYMBOL
)

if (length(missing_targets) > 0) {
  stop(
    paste(
      "Missing target promoter(s):",
      paste(missing_targets, collapse = ", ")
    )
  )
}


# ============================================================
# 8. PROMOTER QUANTIFICATION FUNCTION
#
# Query one chromosome at a time.
# This avoids the multi-chromosome summary issue encountered
# during the earlier global BigWig QC.
# ============================================================

quantify_one_bigwig <- function(
    bw_file,
    promoters
) {

  bw <- rtracklayer::BigWigFile(
    bw_file
  )

  bw_seq <- GenomeInfoDb::seqnames(
    GenomeInfoDb::seqinfo(bw)
  )

  output <- rep(
    NA_real_,
    length(promoters)
  )

  chromosomes <- unique(
    as.character(
      seqnames(promoters)
    )
  )

  for (chr in chromosomes) {

    idx <- which(
      as.character(
        seqnames(promoters)
      ) == chr
    )

    if (!(chr %in% bw_seq)) {
      warning(
        paste(
          "Chromosome",
          chr,
          "absent from",
          bw_file
        )
      )

      output[idx] <- NA_real_
      next
    }

    ranges_chr <- promoters[idx]

    x <- rtracklayer::summary(
      bw,
      ranges_chr,
      size = 1L,
      type = "mean",
      defaultValue = 0,
      as = "matrix"
    )

    x <- as.numeric(x)

    if (length(x) != length(idx)) {
      stop(
        paste(
          "BigWig summary length mismatch:",
          basename(bw_file),
          chr,
          "expected",
          length(idx),
          "but obtained",
          length(x)
        )
      )
    }

    output[idx] <- x
  }

  output
}


# ============================================================
# 9. QUANTIFY ALL 18 TRACKS
# ============================================================

signal_matrix <- matrix(
  NA_real_,
  nrow = length(promoter_gr),
  ncol = nrow(sample_manifest)
)

rownames(signal_matrix) <-
  promoter_manifest$ENTREZID

colnames(signal_matrix) <-
  sample_manifest$sample_id

cat("\n========================================\n")
cat("QUANTIFYING 18 BIGWIG TRACKS\n")
cat("========================================\n\n")

for (i in seq_len(nrow(sample_manifest))) {

  cat(
    sprintf(
      "[%02d/%02d] %s\n",
      i,
      nrow(sample_manifest),
      sample_manifest$sample_id[i]
    )
  )

  signal_matrix[, i] <-
    quantify_one_bigwig(
      sample_manifest$file[i],
      promoter_gr
    )

  cat(
    "    finite promoters:",
    sum(
      is.finite(
        signal_matrix[, i]
      )
    ),
    "/",
    nrow(signal_matrix),
    "\n"
  )
}


# ============================================================
# 10. MATRIX QC
# ============================================================

if (anyNA(signal_matrix)) {

  na_count <- sum(
    is.na(signal_matrix)
  )

  warning(
    paste(
      "Promoter matrix contains",
      na_count,
      "NA values."
    )
  )
}

cat("\nSignal range:\n")
print(
  range(
    signal_matrix,
    na.rm = TRUE
  )
)

cat("\nColumn medians:\n")
print(
  apply(
    signal_matrix,
    2,
    median,
    na.rm = TRUE
  )
)


# ============================================================
# 11. SAVE RAW / AS-DEPOSITED PROMOTER MATRIX
# ============================================================

saveRDS(
  list(
    promoters = promoter_gr,
    promoter_manifest =
      promoter_manifest,
    sample_manifest =
      sample_manifest,
    signal_matrix_raw =
      signal_matrix
  ),
  file.path(
    outdir,
    "12d_day6_promoter_signal_raw.rds"
  )
)

signal_table <- cbind(
  promoter_manifest,
  as.data.frame(
    signal_matrix,
    check.names = FALSE
  )
)

write.csv(
  signal_table,
  file.path(
    outdir,
    "12d_day6_promoter_signal_raw.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. TARGET-GENE RAW SIGNAL AUDIT
# ============================================================

target_idx <- which(
  promoter_manifest$SYMBOL %in%
    targets
)

target_signal <- cbind(
  promoter_manifest[target_idx, ],
  as.data.frame(
    signal_matrix[target_idx, , drop = FALSE],
    check.names = FALSE
  )
)

target_signal <- target_signal[
  match(
    targets,
    target_signal$SYMBOL
  ),
]

write.csv(
  target_signal,
  file.path(
    outdir,
    "12d_target_gene_raw_promoter_signal.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("TARGET RAW PROMOTER SIGNALS\n")
cat("========================================\n\n")

print(
  target_signal,
  row.names = FALSE
)


# ============================================================
# 13. COMPLETE
# ============================================================

cat("\nSaved:\n")

cat(
  " ",
  file.path(
    outdir,
    "12d_hg19_TSS_2kb_promoters.csv"
  ),
  "\n"
)

cat(
  " ",
  file.path(
    outdir,
    "12d_day6_sample_manifest.csv"
  ),
  "\n"
)

cat(
  " ",
  file.path(
    outdir,
    "12d_day6_promoter_signal_raw.rds"
  ),
  "\n"
)

cat(
  " ",
  file.path(
    outdir,
    "12d_day6_promoter_signal_raw.csv"
  ),
  "\n"
)

cat(
  " ",
  file.path(
    outdir,
    "12d_target_gene_raw_promoter_signal.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("12d COMPLETE\n")
cat("========================================\n")
