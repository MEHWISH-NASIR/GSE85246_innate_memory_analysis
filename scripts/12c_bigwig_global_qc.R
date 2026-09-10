# ============================================================
# 12c — DAY-6 BIGWIG GLOBAL SIGNAL QC
#
# Purpose:
#   Compare global signal scale across the 18 Day-6 BigWigs
#   before correcting the three NotNormalized RPMI_rep2 tracks.
#
# Uses rtracklayer BigWig summaries.
# No files are modified.
# ============================================================

rm(list = ls())

required <- c(
  "rtracklayer",
  "GenomicRanges",
  "GenomeInfoDb"
)

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  stop(
    paste(
      "Missing packages:",
      paste(missing, collapse = ", ")
    )
  )
}

library(rtracklayer)
library(GenomicRanges)
library(GenomeInfoDb)

indir <- "data/chip/day6_formal"
outdir <- "results/12_kinase_reevaluation"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

files <- list.files(
  indir,
  pattern = "\\.bw$",
  recursive = TRUE,
  full.names = TRUE
)

files <- sort(files)

cat("\n========================================\n")
cat("12c — BIGWIG GLOBAL SIGNAL QC\n")
cat("========================================\n\n")

cat("BigWigs found:", length(files), "\n\n")

if (length(files) != 18) {
  stop(
    paste(
      "Expected 18 BigWigs but found",
      length(files)
    )
  )
}

results <- vector(
  "list",
  length(files)
)

for (i in seq_along(files)) {

  f <- files[i]

  cat(
    sprintf(
      "[%02d/18] %s\n",
      i,
      f
    )
  )

  bw <- rtracklayer::BigWigFile(f)

  si <- GenomeInfoDb::seqinfo(bw)

  if (length(si) == 0) {
    stop(
      paste(
        "Could not read seqinfo:",
        f
      )
    )
  }

  # Use autosomes + chrX + chrY only
  standard_chr <- c(
    paste0("chr", 1:22),
    "chrX",
    "chrY"
  )

  keep_chr <- intersect(
    standard_chr,
    GenomeInfoDb::seqnames(si)
  )

  si2 <- si[keep_chr]

  gr <- as(
    si2,
    "GRanges"
  )

  # Query each chromosome separately.
  # This avoids matrix-dimension inconsistencies when
  # summary() is called on multiple chromosome ranges at once.

  chr_width <- GenomicRanges::width(gr)

  chr_mean <- numeric(length(gr))
  chr_coverage <- numeric(length(gr))

  for (j in seq_along(gr)) {

    one_gr <- gr[j]

    mean_result <- rtracklayer::summary(
      bw,
      ranges = one_gr,
      size = 1L,
      type = "mean",
      defaultValue = 0,
      as = "matrix"
    )

    coverage_result <- rtracklayer::summary(
      bw,
      ranges = one_gr,
      size = 1L,
      type = "coverage",
      defaultValue = 0,
      as = "matrix"
    )

    chr_mean[j] <- as.numeric(mean_result)[1]
    chr_coverage[j] <- as.numeric(coverage_result)[1]
  }

  if (length(chr_mean) != length(chr_width)) {
    stop(
      paste(
        "Internal chromosome-summary mismatch for",
        f
      )
    )
  }

  global_mean <- sum(
    chr_mean * chr_width,
    na.rm = TRUE
  ) / sum(chr_width)

  weighted_coverage <- sum(
    chr_coverage * chr_width,
    na.rm = TRUE
  ) / sum(chr_width)

  total_reference_bases <- sum(chr_width)

  relative_total_signal <-
    global_mean * total_reference_bases

  mark <- basename(
    dirname(f)
  )

  sample_name <- basename(f)

  condition <- sub(
    "_.*$",
    "",
    sample_name
  )

  replicate <- ifelse(
    grepl("rep1", sample_name),
    "rep1",
    ifelse(
      grepl("rep2", sample_name),
      "rep2",
      NA_character_
    )
  )

  normalization_label <- ifelse(
    grepl(
      "NotNormalized",
      sample_name
    ),
    "NotNormalized",
    "Normalized"
  )

  results[[i]] <- data.frame(
    file = f,
    mark = mark,
    sample = sample_name,
    condition = condition,
    replicate = replicate,
    normalization_label =
      normalization_label,
    file_size_MB =
      file.info(f)$size / 1024^2,
    chromosomes_used =
      length(keep_chr),
    global_mean =
      global_mean,
    relative_total_signal =
      relative_total_signal,
    genome_coverage =
      weighted_coverage,
    stringsAsFactors = FALSE
  )
}

qc <- do.call(
  rbind,
  results
)

write.csv(
  qc,
  file.path(
    outdir,
    "12c_bigwig_global_qc.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("GLOBAL SIGNAL SUMMARY\n")
cat("========================================\n\n")

print(
  qc[
    ,
    c(
      "mark",
      "sample",
      "normalization_label",
      "global_mean",
      "relative_total_signal",
      "genome_coverage"
    )
  ],
  row.names = FALSE
)

cat("\n========================================\n")
cat("RPMI REP2 vs RPMI REP1 RATIOS\n")
cat("========================================\n\n")

marks <- unique(qc$mark)

ratio_list <- list()

for (m in marks) {

  x <- qc[qc$mark == m, ]

  r1 <- x[
    grepl("RPMI_rep1", x$sample),
  ]

  r2 <- x[
    grepl("RPMI_rep2", x$sample),
  ]

  if (
    nrow(r1) == 1 &&
    nrow(r2) == 1
  ) {

    ratio_list[[m]] <- data.frame(
      mark = m,
      RPMI_rep1_global_mean =
        r1$global_mean,
      RPMI_rep2_global_mean =
        r2$global_mean,
      rep2_over_rep1 =
        r2$global_mean /
        r1$global_mean,
      rep1_over_rep2_scaling_factor =
        r1$global_mean /
        r2$global_mean,
      stringsAsFactors = FALSE
    )
  }
}

ratios <- do.call(
  rbind,
  ratio_list
)

print(
  ratios,
  row.names = FALSE
)

write.csv(
  ratios,
  file.path(
    outdir,
    "12c_RPMI_rep2_scaling_audit.csv"
  ),
  row.names = FALSE
)

cat("\nSaved:\n")
cat(
  " ",
  file.path(
    outdir,
    "12c_bigwig_global_qc.csv"
  ),
  "\n"
)

cat(
  " ",
  file.path(
    outdir,
    "12c_RPMI_rep2_scaling_audit.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("12c COMPLETE\n")
cat("========================================\n")
