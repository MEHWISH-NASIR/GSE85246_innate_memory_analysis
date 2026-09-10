# ============================================================
# 12i FIX — robust exact coordinate matching
# ============================================================

rm(list = ls())

outdir <- "results/12_kinase_reevaluation"

old_file <- file.path(
  outdir,
  "12d_hg19_TSS_2kb_promoters.csv"
)

new_file <- file.path(
  outdir,
  "12i_target_transcript_promoters.csv"
)

effects_file <- file.path(
  outdir,
  "12i_target_transcript_TSS_effects.csv"
)

old <- read.csv(
  old_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

new <- read.csv(
  new_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# Standardise fields
# ------------------------------------------------------------

old$SYMBOL <- trimws(as.character(old$SYMBOL))
new$SYMBOL <- trimws(as.character(new$SYMBOL))

old$chr <- trimws(as.character(old$chr))
new$chr <- trimws(as.character(new$chr))

old$start <- as.integer(old$start)
old$end <- as.integer(old$end)

new$promoter_start <- as.integer(new$promoter_start)
new$promoter_end <- as.integer(new$promoter_end)

# ------------------------------------------------------------
# Exact coordinate keys
# ------------------------------------------------------------

old_key <- paste(
  old$SYMBOL,
  old$chr,
  old$start,
  old$end,
  sep = "|"
)

new_key <- paste(
  new$SYMBOL,
  new$chr,
  new$promoter_start,
  new$promoter_end,
  sep = "|"
)

exact_idx <- match(
  new_key,
  old_key
)

# Gene-level lookup for old coordinates
gene_idx <- match(
  new$SYMBOL,
  old$SYMBOL
)

new$Matches_12d_gene_promoter <-
  !is.na(exact_idx)

new$Old_12d_start <-
  old$start[gene_idx]

new$Old_12d_end <-
  old$end[gene_idx]

# ------------------------------------------------------------
# Validate target genes
# ------------------------------------------------------------

targets <- c(
  "JAK3",
  "EPHB2",
  "MET",
  "MAP3K8",
  "BMPR1A"
)

target_matches <- new[
  new$SYMBOL %in% targets &
  new$Matches_12d_gene_promoter,
]

cat("\n========================================\n")
cat("EXACT 12d PROMOTER MATCHES\n")
cat("========================================\n\n")

print(
  target_matches[
    ,
    c(
      "SYMBOL",
      "Promoter_ID",
      "chr",
      "promoter_start",
      "promoter_end",
      "TSS",
      "tx_name"
    )
  ],
  row.names = FALSE
)

cat(
  "\nTarget exact matches:",
  nrow(target_matches),
  "\n"
)

if (nrow(target_matches) != 5) {
  stop(
    paste(
      "Expected 5 target exact matches but found",
      nrow(target_matches)
    )
  )
}

# ------------------------------------------------------------
# Save corrected promoter manifest
# ------------------------------------------------------------

write.csv(
  new,
  new_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# Also repair corresponding annotation columns in effects file
# ------------------------------------------------------------

if (file.exists(effects_file)) {

  e <- read.csv(
    effects_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  idx <- match(
    e$Promoter_ID,
    new$Promoter_ID
  )

  if (anyNA(idx)) {
    stop(
      "Some effect rows could not be matched by Promoter_ID."
    )
  }

  e$Matches_12d_gene_promoter <-
    new$Matches_12d_gene_promoter[idx]

  e$Old_12d_start <-
    new$Old_12d_start[idx]

  e$Old_12d_end <-
    new$Old_12d_end[idx]

  write.csv(
    e,
    effects_file,
    row.names = FALSE
  )
}

cat("\nCorrected files saved.\n")

cat("\n========================================\n")
cat("12i MATCH FIX COMPLETE\n")
cat("========================================\n")
