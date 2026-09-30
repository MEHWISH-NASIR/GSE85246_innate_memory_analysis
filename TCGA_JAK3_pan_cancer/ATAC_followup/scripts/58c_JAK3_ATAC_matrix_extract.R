# ============================================================
# 58c — Extract JAK3-region peaks from normalized TCGA ATAC matrix
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

rds_file <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/data/TCGA-ATAC_PanCan_Log2Norm_Counts.rds"

peak_file <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results/JAK3_ATAC_peaks_gene_plusminus50kb.csv"

outdir <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

priority_cancers <- c(
  "CHOL","KIRC","KIRP","THCA","HNSC",
  "STAD","LUAD","COAD","LIHC"
)

cat("\n========================================\n")
cat("58c JAK3 ATAC MATRIX EXTRACTION\n")
cat("========================================\n\n")

cat("Loading normalized RDS...\n")

x <- readRDS(rds_file)

cat("Object class:", paste(class(x), collapse = ", "), "\n")

if (!is.null(dim(x))) {
  cat("Dimensions:", nrow(x), "x", ncol(x), "\n")
}

cat("\nFirst columns:\n")
print(head(names(x), 12))

# ------------------------------------------------------------
# Read candidate JAK3 peaks
# ------------------------------------------------------------

jak <- read.csv(
  peak_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\nJAK3 candidate peaks:", nrow(jak), "\n")

# ------------------------------------------------------------
# Identify metadata columns
# ------------------------------------------------------------

meta_cols <- c(
  "seqnames",
  "start",
  "end",
  "name",
  "score",
  "annotation",
  "GC"
)

missing_meta <- setdiff(
  meta_cols,
  names(x)
)

if (length(missing_meta) > 0) {
  stop(
    paste(
      "Missing expected metadata columns:",
      paste(missing_meta, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# Match 46 JAK3 peaks by peak name
# ------------------------------------------------------------

idx <- match(
  jak$name,
  x$name
)

cat(
  "Matched JAK3 peaks:",
  sum(!is.na(idx)),
  "of",
  nrow(jak),
  "\n"
)

if (any(is.na(idx))) {
  cat("\nUnmatched peak names:\n")
  print(jak$name[is.na(idx)])
}

# ------------------------------------------------------------
# Identify sample columns for 9 priority cancers
#
# Matrix names may use either "-" or "_" after cancer code.
# ------------------------------------------------------------

sample_cols <- setdiff(
  names(x),
  meta_cols
)

cancer_from_col <- sub(
  "^([A-Za-z0-9]+)[_-].*$",
  "\\1",
  sample_cols
)

keep_sample <- cancer_from_col %in% priority_cancers

priority_cols <- sample_cols[keep_sample]
priority_cancer <- cancer_from_col[keep_sample]

cat("\n========================================\n")
cat("TECHNICAL REPLICATE COLUMNS\n")
cat("========================================\n\n")

print(
  table(
    factor(
      priority_cancer,
      levels = priority_cancers
    )
  )
)

cat(
  "\nTotal priority-cancer technical replicate columns:",
  length(priority_cols),
  "\n"
)

# ------------------------------------------------------------
# Extract only JAK3 peaks + priority-cancer ATAC columns
# ------------------------------------------------------------

jak3_matrix <- x[
  idx[!is.na(idx)],
  c(
    meta_cols,
    priority_cols
  ),
  drop = FALSE
]

outfile <-
  file.path(
    outdir,
    "JAK3_ATAC_normalized_46peaks_9cancers.rds"
  )

saveRDS(
  jak3_matrix,
  outfile
)

# Save metadata for sample columns
sample_map <- data.frame(
  matrix_column = priority_cols,
  Cancer = priority_cancer,
  stringsAsFactors = FALSE
)

write.csv(
  sample_map,
  file.path(
    outdir,
    "JAK3_ATAC_priority_sample_columns.csv"
  ),
  row.names = FALSE
)

cat("\nSaved:\n")
cat(" ", outfile, "\n")
cat(
  " ",
  file.path(
    outdir,
    "JAK3_ATAC_priority_sample_columns.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("58c COMPLETE\n")
cat("========================================\n")
