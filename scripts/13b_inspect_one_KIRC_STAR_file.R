# ============================================================
# 13b — INSPECT ONE TCGA-KIRC STAR-COUNTS FILE
#
# Purpose:
#   Confirm the current GDC STAR-count file structure before
#   extracting JAK3 across TCGA.
#
# Nothing is kept permanently except this console output.
# ============================================================

rm(list = ls())

metadata_file <-
  "results/13_TCGA_JAK3/13a_TCGA_STARCounts_metadata.csv"

if (!file.exists(metadata_file)) {
  stop("Missing 13a metadata file.")
}

meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("13b — INSPECT ONE TCGA-KIRC FILE\n")
cat("========================================\n\n")

kirc <- meta[
  meta$TCGA_project == "TCGA-KIRC" &
  meta$sample_type == "Primary Tumor",
]

cat("KIRC Primary Tumor records:", nrow(kirc), "\n")

if (nrow(kirc) == 0) {
  stop("No TCGA-KIRC Primary Tumor file found.")
}

one <- kirc[1, ]

cat("\nSelected record:\n")
print(one)

file_id <- one$file_id

if (
  is.na(file_id) ||
  file_id == ""
) {
  stop("file_id missing.")
}

url <- paste0(
  "https://api.gdc.cancer.gov/data/",
  file_id
)

tmp <- tempfile(
  fileext = ".tsv"
)

cat("\nDownloading temporary file:\n")
cat(url, "\n\n")

download.file(
  url,
  destfile = tmp,
  mode = "wb",
  quiet = FALSE
)

cat("\nTemporary file size:\n")
print(file.info(tmp)$size)

cat("\n========================================\n")
cat("FIRST 15 LINES\n")
cat("========================================\n\n")

lines <- readLines(
  tmp,
  n = 15
)

cat(
  paste(lines, collapse = "\n"),
  "\n"
)

cat("\n========================================\n")
cat("READING TABLE\n")
cat("========================================\n\n")

dat <- read.delim(
  tmp,
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Rows:", nrow(dat), "\n")
cat("Columns:", ncol(dat), "\n\n")

cat("Column names:\n")
print(colnames(dat))

cat("\n========================================\n")
cat("JAK3 SEARCH\n")
cat("========================================\n\n")

jak3_symbol <- which(
  dat == "JAK3",
  arr.ind = TRUE
)

jak3_ens <- which(
  grepl(
    "^ENSG00000105639",
    as.matrix(dat)
  ),
  arr.ind = TRUE
)

cat(
  "Cells matching symbol JAK3:",
  nrow(jak3_symbol),
  "\n"
)

cat(
  "Cells matching ENSG00000105639:",
  nrow(jak3_ens),
  "\n\n"
)

if (nrow(jak3_symbol) > 0) {

  rows <- unique(
    jak3_symbol[, "row"]
  )

  cat("JAK3 row:\n\n")

  print(
    dat[rows, , drop = FALSE],
    row.names = FALSE
  )

} else if (nrow(jak3_ens) > 0) {

  rows <- unique(
    jak3_ens[, "row"]
  )

  cat("JAK3 Ensembl row:\n\n")

  print(
    dat[rows, , drop = FALSE],
    row.names = FALSE
  )

} else {

  cat("JAK3 was not found in this file.\n")
}

unlink(tmp)

cat("\nTemporary file deleted.\n")

cat("\n========================================\n")
cat("13b COMPLETE\n")
cat("========================================\n")
