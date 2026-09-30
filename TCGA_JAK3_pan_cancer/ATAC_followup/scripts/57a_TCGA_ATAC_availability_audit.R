# ============================================================
# 57a — TCGA ATAC availability audit for prioritized JAK3 tumors
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(
  library(readxl)
)

xlsx <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/data/TCGA-ATAC_DataS1_DonorsAndStats_v4.xlsx"

ranking_file <-
  "TCGA_JAK3_pan_cancer/results/release/JAK3_tumor_priority_ranking.csv"

outdir <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

priority_cancers <- c(
  "CHOL",
  "KIRC",
  "KIRP",
  "THCA",
  "HNSC",
  "STAD",
  "LUAD",
  "COAD",
  "LIHC"
)

# ------------------------------------------------------------
# Read first sheet without assuming where header starts
# ------------------------------------------------------------

raw <- read_excel(
  xlsx,
  sheet = "Cancer Types and Pass-Fail Rate",
  col_names = FALSE
)

# Find row containing actual "Abbreviation" header
header_row <- which(
  apply(
    raw,
    1,
    function(z) {
      any(
        trimws(as.character(z)) == "Abbreviation",
        na.rm = TRUE
      )
    }
  )
)[1]

if (is.na(header_row)) {
  stop("Could not find the Abbreviation header row.")
}

cat("\n========================================\n")
cat("TCGA ATAC AVAILABILITY AUDIT\n")
cat("========================================\n\n")

cat("Detected header row:", header_row, "\n\n")

# Re-read using correct header
tab <- read_excel(
  xlsx,
  sheet = "Cancer Types and Pass-Fail Rate",
  skip = header_row - 1
)

tab <- as.data.frame(
  tab,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Parsed columns:\n\n")
print(names(tab))

# ------------------------------------------------------------
# Identify abbreviation column
# ------------------------------------------------------------

abbr_col <- grep(
  "Abbreviation",
  names(tab),
  ignore.case = TRUE,
  value = TRUE
)[1]

if (is.na(abbr_col)) {
  stop("Abbreviation column not found after parsing.")
}

tab[[abbr_col]] <- trimws(
  as.character(tab[[abbr_col]])
)

# ------------------------------------------------------------
# Show all donor/sample-related columns
# ------------------------------------------------------------

count_cols <- grep(
  "sample|donor",
  names(tab),
  ignore.case = TRUE,
  value = TRUE
)

cat("\nDonor/sample-related columns:\n\n")
print(count_cols)

# ------------------------------------------------------------
# Extract nine JAK3-priority cancers
# ------------------------------------------------------------

audit <- tab[
  tab[[abbr_col]] %in% priority_cancers,
  unique(
    c(
      abbr_col,
      count_cols
    )
  ),
  drop = FALSE
]

# Preserve RNA priority ordering
audit$RNA_priority_rank <- match(
  audit[[abbr_col]],
  priority_cancers
)

audit <- audit[
  order(audit$RNA_priority_rank),
]

audit <- audit[
  c(
    "RNA_priority_rank",
    setdiff(
      names(audit),
      "RNA_priority_rank"
    )
  )
]

cat("\n========================================\n")
cat("NINE PRIORITIZED JAK3 CANCERS\n")
cat("========================================\n\n")

print(
  audit,
  row.names = FALSE
)

write.csv(
  audit,
  file.path(
    outdir,
    "JAK3_priority_TCGA_ATAC_availability_raw.csv"
  ),
  row.names = FALSE
)

cat("\nSaved:\n")
cat(
  " ",
  file.path(
    outdir,
    "JAK3_priority_TCGA_ATAC_availability_raw.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("57a COMPLETE\n")
cat("========================================\n")
