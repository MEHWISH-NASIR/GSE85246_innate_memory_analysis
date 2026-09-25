library(dplyr)

# ============================================================
# JAK3 PRIMARY MALIGNANT EXPRESSION
# ============================================================

expr <- read.csv(
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_expression_clinical_complete_33projects.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr <- expr %>%
  mutate(
    primary_malignant = case_when(
      project == "TCGA-LAML" &
        sample_type ==
        "Primary Blood Derived Cancer - Peripheral Blood" ~ TRUE,

      project != "TCGA-LAML" &
        sample_type == "Primary Tumor" ~ TRUE,

      TRUE ~ FALSE
    ),

    cancer = sub("^TCGA-", "", project),

    sample_barcode =
      substr(sample.submitter_id, 1, 16)
  ) %>%
  filter(primary_malignant)

# ============================================================
# CIBERSORT
# ============================================================

ciber <- read.delim(
  "TCGA_JAK3_pan_cancer/data/PanImmune/TCGA.Kallisto.fullIDs.cibersort.relative.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cell_cols <- names(ciber)[3:24]

ciber <- ciber %>%
  mutate(
    SampleID_hyphen =
      gsub("\\.", "-", SampleID),

    sample_barcode =
      substr(SampleID_hyphen, 1, 16),

    quality_pass =
      !is.na(P.value) &
      P.value < 0.05
  )

cat("\n========================================\n")
cat("JAK3 CIBERSORT MATCH AUDIT\n")
cat("========================================\n")

cat("\nExpression primary-malignant samples:",
    nrow(expr), "\n")

cat("Expression patients:",
    n_distinct(
      paste(
        expr$project,
        expr$cases.submitter_id
      )
    ),
    "\n")

cat("\nCIBERSORT rows:", nrow(ciber), "\n")
cat("CIBERSORT cancers:",
    n_distinct(ciber$CancerType), "\n")

cat("CIBERSORT P < 0.05:",
    sum(ciber$quality_pass), "\n")

# ============================================================
# DUPLICATE AUDIT
# ============================================================

full_id_dup <- ciber %>%
  count(
    CancerType,
    SampleID,
    name = "records"
  ) %>%
  filter(records > 1)

barcode_dup <- ciber %>%
  count(
    CancerType,
    sample_barcode,
    name = "records"
  ) %>%
  filter(records > 1)

cat(
  "\nDuplicated exact CIBERSORT SampleIDs:",
  nrow(full_id_dup),
  "\n"
)

cat(
  "Duplicated 16-char sample barcodes:",
  nrow(barcode_dup),
  "\n"
)

write.csv(
  full_id_dup,
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_duplicate_fullIDs.csv",
  row.names = FALSE
)

write.csv(
  barcode_dup,
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_duplicate_sample_barcodes.csv",
  row.names = FALSE
)

# ============================================================
# CREATE ONE CIBERSORT RECORD PER SAMPLE BARCODE
#
# If multiple records map to one TCGA sample:
#   1. lower P-value preferred
#   2. higher correlation preferred
#   3. lower RMSE preferred
# ============================================================

ciber_clean <- ciber %>%
  arrange(
    CancerType,
    sample_barcode,
    P.value,
    desc(Correlation),
    RMSE
  ) %>%
  group_by(
    CancerType,
    sample_barcode
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

cat(
  "\nUnique CIBERSORT cancer/sample barcodes after cleanup:",
  nrow(ciber_clean),
  "\n"
)

cat(
  "Quality-passing unique samples:",
  sum(ciber_clean$quality_pass),
  "\n"
)

# ============================================================
# MATCH ALL CIBERSORT SAMPLES
# ============================================================

all_match <- expr %>%
  inner_join(
    ciber_clean,
    by = c(
      "cancer" = "CancerType",
      "sample_barcode" = "sample_barcode"
    )
  )

# ============================================================
# QUALITY-FILTERED MATCH
# ============================================================

quality_match <- all_match %>%
  filter(quality_pass)

# ============================================================
# PROJECT-LEVEL COVERAGE
# ============================================================

all_counts <- all_match %>%
  group_by(project) %>%
  summarise(
    CIBERSORT_matched_samples = n(),
    CIBERSORT_matched_patients =
      n_distinct(cases.submitter_id),
    .groups = "drop"
  )

quality_counts <- quality_match %>%
  group_by(project) %>%
  summarise(
    CIBERSORT_P005_samples = n(),
    CIBERSORT_P005_patients =
      n_distinct(cases.submitter_id),
    .groups = "drop"
  )

audit <- expr %>%
  group_by(project) %>%
  summarise(
    Expression_samples = n(),
    Expression_patients =
      n_distinct(cases.submitter_id),
    .groups = "drop"
  ) %>%

  left_join(
    all_counts,
    by = "project"
  ) %>%

  left_join(
    quality_counts,
    by = "project"
  ) %>%

  mutate(
    across(
      c(
        CIBERSORT_matched_samples,
        CIBERSORT_matched_patients,
        CIBERSORT_P005_samples,
        CIBERSORT_P005_patients
      ),
      ~ coalesce(.x, 0L)
    ),

    Match_percent =
      round(
        100 *
        CIBERSORT_matched_patients /
        Expression_patients,
        2
      ),

    Quality_match_percent =
      round(
        100 *
        CIBERSORT_P005_patients /
        Expression_patients,
        2
      ),

    Low_quality_filtered_N =
      CIBERSORT_P005_patients < 30
  ) %>%
  arrange(project)

# ============================================================
# SAVE
# ============================================================

write.csv(
  audit,
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_match_audit.csv",
  row.names = FALSE
)

write.csv(
  quality_match,
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_quality_matched_samples.csv",
  row.names = FALSE
)

# ============================================================
# CONSOLE
# ============================================================

cat("\n===== CIBERSORT PROJECT COVERAGE =====\n")

print(
  as.data.frame(audit),
  row.names = FALSE
)

cat("\n===== OVERALL =====\n")

cat(
  "Matched samples before quality filter:",
  nrow(all_match),
  "/",
  nrow(expr),
  "\n"
)

cat(
  "Matched samples with CIBERSORT P < 0.05:",
  nrow(quality_match),
  "/",
  nrow(expr),
  "\n"
)

cat(
  "Projects represented before filter:",
  n_distinct(all_match$project),
  "\n"
)

cat(
  "Projects represented after P < 0.05:",
  n_distinct(quality_match$project),
  "\n"
)

cat(
  "Projects with >=30 quality-filtered patients:",
  sum(audit$CIBERSORT_P005_patients >= 30),
  "\n"
)

cat("\n52c CIBERSORT MATCH AUDIT COMPLETE\n")
