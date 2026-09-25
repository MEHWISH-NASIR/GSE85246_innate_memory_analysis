library(dplyr)

# ============================================================
# COMPLETE JAK3 EXPRESSION + CLINICAL DATA
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

    sample_barcode =
      substr(sample.submitter_id, 1, 16),

    cancer =
      sub("^TCGA-", "", project)
  ) %>%
  filter(primary_malignant)

# ============================================================
# PANIMMUNE LEUKOCYTE FRACTION
# ============================================================

leuk <- read.delim(
  paste0(
    "TCGA_JAK3_pan_cancer/data/PanImmune/",
    "TCGA_all_leuk_estimate.masked.20170107.tsv"
  ),
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(leuk) <- c(
  "cancer",
  "panimmune_barcode",
  "leukocyte_fraction"
)

leuk <- leuk %>%
  mutate(
    sample_barcode =
      substr(panimmune_barcode, 1, 16)
  )

cat("\n========================================\n")
cat("JAK3 PANIMMUNE MATCH AUDIT\n")
cat("========================================\n")

cat("\nExpression primary-malignant rows:",
    nrow(expr), "\n")

cat("Expression projects:",
    n_distinct(expr$project), "\n")

cat("\nPanImmune rows:",
    nrow(leuk), "\n")

cat("PanImmune cancers:",
    n_distinct(leuk$cancer), "\n")

# ============================================================
# DUPLICATE PANIMMUNE SAMPLE-BARCODE AUDIT
# ============================================================

leuk_dup <- leuk %>%
  count(
    cancer,
    sample_barcode,
    name = "records"
  ) %>%
  filter(records > 1)

cat(
  "\nDuplicated PanImmune cancer/sample barcodes:",
  nrow(leuk_dup),
  "\n"
)

# Collapse duplicate records mapping to same 16-char sample barcode
leuk_clean <- leuk %>%
  group_by(
    cancer,
    sample_barcode
  ) %>%
  summarise(
    leukocyte_fraction =
      mean(
        leukocyte_fraction,
        na.rm = TRUE
      ),

    panimmune_records = n(),

    .groups = "drop"
  )

# ============================================================
# SAMPLE-LEVEL MATCH
# ============================================================

expr <- expr %>%
  mutate(
    immune_match =
      paste(cancer, sample_barcode) %in%
      paste(
        leuk_clean$cancer,
        leuk_clean$sample_barcode
      )
  )

sample_audit <- expr %>%
  group_by(project) %>%
  summarise(
    Expression_samples = n(),

    Expression_patients =
      n_distinct(cases.submitter_id),

    PanImmune_matched_samples =
      sum(immune_match),

    PanImmune_unmatched_samples =
      sum(!immune_match),

    Sample_match_percent =
      round(
        100 *
        PanImmune_matched_samples /
        Expression_samples,
        2
      ),

    .groups = "drop"
  ) %>%
  arrange(project)

# ============================================================
# ACTUAL INNER JOIN
# ============================================================

matched <- expr %>%
  inner_join(
    leuk_clean,
    by = c(
      "cancer",
      "sample_barcode"
    )
  )

# ============================================================
# COLLAPSE TO PATIENT LEVEL
# ============================================================

patient <- matched %>%
  group_by(
    project,
    cases.submitter_id
  ) %>%
  summarise(
    JAK3_TPM =
      mean(TPM, na.rm = TRUE),

    JAK3_log2TPM1 =
      mean(
        log2(TPM + 1),
        na.rm = TRUE
      ),

    leukocyte_fraction =
      mean(
        leukocyte_fraction,
        na.rm = TRUE
      ),

    stage = {
      z <- unique(
        na.omit(
          ajcc_pathologic_stage
        )
      )

      if (length(z) == 0)
        NA_character_
      else
        z[1]
    },

    grade = {
      z <- unique(
        na.omit(tumor_grade)
      )

      if (length(z) == 0)
        NA_character_
      else
        z[1]
    },

    expression_samples = n(),

    .groups = "drop"
  )

patient_audit <- patient %>%
  count(
    project,
    name = "PanImmune_matched_patients"
  )

final_audit <- sample_audit %>%
  left_join(
    patient_audit,
    by = "project"
  ) %>%
  mutate(
    PanImmune_matched_patients =
      coalesce(
        PanImmune_matched_patients,
        0L
      ),

    Patient_match_percent =
      round(
        100 *
        PanImmune_matched_patients /
        Expression_patients,
        2
      )
  )

# ============================================================
# IDENTIFY PROJECTS WITHOUT PANIMMUNE LF
# ============================================================

missing_projects <- setdiff(
  sort(unique(expr$project)),
  paste0(
    "TCGA-",
    sort(unique(leuk$cancer))
  )
)

# ============================================================
# SAVE
# ============================================================

write.csv(
  final_audit,
  paste0(
    "TCGA_JAK3_pan_cancer/results/",
    "JAK3_PanImmune_match_audit_complete.csv"
  ),
  row.names = FALSE
)

write.csv(
  patient,
  paste0(
    "TCGA_JAK3_pan_cancer/results/",
    "JAK3_PanImmune_patient_data_complete.csv"
  ),
  row.names = FALSE
)

write.csv(
  leuk_dup,
  paste0(
    "TCGA_JAK3_pan_cancer/results/",
    "JAK3_PanImmune_duplicate_sample_barcodes.csv"
  ),
  row.names = FALSE
)

# ============================================================
# CONSOLE
# ============================================================

cat("\n===== PROJECT MATCH COVERAGE =====\n")

print(
  as.data.frame(final_audit),
  row.names = FALSE
)

cat("\n===== PROJECTS WITHOUT PANIMMUNE LF =====\n")

print(missing_projects)

cat("\n===== OVERALL MATCH =====\n")

cat(
  "Matched expression samples:",
  nrow(matched),
  "/",
  nrow(expr),
  "\n"
)

cat(
  "Matched patients:",
  n_distinct(
    paste(
      patient$project,
      patient$cases.submitter_id
    )
  ),
  "/",
  n_distinct(
    paste(
      expr$project,
      expr$cases.submitter_id
    )
  ),
  "\n"
)

cat(
  "Projects with matched PanImmune patients:",
  n_distinct(patient$project),
  "\n"
)

cat("\n52a PANIMMUNE MATCH AUDIT COMPLETE\n")
