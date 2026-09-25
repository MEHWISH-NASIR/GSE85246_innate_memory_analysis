library(dplyr)

d <- read.csv(
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_expression_clinical_complete_33projects.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

first_nonmissing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & trimws(x) != ""]
  if (length(x) == 0) NA_character_ else x[1]
}

clean_missing <- function(x) {
  x <- trimws(as.character(x))

  bad <- tolower(x) %in% c(
    "",
    "na",
    "n/a",
    "unknown",
    "not reported",
    "not available",
    "[not available]",
    "--"
  )

  x[bad] <- NA_character_
  x
}

# ------------------------------------------------------------
# Canonical primary malignant samples
# ------------------------------------------------------------

d <- d %>%
  mutate(
    primary_malignant = case_when(

      project == "TCGA-LAML" &
        sample_type ==
        "Primary Blood Derived Cancer - Peripheral Blood" ~ TRUE,

      project != "TCGA-LAML" &
        sample_type == "Primary Tumor" ~ TRUE,

      TRUE ~ FALSE
    )
  )

# ------------------------------------------------------------
# Patient-level collapse
# ------------------------------------------------------------

patient <- d %>%
  filter(primary_malignant) %>%
  group_by(project, cases.submitter_id) %>%
  summarise(
    TPM = mean(TPM, na.rm = TRUE),

    stage_raw =
      first_nonmissing(ajcc_pathologic_stage),

    grade_raw =
      first_nonmissing(tumor_grade),

    expression_samples = n(),

    .groups = "drop"
  ) %>%
  mutate(
    stage_raw = clean_missing(stage_raw),
    grade_raw = clean_missing(grade_raw),

    stage_upper = toupper(stage_raw),

    # Order is important: IV before III before II before I
    stage_simple = case_when(
      grepl("^STAGE IV", stage_upper)  ~ "Stage IV",
      grepl("^STAGE III", stage_upper) ~ "Stage III",
      grepl("^STAGE II", stage_upper)  ~ "Stage II",
      grepl("^STAGE I", stage_upper)   ~ "Stage I",
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------
# Raw stage labels
# ------------------------------------------------------------

stage_raw_counts <- patient %>%
  filter(!is.na(stage_raw)) %>%
  count(
    project,
    stage_raw,
    stage_simple,
    name = "patients"
  ) %>%
  arrange(project, stage_simple, stage_raw)

# ------------------------------------------------------------
# Simplified stage counts
# ------------------------------------------------------------

stage_simple_counts <- patient %>%
  filter(!is.na(stage_simple)) %>%
  count(
    project,
    stage_simple,
    name = "patients"
  ) %>%
  arrange(project, stage_simple)

# ------------------------------------------------------------
# Raw grade labels
# ------------------------------------------------------------

grade_counts <- patient %>%
  filter(!is.na(grade_raw)) %>%
  count(
    project,
    grade_raw,
    name = "patients"
  ) %>%
  arrange(project, grade_raw)

# ------------------------------------------------------------
# Project-level availability
# ------------------------------------------------------------

availability <- patient %>%
  group_by(project) %>%
  summarise(
    Primary_patients = n(),

    Stage_available_n =
      sum(!is.na(stage_simple)),

    Stage_groups =
      n_distinct(
        stage_simple[!is.na(stage_simple)]
      ),

    Grade_available_n =
      sum(!is.na(grade_raw)),

    Grade_raw_groups =
      n_distinct(
        grade_raw[!is.na(grade_raw)]
      ),

    .groups = "drop"
  ) %>%
  arrange(project)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  patient,
  "TCGA_JAK3_pan_cancer/results/JAK3_stage_grade_patient_level.csv",
  row.names = FALSE
)

write.csv(
  stage_raw_counts,
  "TCGA_JAK3_pan_cancer/results/JAK3_stage_raw_labels.csv",
  row.names = FALSE
)

write.csv(
  stage_simple_counts,
  "TCGA_JAK3_pan_cancer/results/JAK3_stage_simple_counts.csv",
  row.names = FALSE
)

write.csv(
  grade_counts,
  "TCGA_JAK3_pan_cancer/results/JAK3_grade_raw_labels.csv",
  row.names = FALSE
)

write.csv(
  availability,
  "TCGA_JAK3_pan_cancer/results/JAK3_stage_grade_availability_final.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("JAK3 STAGE / GRADE LABEL AUDIT\n")
cat("========================================\n")

cat("\n===== STAGE SIMPLE COUNTS =====\n")
print(
  as.data.frame(stage_simple_counts),
  row.names = FALSE
)

cat("\n===== RAW GRADE LABELS =====\n")
print(
  as.data.frame(grade_counts),
  row.names = FALSE
)

cat("\n===== PROJECT AVAILABILITY =====\n")
print(
  as.data.frame(availability),
  row.names = FALSE
)

cat("\n50a STAGE/GRADE LABEL AUDIT COMPLETE\n")
