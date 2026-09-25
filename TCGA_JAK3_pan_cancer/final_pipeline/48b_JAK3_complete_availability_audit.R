library(dplyr)

infile <- paste0(
  "TCGA_JAK3_pan_cancer/results/expression/",
  "JAK3_expression_clinical_complete_33projects.csv"
)

d <- read.csv(
  infile,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & trimws(as.character(x)) != ""]
  if (length(x) == 0) NA_character_ else as.character(x[1])
}

first_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else x[1]
}

clean_clinical_text <- function(x) {
  x <- trimws(as.character(x))

  bad <- tolower(x) %in% c(
    "",
    "na",
    "n/a",
    "not reported",
    "not available",
    "unknown",
    "--",
    "[not available]"
  )

  x[bad] <- NA_character_
  x
}

# ------------------------------------------------------------
# Define canonical analysis sample classes
# ------------------------------------------------------------

d <- d %>%
  mutate(
    analysis_class = case_when(

      project == "TCGA-LAML" &
        sample_type ==
        "Primary Blood Derived Cancer - Peripheral Blood" ~
        "Primary malignant",

      project != "TCGA-LAML" &
        sample_type == "Primary Tumor" ~
        "Primary malignant",

      sample_type == "Solid Tissue Normal" ~
        "Normal",

      TRUE ~ "Other"
    )
  )

# ------------------------------------------------------------
# Collapse primary malignant samples to patient level
# ------------------------------------------------------------

primary <- d %>%
  filter(analysis_class == "Primary malignant") %>%
  group_by(project, cases.submitter_id) %>%
  summarise(

    JAK3_TPM = mean(TPM, na.rm = TRUE),

    stage = first_nonmissing(
      ajcc_pathologic_stage
    ),

    grade = first_nonmissing(
      tumor_grade
    ),

    vital_status = first_nonmissing(
      vital_status
    ),

    days_to_death = first_numeric(
      days_to_death
    ),

    days_to_last_follow_up = first_numeric(
      days_to_last_follow_up
    ),

    expression_samples = n(),

    .groups = "drop"
  ) %>%
  mutate(
    stage = clean_clinical_text(stage),
    grade = clean_clinical_text(grade),

    stage_upper = toupper(stage),

    is_stage_I = !is.na(stage_upper) &
      grepl(
        "^STAGE I([ABC])?$",
        stage_upper
      ),

    survival_time = case_when(
      toupper(vital_status) == "DEAD" ~
        days_to_death,

      toupper(vital_status) == "ALIVE" ~
        days_to_last_follow_up,

      TRUE ~ NA_real_
    ),

    survival_ready =
      !is.na(vital_status) &
      !is.na(survival_time) &
      survival_time >= 0
  )

# ------------------------------------------------------------
# Collapse normals to patient level
# ------------------------------------------------------------

normal <- d %>%
  filter(analysis_class == "Normal") %>%
  group_by(project, cases.submitter_id) %>%
  summarise(
    JAK3_TPM = mean(TPM, na.rm = TRUE),
    normal_expression_samples = n(),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Core counts
# ------------------------------------------------------------

primary_counts <- primary %>%
  count(
    project,
    name = "Primary_malignant_patients"
  )

normal_counts <- normal %>%
  count(
    project,
    name = "Normal_patients"
  )

matched_counts <- inner_join(
  primary %>%
    select(project, cases.submitter_id),

  normal %>%
    select(project, cases.submitter_id),

  by = c(
    "project",
    "cases.submitter_id"
  )
) %>%
  count(
    project,
    name = "Matched_pairs"
  )

# ------------------------------------------------------------
# Stage availability
# ------------------------------------------------------------

stage_summary <- primary %>%
  group_by(project) %>%
  summarise(

    Patients_with_stage =
      sum(!is.na(stage)),

    Stage_I_patients =
      sum(is_stage_I, na.rm = TRUE),

    Stage_groups =
      n_distinct(
        stage[!is.na(stage)]
      ),

    .groups = "drop"
  )

# ------------------------------------------------------------
# Grade availability
# ------------------------------------------------------------

grade_summary <- primary %>%
  group_by(project) %>%
  summarise(

    Patients_with_grade =
      sum(!is.na(grade)),

    Grade_groups =
      n_distinct(
        grade[!is.na(grade)]
      ),

    .groups = "drop"
  )

# ------------------------------------------------------------
# Survival availability
# ------------------------------------------------------------

survival_summary <- primary %>%
  group_by(project) %>%
  summarise(

    Patients_with_vital_status =
      sum(!is.na(vital_status)),

    Survival_ready_patients =
      sum(survival_ready),

    Deaths_ready =
      sum(
        survival_ready &
        toupper(vital_status) == "DEAD"
      ),

    .groups = "drop"
  )

# ------------------------------------------------------------
# Build full 33-project audit
# ------------------------------------------------------------

projects <- data.frame(
  project = sort(unique(d$project)),
  stringsAsFactors = FALSE
)

audit <- projects %>%

  left_join(
    primary_counts,
    by = "project"
  ) %>%

  left_join(
    normal_counts,
    by = "project"
  ) %>%

  left_join(
    matched_counts,
    by = "project"
  ) %>%

  left_join(
    stage_summary,
    by = "project"
  ) %>%

  left_join(
    grade_summary,
    by = "project"
  ) %>%

  left_join(
    survival_summary,
    by = "project"
  ) %>%

  mutate(
    across(
      c(
        Primary_malignant_patients,
        Normal_patients,
        Matched_pairs,
        Patients_with_stage,
        Stage_I_patients,
        Stage_groups,
        Patients_with_grade,
        Grade_groups,
        Patients_with_vital_status,
        Survival_ready_patients,
        Deaths_ready
      ),
      ~ coalesce(.x, 0L)
    ),

    Tumor_vs_normal_available =
      Primary_malignant_patients > 0 &
      Normal_patients > 0,

    Matched_available =
      Matched_pairs > 0,

    Early_stage_vs_normal_available =
      Stage_I_patients > 0 &
      Normal_patients > 0,

    Stage_trend_available =
      Patients_with_stage > 0 &
      Stage_groups >= 2,

    Grade_trend_available =
      Patients_with_grade > 0 &
      Grade_groups >= 2,

    Prevalence_available =
      Normal_patients > 0,

    Survival_available =
      Survival_ready_patients > 0,

    Small_normal_n =
      Normal_patients > 0 &
      Normal_patients < 10,

    Small_matched_n =
      Matched_pairs > 0 &
      Matched_pairs < 10,

    Small_stageI_n =
      Stage_I_patients > 0 &
      Stage_I_patients < 10
  ) %>%
  arrange(project)

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

cat("\n========================================\n")
cat("JAK3 COMPLETE 33-CANCER AVAILABILITY AUDIT\n")
cat("========================================\n\n")

print(
  as.data.frame(audit),
  row.names = FALSE
)

cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n")

cat("TCGA projects:", nrow(audit), "\n")

cat(
  "Tumor-normal available:",
  sum(audit$Tumor_vs_normal_available),
  "\n"
)

cat(
  "Matched-pair available:",
  sum(audit$Matched_available),
  "\n"
)

cat(
  "Early-stage vs normal available:",
  sum(audit$Early_stage_vs_normal_available),
  "\n"
)

cat(
  "Stage trend available:",
  sum(audit$Stage_trend_available),
  "\n"
)

cat(
  "Grade trend available:",
  sum(audit$Grade_trend_available),
  "\n"
)

cat(
  "Prevalence available:",
  sum(audit$Prevalence_available),
  "\n"
)

cat(
  "Survival available:",
  sum(audit$Survival_available),
  "\n"
)

write.csv(
  audit,
  paste0(
    "TCGA_JAK3_pan_cancer/results/",
    "JAK3_complete_33cancer_analysis_availability.csv"
  ),
  row.names = FALSE
)

cat("\n48b COMPLETE AVAILABILITY AUDIT FINISHED\n")
