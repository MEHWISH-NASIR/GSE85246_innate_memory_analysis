library(dplyr)

expr <- read.csv(
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_final_TPM_metadata_complete.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

clinical <- read.csv(
  "TCGA_JAK3_pan_cancer/results/JAK3_TCGA_clinical_complete_33projects.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("JAK3 COMPLETE EXPRESSION + CLINICAL MERGE\n")
cat("========================================\n")

cat("\nExpression rows:", nrow(expr), "\n")
cat("Expression projects:", n_distinct(expr$project), "\n")

cat("Clinical rows:", nrow(clinical), "\n")
cat("Clinical projects:", n_distinct(clinical$project), "\n")

# ------------------------------------------------------------
# Clinical uniqueness
# ------------------------------------------------------------

dup <- clinical %>%
  count(project, submitter_id, name = "n") %>%
  filter(n > 1)

cat(
  "Duplicated clinical project/patient IDs:",
  nrow(dup),
  "\n"
)

if (nrow(dup) > 0) {
  stop("Clinical table contains duplicated project/patient IDs.")
}

# ------------------------------------------------------------
# Expression-clinical overlap before merge
# ------------------------------------------------------------

expr_patients <- expr %>%
  distinct(project, cases.submitter_id)

expr_keys <- paste(
  expr_patients$project,
  expr_patients$cases.submitter_id,
  sep = "|"
)

clinical_keys <- paste(
  clinical$project,
  clinical$submitter_id,
  sep = "|"
)

cat(
  "Unique expression patients:",
  nrow(expr_patients),
  "\n"
)

cat(
  "Expression patients with clinical match:",
  sum(expr_keys %in% clinical_keys),
  "/",
  nrow(expr_patients),
  "\n"
)

# ------------------------------------------------------------
# Retain useful clinical variables
# ------------------------------------------------------------

wanted <- c(
  "submitter_id",
  "project",
  "ajcc_pathologic_stage",
  "tumor_grade",
  "ajcc_pathologic_t",
  "ajcc_pathologic_n",
  "ajcc_pathologic_m",
  "vital_status",
  "days_to_death",
  "days_to_last_follow_up",
  "days_to_last_followup",
  "age_at_diagnosis",
  "gender"
)

clinical_small <- clinical %>%
  select(any_of(wanted))

cat("\nClinical fields retained:\n")
print(names(clinical_small))

# ------------------------------------------------------------
# Merge by project + patient
# ------------------------------------------------------------

merged <- expr %>%
  left_join(
    clinical_small,
    by = c(
      "cases.submitter_id" = "submitter_id",
      "project" = "project"
    )
  )

cat("\n========================================\n")
cat("MERGE AUDIT\n")
cat("========================================\n")

cat("Rows before merge:", nrow(expr), "\n")
cat("Rows after merge:", nrow(merged), "\n")
cat("Row difference:", nrow(merged) - nrow(expr), "\n")
cat("Projects after merge:", n_distinct(merged$project), "\n")

cat(
  "Unique expression patients after merge:",
  n_distinct(
    paste(
      merged$project,
      merged$cases.submitter_id,
      sep = "|"
    )
  ),
  "\n"
)

if (nrow(merged) != nrow(expr)) {
  stop("Clinical join changed expression row count.")
}

if (n_distinct(merged$project) != 33) {
  stop(
    "Expected 33 projects after merge; found ",
    n_distinct(merged$project)
  )
}

# ------------------------------------------------------------
# Clinical availability
# ------------------------------------------------------------

if ("ajcc_pathologic_stage" %in% names(merged)) {
  cat(
    "Rows with AJCC stage:",
    sum(
      !is.na(merged$ajcc_pathologic_stage) &
      merged$ajcc_pathologic_stage != ""
    ),
    "\n"
  )
}

if ("tumor_grade" %in% names(merged)) {
  cat(
    "Rows with tumor grade:",
    sum(
      !is.na(merged$tumor_grade) &
      merged$tumor_grade != ""
    ),
    "\n"
  )
}

if ("vital_status" %in% names(merged)) {
  cat(
    "Rows with vital status:",
    sum(
      !is.na(merged$vital_status) &
      merged$vital_status != ""
    ),
    "\n"
  )
}

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

outfile <- paste0(
  "TCGA_JAK3_pan_cancer/results/expression/",
  "JAK3_expression_clinical_complete_33projects.csv"
)

write.csv(
  merged,
  outfile,
  row.names = FALSE
)

cat("\nSaved:\n", outfile, "\n")
cat("\nCOMPLETE EXPRESSION-CLINICAL MERGE FINISHED\n")
