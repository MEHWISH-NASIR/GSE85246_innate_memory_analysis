library(readxl)
library(dplyr)

# ------------------------------------------------------------
# JAK3 patient-level primary malignant expression
# ------------------------------------------------------------

expr <- read.csv(
  "TCGA_JAK3_pan_cancer/results/JAK3_primary_tumor_patient_level.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    log2TPM1 = log2(TPM + 1)
  )

# ------------------------------------------------------------
# TCGA-CDR
# ------------------------------------------------------------

cdr_file <- list.files(
  "D:/TCGA_CDR",
  pattern = "\\.xlsx$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(cdr_file) != 1) {
  stop("Expected exactly one TCGA-CDR workbook.")
}

cdr <- read_excel(
  cdr_file[1],
  sheet = "TCGA-CDR"
) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  mutate(
    project = paste0("TCGA-", type),

    Redaction_clean = trimws(as.character(Redaction)),
    Redacted =
      !is.na(Redaction_clean) &
      Redaction_clean != "",

    across(
      c(OS, OS.time,
        DSS, DSS.time,
        DFI, DFI.time,
        PFI, PFI.time),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

# ------------------------------------------------------------
# CDR uniqueness
# ------------------------------------------------------------

dup <- cdr %>%
  count(project, bcr_patient_barcode) %>%
  filter(n > 1)

cat("\n========================================\n")
cat("TCGA-CDR MATCH AUDIT\n")
cat("========================================\n")

cat("\nExpression patients:", nrow(expr), "\n")
cat("Expression projects:", n_distinct(expr$project), "\n")

cat("CDR rows:", nrow(cdr), "\n")
cat("CDR projects:", n_distinct(cdr$project), "\n")
cat("Duplicated CDR project/patient IDs:", nrow(dup), "\n")
cat("Redacted CDR cases:", sum(cdr$Redacted), "\n")

if (nrow(dup) > 0) {
  stop("TCGA-CDR contains duplicated project/patient IDs.")
}

# ------------------------------------------------------------
# Exact patient + project join
# ------------------------------------------------------------

m <- expr %>%
  left_join(
    cdr,
    by = c(
      "project" = "project",
      "cases.submitter_id" = "bcr_patient_barcode"
    )
  )

cat("\nMatched expression patients:",
    sum(!is.na(m$type)),
    "/",
    nrow(m),
    "\n")

cat("Unmatched expression patients:",
    sum(is.na(m$type)),
    "\n")

cat("Matched but redacted:",
    sum(m$Redacted %in% TRUE, na.rm = TRUE),
    "\n")

# ------------------------------------------------------------
# Endpoint readiness
# A survival record requires:
#   event indicator = 0/1
#   positive/non-negative time
# ------------------------------------------------------------

m <- m %>%
  mutate(
    OS_ready =
      !is.na(OS) &
      !is.na(OS.time) &
      OS.time >= 0,

    DSS_ready =
      !is.na(DSS) &
      !is.na(DSS.time) &
      DSS.time >= 0,

    DFI_ready =
      !is.na(DFI) &
      !is.na(DFI.time) &
      DFI.time >= 0,

    PFI_ready =
      !is.na(PFI) &
      !is.na(PFI.time) &
      PFI.time >= 0,

    survival_eligible =
      !Redacted %in% TRUE
  )

# ------------------------------------------------------------
# Project-level audit
# ------------------------------------------------------------

audit <- m %>%
  group_by(project) %>%
  summarise(
    Expression_patients = n(),

    CDR_matched =
      sum(!is.na(type)),

    Redacted =
      sum(Redacted %in% TRUE, na.rm = TRUE),

    OS_ready =
      sum(OS_ready & survival_eligible, na.rm = TRUE),

    OS_events =
      sum(OS == 1 &
          OS_ready &
          survival_eligible,
          na.rm = TRUE),

    PFI_ready =
      sum(PFI_ready & survival_eligible, na.rm = TRUE),

    PFI_events =
      sum(PFI == 1 &
          PFI_ready &
          survival_eligible,
          na.rm = TRUE),

    DFI_ready =
      sum(DFI_ready & survival_eligible, na.rm = TRUE),

    DFI_events =
      sum(DFI == 1 &
          DFI_ready &
          survival_eligible,
          na.rm = TRUE),

    DSS_ready =
      sum(DSS_ready & survival_eligible, na.rm = TRUE),

    DSS_events =
      sum(DSS == 1 &
          DSS_ready &
          survival_eligible,
          na.rm = TRUE),

    .groups = "drop"
  ) %>%
  arrange(project)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  m,
  "TCGA_JAK3_pan_cancer/results/JAK3_TCGA_CDR_patient_data.csv",
  row.names = FALSE
)

write.csv(
  audit,
  "TCGA_JAK3_pan_cancer/results/JAK3_TCGA_CDR_endpoint_availability.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n===== ENDPOINT AVAILABILITY =====\n")
print(
  as.data.frame(audit),
  row.names = FALSE
)

cat("\n===== TOTALS =====\n")

cat(
  "Projects with >=10 PFI events:",
  sum(audit$PFI_events >= 10),
  "\n"
)

cat(
  "Projects with >=10 OS events:",
  sum(audit$OS_events >= 10),
  "\n"
)

cat(
  "Projects with >=10 DFI events:",
  sum(audit$DFI_events >= 10),
  "\n"
)

cat(
  "Projects with >=10 DSS events:",
  sum(audit$DSS_events >= 10),
  "\n"
)

cat("\n51b TCGA-CDR MATCH AUDIT COMPLETE\n")
