# ============================================================
# 60a — JAK3 matched ATAC-RNA patient audit
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

atac_file <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results/JAK3_linked_ATAC_tissue_values.csv"

rna_file <-
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_expression_clinical_complete_33projects.csv"

outdir <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results"

priority_cancers <- c(
  "CHOL","KIRC","KIRP","THCA","HNSC",
  "STAD","LUAD","COAD","LIHC"
)

cat("\n========================================\n")
cat("60a JAK3 ATAC-RNA MATCH AUDIT\n")
cat("========================================\n\n")

atac <- read.csv(
  atac_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rna <- read.csv(
  rna_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("ATAC rows:", nrow(atac), "\n")
cat("ATAC unique patients:", length(unique(atac$patient_id)), "\n")
cat("RNA rows:", nrow(rna), "\n\n")

# ------------------------------------------------------------
# Identify required RNA columns
# ------------------------------------------------------------

patient_candidates <- c(
  "cases.submitter_id",
  "patient_id",
  "case_submitter_id"
)

patient_col <- patient_candidates[
  patient_candidates %in% names(rna)
][1]

if (is.na(patient_col)) {
  stop("Could not identify RNA patient-ID column.")
}

if (!"project" %in% names(rna)) {
  stop("RNA project column not found.")
}

if (!"sample_type" %in% names(rna)) {
  stop("RNA sample_type column not found.")
}

cat("RNA patient column:", patient_col, "\n")

# ------------------------------------------------------------
# JAK3 expression field
# ------------------------------------------------------------

if ("log2TPM1" %in% names(rna)) {

  expression_col <- "log2TPM1"

} else if ("TPM" %in% names(rna)) {

  rna$log2TPM1 <- log2(rna$TPM + 1)
  expression_col <- "log2TPM1"

} else {

  stop("Neither log2TPM1 nor TPM was found in RNA table.")
}

cat("RNA expression column:", expression_col, "\n\n")

# ------------------------------------------------------------
# Restrict to nine priority cancers and primary tumors
# ------------------------------------------------------------

rna$Cancer <- sub(
  "^TCGA-",
  "",
  rna$project
)

rna$patient_match_id <-
  as.character(
    rna[[patient_col]]
  )

rna_primary <- rna[
  rna$Cancer %in% priority_cancers &
  rna$sample_type == "Primary Tumor",
]

# ------------------------------------------------------------
# Collapse RNA to patient level
# ------------------------------------------------------------

rna_patient <- aggregate(
  rna_primary[[expression_col]],
  by = list(
    Cancer = rna_primary$Cancer,
    patient_match_id = rna_primary$patient_match_id
  ),
  FUN = mean,
  na.rm = TRUE
)

names(rna_patient)[3] <- "JAK3_log2TPM1"

# ------------------------------------------------------------
# Unique ATAC patients
# ------------------------------------------------------------

atac_patients <- unique(
  atac[
    ,
    c(
      "Cancer",
      "patient_id"
    )
  ]
)

# ------------------------------------------------------------
# Match audit by cancer
# ------------------------------------------------------------

rows <- list()

for (ca in priority_cancers) {

  a <- atac_patients[
    atac_patients$Cancer == ca,
    ,
    drop = FALSE
  ]

  r <- rna_patient[
    rna_patient$Cancer == ca,
    ,
    drop = FALSE
  ]

  matched <- intersect(
    a$patient_id,
    r$patient_match_id
  )

  rows[[length(rows) + 1]] <- data.frame(

    Cancer = ca,

    ATAC_patients =
      length(unique(a$patient_id)),

    RNA_primary_patients =
      length(unique(r$patient_match_id)),

    Matched_ATAC_RNA_patients =
      length(matched),

    Match_percent =
      ifelse(
        nrow(a) > 0,
        100 * length(matched) /
          length(unique(a$patient_id)),
        NA_real_
      ),

    stringsAsFactors = FALSE
  )
}

audit <- do.call(
  rbind,
  rows
)

write.csv(
  audit,
  file.path(
    outdir,
    "JAK3_ATAC_RNA_match_audit.csv"
  ),
  row.names = FALSE
)

cat("========================================\n")
cat("MATCHED PATIENTS BY CANCER\n")
cat("========================================\n\n")

print(
  audit,
  row.names = FALSE
)

cat(
  "\nTotal matched ATAC-RNA patients:",
  sum(audit$Matched_ATAC_RNA_patients),
  "\n"
)

cat("\n========================================\n")
cat("60a COMPLETE\n")
cat("========================================\n")
