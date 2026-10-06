
library(data.table)

# ==================================================
# 03 — MULTIKINASE CLINICAL METADATA INTEGRATION
# ==================================================

root <- "TCGA_multikinase_analysis"

outdir <- file.path(root, "results/expression")
auditdir <- file.path(root, "audits")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(auditdir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------
# 1. Load expression and original clinical data
# --------------------------------------------------

expr <- fread(
  file.path(outdir, "01_multikinase_TPM_complete.csv")
)

original <- fread(
  paste0(
    "TCGA_JAK3_pan_cancer/results/expression/",
    "JAK3_expression_clinical_complete_33projects.csv"
  )
)

# --------------------------------------------------
# 2. Validate the original dataset
# --------------------------------------------------

stopifnot(
  nrow(expr) == 115050,
  uniqueN(expr$file_id) == 11505,
  uniqueN(expr$gene) == 10,
  nrow(original) == 11505,
  !anyDuplicated(original$file_id),
  !anyDuplicated(expr[, .(file_id, gene)])
)

# --------------------------------------------------
# 3. Retain original metadata
# --------------------------------------------------

# Remove JAK3-specific expression columns.
# All other fields come from the completed study.

metadata <- copy(original)

metadata[, c(
  "gene",
  "counts",
  "TPM"
) := NULL]

# Keep GDC file_id as the matching identifier.
# Check complete correspondence before merging.

stopifnot(
  setequal(expr$file_id, metadata$file_id)
)

# --------------------------------------------------
# 4. Merge all ten genes
# --------------------------------------------------

combined <- merge(
  expr,
  metadata,
  by = "file_id",
  all.x = TRUE,
  sort = FALSE
)

# --------------------------------------------------
# 5. Validate merged data
# --------------------------------------------------

stopifnot(
  nrow(combined) == 115050,
  !anyDuplicated(combined[, .(file_id, gene)]),
  !anyNA(combined$file_id),
  !anyNA(combined$project),
  !anyNA(combined$sample_type),
  !anyNA(combined$cases.submitter_id)
)

combined[, patient_id := cases.submitter_id]

combined[, log2TPM1 := log2(TPM + 1)]

# --------------------------------------------------
# 6. Verify original JAK3 metadata and expression
# --------------------------------------------------

reference <- original[, .(
  file_id,
  reference_TPM = TPM,
  reference_project = project,
  reference_sample_type = sample_type,
  reference_patient = cases.submitter_id
)]

validation <- merge(
  combined[gene == "JAK3", .(
    file_id,
    TPM,
    project,
    sample_type,
    patient_id
  )],
  reference,
  by = "file_id"
)

stopifnot(
  nrow(validation) == 11505,
  all(validation$TPM == validation$reference_TPM),
  all(validation$project == validation$reference_project),
  all(validation$sample_type ==
      validation$reference_sample_type),
  all(validation$patient_id ==
      validation$reference_patient)
)

# --------------------------------------------------
# 7. Save complete dataset
# --------------------------------------------------

fwrite(
  combined,
  file.path(
    outdir,
    "03_multikinase_expression_clinical_complete.csv"
  )
)

# --------------------------------------------------
# 8. Cancer-specific availability audit
# --------------------------------------------------

availability <- combined[, .(
  Files = uniqueN(file_id),
  Patients = uniqueN(patient_id),
  Missing_TPM = sum(is.na(TPM))
), by = .(gene, project, sample_type)]

setorder(availability, gene, project, sample_type)

fwrite(
  availability,
  file.path(
    auditdir,
    "03_multikinase_sample_availability.csv"
  )
)

# --------------------------------------------------
# 9. Final report
# --------------------------------------------------

cat("\n========================================\n")
cat("MULTIKINASE CLINICAL INTEGRATION\n")
cat("========================================\n")

cat("Total rows:", nrow(combined), "\n")
cat("Genes:", uniqueN(combined$gene), "\n")
cat("Cancer projects:", uniqueN(combined$project), "\n")
cat("Unique GDC files:", uniqueN(combined$file_id), "\n")
cat("Missing TPM:", sum(is.na(combined$TPM)), "\n")

cat("\nJAK3 reference validation: PASSED\n")

cat("\n===== SAMPLE TYPE SUMMARY =====\n")

print(
  unique(combined[, .(
    file_id,
    project,
    sample_type,
    patient_id
  )])[, .(
    Files = .N,
    Patients = uniqueN(patient_id)
  ), by = sample_type]
)

cat("\n03 COMPLETE\n")

