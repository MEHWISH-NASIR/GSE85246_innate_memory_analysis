library(dplyr)

tpm_file <- paste0(
  "TCGA_JAK3_pan_cancer/results/expression/",
  "JAK3_TPM_complete.csv"
)

meta_file <- paste0(
  "TCGA_JAK3_pan_cancer/results/expression/",
  "JAK3_GDC_metadata.csv"
)

out_file <- paste0(
  "TCGA_JAK3_pan_cancer/results/expression/",
  "JAK3_final_TPM_metadata_complete.csv"
)

tpm <- read.csv(
  tpm_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta <- read.csv(
  meta_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("JAK3 COMPLETE TPM + METADATA MERGE\n")
cat("========================================\n")

cat("\nTPM rows:", nrow(tpm), "\n")
cat("Unique TPM samples:", n_distinct(tpm$sample), "\n")

cat("\nMetadata rows:", nrow(meta), "\n")
cat("Unique metadata filenames:", n_distinct(meta$file_name), "\n")
cat("Metadata projects:", n_distinct(meta$project), "\n")

# ------------------------------------------------------------
# Mandatory uniqueness checks
# ------------------------------------------------------------

if (nrow(tpm) != 11505) {
  stop("Expected 11505 TPM rows; found ", nrow(tpm))
}

if (n_distinct(tpm$sample) != 11505) {
  stop("TPM sample filenames are not unique.")
}

if (nrow(meta) != 11505) {
  stop("Expected 11505 metadata rows; found ", nrow(meta))
}

if (n_distinct(meta$file_name) != 11505) {
  stop("GDC metadata file_name values are not unique.")
}

# ------------------------------------------------------------
# Pre-merge overlap audit
# ------------------------------------------------------------

tpm_not_meta <- setdiff(
  tpm$sample,
  meta$file_name
)

meta_not_tpm <- setdiff(
  meta$file_name,
  tpm$sample
)

cat("\nTPM files absent from metadata:",
    length(tpm_not_meta), "\n")

cat("Metadata files absent from TPM:",
    length(meta_not_tpm), "\n")

if (length(tpm_not_meta) > 0 ||
    length(meta_not_tpm) > 0) {
  stop("TPM and metadata file sets do not match exactly.")
}

# ------------------------------------------------------------
# Exact filename merge
# ------------------------------------------------------------

merged <- tpm %>%
  left_join(
    meta,
    by = c("sample" = "file_name")
  )

cat("\nMerged rows:", nrow(merged), "\n")
cat("Unique files:", n_distinct(merged$sample), "\n")
cat("Unique projects:", n_distinct(merged$project), "\n")

cat("\nMissing project:", sum(is.na(merged$project)), "\n")
cat("Missing file_id:", sum(is.na(merged$file_id)), "\n")
cat("Missing case ID:", sum(is.na(merged$cases.submitter_id)), "\n")
cat("Missing sample ID:", sum(is.na(merged$sample.submitter_id)), "\n")
cat("Missing sample_type:", sum(is.na(merged$sample_type)), "\n")
cat("Missing TPM:", sum(is.na(merged$TPM)), "\n")

if (nrow(merged) != 11505) {
  stop("Merged table does not contain 11505 rows.")
}

if (n_distinct(merged$project) != 33) {
  stop(
    "Expected 33 TCGA projects; found ",
    n_distinct(merged$project)
  )
}

# ------------------------------------------------------------
# Project-level summary
# ------------------------------------------------------------

project_summary <- merged %>%
  group_by(project) %>%
  summarise(
    expression_files = n(),
    unique_cases = n_distinct(cases.submitter_id),
    unique_samples = n_distinct(sample.submitter_id),
    .groups = "drop"
  ) %>%
  arrange(project)

cat("\n===== 33-PROJECT SUMMARY =====\n")
print(project_summary, row.names = FALSE)

# ------------------------------------------------------------
# Sample-type summary
# ------------------------------------------------------------

sample_summary <- merged %>%
  count(
    project,
    sample_type,
    name = "n"
  ) %>%
  arrange(project, desc(n))

cat("\n===== SAMPLE TYPES =====\n")
print(sample_summary, row.names = FALSE)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  merged,
  out_file,
  row.names = FALSE
)

write.csv(
  project_summary,
  paste0(
    "TCGA_JAK3_pan_cancer/results/",
    "JAK3_complete_project_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  sample_summary,
  paste0(
    "TCGA_JAK3_pan_cancer/results/",
    "JAK3_complete_sampletype_summary.csv"
  ),
  row.names = FALSE
)

cat("\nCOMPLETE METADATA MERGE FINISHED\n")
cat("Saved:\n", out_file, "\n")
