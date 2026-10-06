
library(data.table)

root <- "TCGA_multikinase_analysis"
outdir <- file.path(root, "audits")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Original JAK3 extraction
old <- fread(
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_TPM_complete.csv"
)

# New multikinase extraction
new <- fread(
  file.path(
    root,
    "results/expression/01_multikinase_TPM_complete.csv"
  )
)

# Original GDC metadata
meta <- fread(
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_GDC_metadata.csv"
)

# Retain JAK3
new <- new[gene == "JAK3"]

# Match original filenames to GDC file IDs
mapping <- unique(
  meta[, .(file_id, file_name)]
)

stopifnot(
  !anyDuplicated(mapping$file_id),
  !anyDuplicated(mapping$file_name),
  !anyDuplicated(old$sample),
  !anyDuplicated(new$file_id)
)

old <- merge(
  old,
  mapping,
  by.x = "sample",
  by.y = "file_name",
  all.x = TRUE
)

# Match by exact GDC file ID
comparison <- merge(
  old[, .(
    file_id,
    old_counts = counts,
    old_TPM = TPM
  )],
  new[, .(
    file_id,
    new_counts = counts,
    new_TPM = TPM
  )],
  by = "file_id",
  all = TRUE
)

comparison[, Count_difference :=
  new_counts - old_counts
]

comparison[, TPM_difference :=
  new_TPM - old_TPM
]

# Audit
cat("\n================================\n")
cat("JAK3 COMPLETE VALIDATION\n")
cat("================================\n")

cat("Original rows:", nrow(old), "\n")
cat("New JAK3 rows:", nrow(new), "\n")

cat(
  "Original files without metadata:",
  sum(is.na(old$file_id)),
  "\n"
)

cat(
  "Fully matched samples:",
  sum(
    complete.cases(
      comparison[, .(
        old_counts,
        new_counts,
        old_TPM,
        new_TPM
      )]
    )
  ),
  "\n"
)

cat(
  "Count mismatches:",
  sum(
    comparison$Count_difference != 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Maximum absolute TPM difference:",
  max(
    abs(comparison$TPM_difference),
    na.rm = TRUE
  ),
  "\n"
)

fwrite(
  comparison,
  file.path(outdir, "02_JAK3_validation.csv")
)

# Strict validation
stopifnot(
  nrow(comparison) == 11505,
  !anyNA(comparison$old_counts),
  !anyNA(comparison$new_counts),
  !anyNA(comparison$old_TPM),
  !anyNA(comparison$new_TPM),
  all(comparison$Count_difference == 0),
  all(abs(comparison$TPM_difference) < 1e-8)
)

cat("\nJAK3 VALIDATION PASSED\n")

