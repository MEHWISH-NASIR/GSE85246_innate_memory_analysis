
library(data.table)
setDTthreads(1)

root <- "TCGA_multikinase_analysis"
auditdir <- file.path(root, "audits")
dir.create(auditdir, recursive = TRUE, showWarnings = FALSE)

new_tn <- fread(
  file.path(root, "results/primary/04_multikinase_tumor_normal.csv")
)[gene == "JAK3"]

new_mp <- fread(
  file.path(root, "results/primary/04_multikinase_matched_tumor_normal.csv")
)[gene == "JAK3"]

old_tn <- fread(
  "TCGA_JAK3_pan_cancer/results/JAK3_complete_tumor_normal.csv"
)

old_mp <- fread(
  "TCGA_JAK3_pan_cancer/results/JAK3_complete_matched_tumor_normal.csv"
)

v_tn <- merge(
  old_tn,
  new_tn,
  by = "project",
  suffixes = c("_old", "_new")
)

v_mp <- merge(
  old_mp,
  new_mp,
  by = "project",
  suffixes = c("_old", "_new")
)

# --------------------------------------------------
# Primary tumor-normal validation
# --------------------------------------------------

v_tn[, diff_tumor_n := tumor_n_new - tumor_n_old]
v_tn[, diff_normal_n := normal_n_new - normal_n_old]

v_tn[, diff_tumor_median :=
       tumor_median_TPM_new - tumor_median_TPM_old]

v_tn[, diff_normal_median :=
       normal_median_TPM_new - normal_median_TPM_old]

v_tn[, diff_effect :=
       Median_log2TPM1_difference_new -
       Median_log2TPM1_difference_old]

v_tn[, diff_P :=
       Wilcoxon_P_new - Wilcoxon_P_old]

# Stored FDR comparison retained only as an audit
v_tn[, stored_FDR_difference :=
       Wilcoxon_FDR_new - Wilcoxon_FDR_old]

# --------------------------------------------------
# Matched validation
# --------------------------------------------------

v_mp[, diff_pairs :=
       matched_pairs_new - matched_pairs_old]

v_mp[, diff_effect :=
       median_log2_difference_new -
       median_log2_difference_old]

v_mp[, diff_P :=
       Paired_Wilcoxon_P_new -
       Paired_Wilcoxon_P_old]

v_mp[, diff_FDR :=
       Paired_Wilcoxon_FDR_new -
       Paired_Wilcoxon_FDR_old]

# --------------------------------------------------
# Save audits
# --------------------------------------------------

fwrite(
  v_tn,
  file.path(
    auditdir,
    "04b_JAK3_tumor_normal_validation.csv"
  )
)

fwrite(
  v_mp,
  file.path(
    auditdir,
    "04b_JAK3_matched_validation.csv"
  )
)

# --------------------------------------------------
# Report
# --------------------------------------------------

cat("\n========================================\n")
cat("JAK3 PRIMARY ANALYSIS VALIDATION\n")
cat("========================================\n")

cat("\nTumor-normal cancers:", nrow(v_tn), "\n")

cat("Max tumor-N difference:",
    max(abs(v_tn$diff_tumor_n)), "\n")

cat("Max normal-N difference:",
    max(abs(v_tn$diff_normal_n)), "\n")

cat("Max tumor median difference:",
    max(abs(v_tn$diff_tumor_median)), "\n")

cat("Max normal median difference:",
    max(abs(v_tn$diff_normal_median)), "\n")

cat("Max effect difference:",
    max(abs(v_tn$diff_effect)), "\n")

cat("Max raw P difference:",
    max(abs(v_tn$diff_P), na.rm = TRUE), "\n")

cat("Historical stored FDR difference:",
    max(abs(v_tn$stored_FDR_difference),
        na.rm = TRUE), "\n")

cat("\nMatched cancers:", nrow(v_mp), "\n")

cat("Max pair difference:",
    max(abs(v_mp$diff_pairs)), "\n")

cat("Max matched effect difference:",
    max(abs(v_mp$diff_effect)), "\n")

cat("Max matched P difference:",
    max(abs(v_mp$diff_P), na.rm = TRUE), "\n")

cat("Max matched FDR difference:",
    max(abs(v_mp$diff_FDR), na.rm = TRUE), "\n")

# --------------------------------------------------
# Strict validation of underlying analysis
# --------------------------------------------------

stopifnot(
  nrow(v_tn) == 24,
  nrow(v_mp) == 22,

  all(v_tn$diff_tumor_n == 0),
  all(v_tn$diff_normal_n == 0),

  max(abs(v_tn$diff_tumor_median)) < 1e-12,
  max(abs(v_tn$diff_normal_median)) < 1e-12,
  max(abs(v_tn$diff_effect)) < 1e-12,
  max(abs(v_tn$diff_P), na.rm = TRUE) < 1e-12,

  all(v_mp$diff_pairs == 0),
  max(abs(v_mp$diff_effect)) < 1e-12,
  max(abs(v_mp$diff_P), na.rm = TRUE) < 1e-12,
  max(abs(v_mp$diff_FDR), na.rm = TRUE) < 1e-12
)

cat("\nJAK3 UNDERLYING ANALYSIS REPRODUCED EXACTLY\n")
cat("Historical tumor-normal FDR difference retained as audit only.\n")
cat("04b COMPLETE\n")

