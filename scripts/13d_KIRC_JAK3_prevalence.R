# ============================================================
# 13d — TCGA-KIRC JAK3 TUMOR/NORMAL + PREVALENCE
# ============================================================

rm(list = ls())

infile <-
  "results/13_TCGA_JAK3/13c_KIRC_JAK3_expression.csv"

outdir <-
  "results/13_TCGA_JAK3"

if (!file.exists(infile)) {
  stop("Missing 13c KIRC JAK3 expression file.")
}

dat <- read.csv(
  infile,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dat <- dat[
  dat$extraction_status == "OK" &
  is.finite(dat$TPM),
]

cat("\n========================================\n")
cat("13d — KIRC JAK3 PREVALENCE\n")
cat("========================================\n\n")

cat("Raw records:", nrow(dat), "\n")

# ============================================================
# 1. PATIENT-LEVEL COLLAPSE
#
# If multiple aliquots exist for the same patient/sample type,
# use median TPM so a patient is not counted twice.
# ============================================================

patient_dat <- aggregate(
  TPM ~ patient_id + sample_type,
  data = dat,
  FUN = median
)

cat(
  "Patient-level records:",
  nrow(patient_dat),
  "\n\n"
)

cat("Patient sample types:\n")
print(
  table(patient_dat$sample_type)
)

tumor <- patient_dat$TPM[
  patient_dat$sample_type == "Primary Tumor"
]

normal <- patient_dat$TPM[
  patient_dat$sample_type == "Solid Tissue Normal"
]

# ============================================================
# 2. DESCRIPTIVE STATISTICS
# ============================================================

tumor_median <- median(tumor)
normal_median <- median(normal)

tumor_mean <- mean(tumor)
normal_mean <- mean(normal)

median_ratio <-
  tumor_median / normal_median

log2_median_ratio <-
  log2(median_ratio)

# ============================================================
# 3. UNPAIRED TUMOR vs NORMAL TEST
# ============================================================

wilcox_all <- wilcox.test(
  tumor,
  normal,
  alternative = "two.sided",
  exact = FALSE
)

# ============================================================
# 4. JAK3-HIGH PREVALENCE
#
# Primary threshold:
# 95th percentile of normal JAK3 TPM
# ============================================================

normal_q95 <- as.numeric(
  quantile(
    normal,
    probs = 0.95,
    na.rm = TRUE
  )
)

n_tumor_high <- sum(
  tumor > normal_q95
)

pct_tumor_high <-
  100 * n_tumor_high / length(tumor)

# Additional descriptive threshold:
# >= 2 x normal median

twofold_threshold <-
  2 * normal_median

n_tumor_2x <- sum(
  tumor >= twofold_threshold
)

pct_tumor_2x <-
  100 * n_tumor_2x / length(tumor)

# ============================================================
# 5. MATCHED PATIENT ANALYSIS
# ============================================================

tumor_df <- patient_dat[
  patient_dat$sample_type == "Primary Tumor",
  c("patient_id", "TPM")
]

normal_df <- patient_dat[
  patient_dat$sample_type == "Solid Tissue Normal",
  c("patient_id", "TPM")
]

colnames(tumor_df)[2] <- "Tumor_TPM"
colnames(normal_df)[2] <- "Normal_TPM"

paired <- merge(
  tumor_df,
  normal_df,
  by = "patient_id"
)

paired$log2FC_Tumor_vs_Normal <-
  log2(
    (paired$Tumor_TPM + 1e-6) /
    (paired$Normal_TPM + 1e-6)
  )

paired$Tumor_higher <-
  paired$Tumor_TPM >
  paired$Normal_TPM

paired$Tumor_at_least_2fold <-
  paired$Tumor_TPM >=
  2 * paired$Normal_TPM

if (nrow(paired) > 0) {

  paired_test <- wilcox.test(
    paired$Tumor_TPM,
    paired$Normal_TPM,
    paired = TRUE,
    alternative = "two.sided",
    exact = FALSE
  )

  paired_p <- paired_test$p.value

  pct_paired_higher <-
    100 * mean(
      paired$Tumor_higher
    )

  pct_paired_2fold <-
    100 * mean(
      paired$Tumor_at_least_2fold
    )

} else {

  paired_p <- NA_real_
  pct_paired_higher <- NA_real_
  pct_paired_2fold <- NA_real_
}

# ============================================================
# 6. SAVE SUMMARY
# ============================================================

summary_table <- data.frame(

  Cancer = "TCGA-KIRC",

  Tumor_N =
    length(tumor),

  Normal_N =
    length(normal),

  Tumor_median_TPM =
    tumor_median,

  Normal_median_TPM =
    normal_median,

  Tumor_mean_TPM =
    tumor_mean,

  Normal_mean_TPM =
    normal_mean,

  Median_ratio_Tumor_Normal =
    median_ratio,

  Log2_median_ratio =
    log2_median_ratio,

  Tumor_vs_Normal_Wilcox_P =
    wilcox_all$p.value,

  Normal_95th_percentile_TPM =
    normal_q95,

  JAK3_high_tumors =
    n_tumor_high,

  Percent_JAK3_high =
    pct_tumor_high,

  Twofold_normal_median_threshold =
    twofold_threshold,

  Tumors_ge_2x_normal_median =
    n_tumor_2x,

  Percent_ge_2x_normal_median =
    pct_tumor_2x,

  Matched_pairs =
    nrow(paired),

  Matched_Wilcox_P =
    paired_p,

  Percent_matched_tumor_higher =
    pct_paired_higher,

  Percent_matched_ge_2fold =
    pct_paired_2fold,

  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(
    outdir,
    "13d_KIRC_JAK3_prevalence_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  patient_dat,
  file.path(
    outdir,
    "13d_KIRC_JAK3_patient_level.csv"
  ),
  row.names = FALSE
)

write.csv(
  paired,
  file.path(
    outdir,
    "13d_KIRC_JAK3_matched_pairs.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 7. CONSOLE REPORT
# ============================================================

cat("\n========================================\n")
cat("TUMOR vs NORMAL\n")
cat("========================================\n\n")

cat("Tumor N:", length(tumor), "\n")
cat("Normal N:", length(normal), "\n\n")

cat(
  "Tumor median TPM:",
  round(tumor_median, 4),
  "\n"
)

cat(
  "Normal median TPM:",
  round(normal_median, 4),
  "\n"
)

cat(
  "Median tumor/normal ratio:",
  round(median_ratio, 3),
  "\n"
)

cat(
  "log2 median ratio:",
  round(log2_median_ratio, 3),
  "\n"
)

cat(
  "Wilcoxon P:",
  format(
    wilcox_all$p.value,
    scientific = TRUE
  ),
  "\n"
)

cat("\n========================================\n")
cat("JAK3-HIGH PREVALENCE\n")
cat("========================================\n\n")

cat(
  "Normal 95th percentile TPM:",
  round(normal_q95, 4),
  "\n"
)

cat(
  "JAK3-high tumors:",
  n_tumor_high,
  "/",
  length(tumor),
  "\n"
)

cat(
  "Prevalence:",
  round(pct_tumor_high, 2),
  "%\n"
)

cat(
  "Tumors >=2x normal median:",
  n_tumor_2x,
  "/",
  length(tumor),
  "=",
  round(pct_tumor_2x, 2),
  "%\n"
)

cat("\n========================================\n")
cat("MATCHED PATIENT ANALYSIS\n")
cat("========================================\n\n")

cat(
  "Matched tumor-normal patients:",
  nrow(paired),
  "\n"
)

cat(
  "Paired Wilcoxon P:",
  format(
    paired_p,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Matched patients with higher tumor JAK3:",
  round(pct_paired_higher, 2),
  "%\n"
)

cat(
  "Matched patients with >=2-fold tumor increase:",
  round(pct_paired_2fold, 2),
  "%\n"
)

cat("\n========================================\n")
cat("13d COMPLETE\n")
cat("========================================\n")
