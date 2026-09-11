# ============================================================
# 13g — TCGA-KIRC comparison of priority genes
#        JAK3 + EPHB2
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

jak3_file <-
  "results/13_TCGA_JAK3/13c_KIRC_JAK3_expression.csv"

ephb2_file <-
  "results/13_TCGA_EPHB2/13f_KIRC_EPHB2_expression.csv"

outdir <- "results/13_TCGA_priority_genes"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(jak3_file)) stop("Missing JAK3 file")
if (!file.exists(ephb2_file)) stop("Missing EPHB2 file")

jak3 <- read.csv(
  jak3_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ephb2 <- read.csv(
  ephb2_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

jak3$Gene <- "JAK3"
ephb2$Gene <- "EPHB2"

# Keep only successful extractions
jak3 <- jak3[jak3$extraction_status == "OK", ]
ephb2 <- ephb2[ephb2$extraction_status == "OK", ]

dat <- rbind(jak3, ephb2)

dat$log2TPM1 <- log2(dat$TPM + 1)

write.csv(
  dat,
  file.path(
    outdir,
    "13g_KIRC_JAK3_EPHB2_expression.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Unpaired tumor vs normal analysis
# ------------------------------------------------------------

genes <- c("JAK3", "EPHB2")

summary_list <- list()

for (g in genes) {

  z <- dat[dat$Gene == g, ]

  tumor <- z$TPM[
    z$sample_type == "Primary Tumor"
  ]

  normal <- z$TPM[
    z$sample_type == "Solid Tissue Normal"
  ]

  tumor_log <- log2(tumor + 1)
  normal_log <- log2(normal + 1)

  wt <- wilcox.test(
    tumor_log,
    normal_log,
    exact = FALSE
  )

  normal_p95 <- as.numeric(
    quantile(
      normal,
      probs = 0.95,
      na.rm = TRUE
    )
  )

  n_high <- sum(
    tumor > normal_p95,
    na.rm = TRUE
  )

  summary_list[[g]] <- data.frame(
    Gene = g,
    Tumor_n = length(tumor),
    Normal_n = length(normal),

    Tumor_median_TPM =
      median(tumor, na.rm = TRUE),

    Normal_median_TPM =
      median(normal, na.rm = TRUE),

    Tumor_mean_TPM =
      mean(tumor, na.rm = TRUE),

    Normal_mean_TPM =
      mean(normal, na.rm = TRUE),

    Median_log2TPM1_difference =
      median(tumor_log, na.rm = TRUE) -
      median(normal_log, na.rm = TRUE),

    Unpaired_Wilcoxon_P =
      wt$p.value,

    Normal_95th_percentile_TPM =
      normal_p95,

    Tumors_above_normal_P95 =
      n_high,

    Percent_tumors_above_normal_P95 =
      100 * n_high / length(tumor),

    stringsAsFactors = FALSE
  )
}

summary_tab <- do.call(
  rbind,
  summary_list
)

summary_tab$Unpaired_FDR <- p.adjust(
  summary_tab$Unpaired_Wilcoxon_P,
  method = "BH"
)

# ------------------------------------------------------------
# Matched patient analysis
# ------------------------------------------------------------

pair_results <- list()
paired_summary <- list()

for (g in genes) {

  z <- dat[dat$Gene == g, ]

  # patient-level median in case duplicate aliquots exist
  agg <- aggregate(
    TPM ~ patient_id + sample_type,
    data = z,
    FUN = median
  )

  tumor <- agg[
    agg$sample_type == "Primary Tumor",
    c("patient_id", "TPM")
  ]

  normal <- agg[
    agg$sample_type == "Solid Tissue Normal",
    c("patient_id", "TPM")
  ]

  names(tumor)[2] <- "Tumor_TPM"
  names(normal)[2] <- "Normal_TPM"

  pairs <- merge(
    tumor,
    normal,
    by = "patient_id"
  )

  pairs$Gene <- g

  pairs$Tumor_log2TPM1 <-
    log2(pairs$Tumor_TPM + 1)

  pairs$Normal_log2TPM1 <-
    log2(pairs$Normal_TPM + 1)

  pairs$Difference <-
    pairs$Tumor_log2TPM1 -
    pairs$Normal_log2TPM1

  pair_results[[g]] <- pairs

  if (nrow(pairs) >= 2) {

    wt <- wilcox.test(
      pairs$Tumor_log2TPM1,
      pairs$Normal_log2TPM1,
      paired = TRUE,
      exact = FALSE
    )

    p <- wt$p.value

  } else {

    p <- NA_real_
  }

  paired_summary[[g]] <- data.frame(
    Gene = g,

    Matched_patients =
      nrow(pairs),

    Median_paired_log2_difference =
      if (nrow(pairs) > 0)
        median(pairs$Difference)
      else NA_real_,

    Tumor_higher_n =
      if (nrow(pairs) > 0)
        sum(pairs$Difference > 0)
      else NA_integer_,

    Tumor_lower_n =
      if (nrow(pairs) > 0)
        sum(pairs$Difference < 0)
      else NA_integer_,

    Paired_Wilcoxon_P = p,

    stringsAsFactors = FALSE
  )
}

pairs_all <- do.call(
  rbind,
  pair_results
)

paired_tab <- do.call(
  rbind,
  paired_summary
)

paired_tab$Paired_FDR <- p.adjust(
  paired_tab$Paired_Wilcoxon_P,
  method = "BH"
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  summary_tab,
  file.path(
    outdir,
    "13g_KIRC_unpaired_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  paired_tab,
  file.path(
    outdir,
    "13g_KIRC_matched_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  pairs_all,
  file.path(
    outdir,
    "13g_KIRC_matched_patient_values.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("KIRC JAK3 + EPHB2 ANALYSIS\n")
cat("========================================\n\n")

cat("UNPAIRED TUMOR VS NORMAL\n\n")
print(
  summary_tab,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("MATCHED PATIENT ANALYSIS\n")
cat("========================================\n\n")

print(
  paired_tab,
  row.names = FALSE,
  digits = 5
)

cat("\nFiles saved under:\n")
cat("  results/13_TCGA_priority_genes/\n")

cat("\n========================================\n")
cat("13g COMPLETE\n")
cat("========================================\n")
