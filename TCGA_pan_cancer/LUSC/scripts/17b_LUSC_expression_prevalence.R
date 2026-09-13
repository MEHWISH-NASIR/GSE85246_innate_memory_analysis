# ============================================================
# 17b — TCGA-LUSC
# JAK3 + EPHB2 tumor-normal, matched pairs, prevalence
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

infile <-
  "results/17_TCGA_LUSC/17a_LUSC_JAK3_EPHB2_expression.csv"

outdir <- "results/17_TCGA_LUSC"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

x <- read.csv(
  infile,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

x <- x[
  x$extraction_status == "OK" &
  !is.na(x$TPM),
]

x$log2TPM1 <- log2(x$TPM + 1)

# ------------------------------------------------------------
# Patient-level expression
# Collapse duplicate aliquots using median
# ------------------------------------------------------------

patient <- aggregate(
  TPM ~ gene_name + patient_id + sample_type,
  data = x,
  FUN = median
)

patient$log2TPM1 <-
  log2(patient$TPM + 1)

genes <- c("JAK3", "EPHB2")

unpaired_list <- list()
prevalence_list <- list()
matched_list <- list()
matched_values <- list()

for (g in genes) {

  # ==========================================================
  # Sample-level tumor vs normal
  # ==========================================================

  z <- x[x$gene_name == g, ]

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

  unpaired_list[[g]] <- data.frame(
    Gene = g,

    Tumor_n =
      length(tumor),

    Normal_n =
      length(normal),

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

    Wilcoxon_P =
      wt$p.value,

    stringsAsFactors = FALSE
  )

  # ==========================================================
  # Patient-level prevalence
  # ==========================================================

  zp <- patient[
    patient$gene_name == g,
  ]

  tp <- zp$TPM[
    zp$sample_type == "Primary Tumor"
  ]

  np <- zp$TPM[
    zp$sample_type == "Solid Tissue Normal"
  ]

  tumor_med <-
    median(tp, na.rm = TRUE)

  normal_med <-
    median(np, na.rm = TRUE)

  direction <- if (
    tumor_med > normal_med
  ) {
    "UP"
  } else {
    "DOWN"
  }

  p05 <- as.numeric(
    quantile(
      np,
      0.05,
      na.rm = TRUE
    )
  )

  p95 <- as.numeric(
    quantile(
      np,
      0.95,
      na.rm = TRUE
    )
  )

  if (direction == "UP") {

    threshold <- p95

    aberrant_n <- sum(
      tp > threshold,
      na.rm = TRUE
    )

    definition <-
      "Tumor TPM > normal 95th percentile"

  } else {

    threshold <- p05

    aberrant_n <- sum(
      tp < threshold,
      na.rm = TRUE
    )

    definition <-
      "Tumor TPM < normal 5th percentile"
  }

  prevalence_list[[g]] <- data.frame(
    Gene = g,

    Direction =
      direction,

    Tumor_patient_n =
      length(tp),

    Normal_patient_n =
      length(np),

    Tumor_median_TPM =
      tumor_med,

    Normal_median_TPM =
      normal_med,

    Normal_P05_TPM =
      p05,

    Normal_P95_TPM =
      p95,

    Prevalence_threshold_TPM =
      threshold,

    Prevalence_definition =
      definition,

    Aberrant_tumor_n =
      aberrant_n,

    Aberrant_tumor_percent =
      100 * aberrant_n /
      length(tp),

    stringsAsFactors = FALSE
  )

  # ==========================================================
  # Matched patient analysis
  # ==========================================================

  tum <- zp[
    zp$sample_type == "Primary Tumor",
    c("patient_id", "TPM")
  ]

  nor <- zp[
    zp$sample_type == "Solid Tissue Normal",
    c("patient_id", "TPM")
  ]

  names(tum)[2] <- "Tumor_TPM"
  names(nor)[2] <- "Normal_TPM"

  pairs <- merge(
    tum,
    nor,
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

  matched_values[[g]] <- pairs

  if (nrow(pairs) >= 2) {

    pwt <- wilcox.test(
      pairs$Tumor_log2TPM1,
      pairs$Normal_log2TPM1,
      paired = TRUE,
      exact = FALSE
    )

    pp <- pwt$p.value

  } else {

    pp <- NA_real_
  }

  matched_list[[g]] <- data.frame(
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

    Paired_Wilcoxon_P =
      pp,

    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Combine + FDR
# ------------------------------------------------------------

unpaired <- do.call(
  rbind,
  unpaired_list
)

unpaired$Wilcoxon_FDR <- p.adjust(
  unpaired$Wilcoxon_P,
  method = "BH"
)

prevalence <- do.call(
  rbind,
  prevalence_list
)

matched <- do.call(
  rbind,
  matched_list
)

matched$Paired_FDR <- p.adjust(
  matched$Paired_Wilcoxon_P,
  method = "BH"
)

pairs_all <- do.call(
  rbind,
  matched_values
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  unpaired,
  file.path(
    outdir,
    "17b_LUSC_tumor_normal_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  prevalence,
  file.path(
    outdir,
    "17b_LUSC_prevalence_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  matched,
  file.path(
    outdir,
    "17b_LUSC_matched_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  pairs_all,
  file.path(
    outdir,
    "17b_LUSC_matched_values.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("LUSC TUMOR VS NORMAL\n")
cat("========================================\n\n")

print(
  unpaired,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("LUSC DIRECTIONAL PREVALENCE\n")
cat("========================================\n\n")

print(
  prevalence,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("LUSC MATCHED PATIENTS\n")
cat("========================================\n\n")

print(
  matched,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("17b COMPLETE\n")
cat("========================================\n")
