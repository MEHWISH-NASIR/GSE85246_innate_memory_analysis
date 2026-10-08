
rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

setDTthreads(1)

cat("\n========================================\n")
cat("19a — GSE58310 JAK3 PAIRED VALIDATION\n")
cat("========================================\n\n")

infile <- "data/transcriptomics/GSE58310/GSE58310_GeneExpression.csv.gz"

outdir <- "results/19_JAK3_cross_dataset_validation/GSE58310"
figdir <- file.path(outdir, "figures")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(infile)) {
  stop("Missing input file: ", infile)
}

# ============================================================
# 1. LOAD MATRIX
# ============================================================

x <- fread(
  infile,
  nThread = 1
)

cat("Rows:", nrow(x), "\n")
cat("Columns:", ncol(x), "\n")

if (!"feature_id" %in% names(x)) {
  stop("feature_id column not found.")
}

# ============================================================
# 2. EXTRACT JAK3
# ============================================================

JAK3_ID <- "ENSG00000105639"

jak3 <- x[
  feature_id == JAK3_ID
]

if (nrow(jak3) != 1) {
  stop("Expected exactly one JAK3 row.")
}

cat("\nJAK3 row found:", JAK3_ID, "\n")

# ============================================================
# 3. DEFINE EXACT SAMPLE GROUPS
# ============================================================

donors <- c("BC8", "BC9", "BC11", "BC12")

rpmid6_cols <- c(
  "Monocytes_6d_RPMI_BC8_4640",
  "Monocytes_6d_RPMI_BC9_4818",
  "Monocytes_6d_RPMI_BC11_5384",
  "Monocytes_6d_RPMI_BC12_5388"
)

lpsd6_cols <- c(
  "Monocytes_6d_LPS_BC8_4641",
  "Monocytes_6d_LPS_BC9_4819",
  "Monocytes_6d_LPS_BC11_5385",
  "Monocytes_6d_LPS_BC12_5389"
)

bgd6_cols <- c(
  "Monocytes_6d_BG_BC8_4642",
  "Monocytes_6d_BG_BC9_4820",
  "Monocytes_6d_BG_BC11_5386",
  "Monocytes_6d_BG_BC12_5390"
)

day0_cols <- c(
  "Monocytes_0d_RPMI_BC8_4639",
  "Monocytes_0d_RPMI_BC9_4817",
  "Monocytes_0d_RPMI_BC11_5383",
  "Monocytes_0d_RPMI_BC12_5387"
)

needed_cols <- c(
  day0_cols,
  rpmid6_cols,
  lpsd6_cols,
  bgd6_cols
)

missing_cols <- setdiff(
  needed_cols,
  names(jak3)
)

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing expected columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# ============================================================
# 4. BUILD DONOR-LEVEL TABLE
# ============================================================

dat <- data.frame(
  donor = donors,
  Day0_RPMI = as.numeric(jak3[, ..day0_cols]),
  Day6_RPMI = as.numeric(jak3[, ..rpmid6_cols]),
  Day6_LPS  = as.numeric(jak3[, ..lpsd6_cols]),
  Day6_BG   = as.numeric(jak3[, ..bgd6_cols]),
  stringsAsFactors = FALSE
)

dat$LPS_minus_RPMI <- dat$Day6_LPS - dat$Day6_RPMI
dat$BG_minus_RPMI  <- dat$Day6_BG  - dat$Day6_RPMI

write.csv(
  dat,
  file.path(
    outdir,
    "19a_GSE58310_JAK3_donor_values.csv"
  ),
  row.names = FALSE
)

cat("\n===== DONOR-LEVEL VALUES =====\n")
print(dat, row.names = FALSE)

# ============================================================
# 5. PAIRED TEST FUNCTION
# ============================================================

paired_summary <- function(
  treated,
  control,
  contrast_name
) {

  diff <- treated - control
  n <- length(diff)

  mean_diff <- mean(diff)
  sd_diff <- sd(diff)
  se_diff <- sd_diff / sqrt(n)

  tcrit <- qt(
    0.975,
    df = n - 1
  )

  ci_low <- mean_diff - tcrit * se_diff
  ci_high <- mean_diff + tcrit * se_diff

  ttest <- t.test(
    treated,
    control,
    paired = TRUE
  )

  # exact = FALSE avoids exact-distribution warnings with ties/zeros
  wtest <- suppressWarnings(
    wilcox.test(
      treated,
      control,
      paired = TRUE,
      exact = FALSE
    )
  )

  n_positive <- sum(diff > 0)
  n_negative <- sum(diff < 0)
  n_zero <- sum(diff == 0)

  data.frame(
    dataset = "GSE58310",
    gene = "JAK3",
    gene_id = JAK3_ID,
    contrast = contrast_name,
    n_donors = n,
    mean_treated = mean(treated),
    mean_control = mean(control),
    mean_paired_effect_lnRPKM = mean_diff,
    sd_paired_effect = sd_diff,
    se_paired_effect = se_diff,
    CI95_low = ci_low,
    CI95_high = ci_high,
    paired_t_p = unname(ttest$p.value),
    wilcoxon_p = unname(wtest$p.value),
    positive_donors = n_positive,
    negative_donors = n_negative,
    zero_donors = n_zero,
    direction_consistency = max(
      n_positive,
      n_negative
    ) / n,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 6. PRIMARY COMPARISONS
# ============================================================

res_lps <- paired_summary(
  dat$Day6_LPS,
  dat$Day6_RPMI,
  "Day6_LPS_vs_Day6_RPMI"
)

res_bg <- paired_summary(
  dat$Day6_BG,
  dat$Day6_RPMI,
  "Day6_BG_vs_Day6_RPMI"
)

summary_results <- rbind(
  res_lps,
  res_bg
)

write.csv(
  summary_results,
  file.path(
    outdir,
    "19a_GSE58310_JAK3_paired_summary.csv"
  ),
  row.names = FALSE
)

cat("\n===== PAIRED SUMMARY =====\n")
print(
  summary_results,
  row.names = FALSE
)

# ============================================================
# 7. DESCRIPTIVE CLASSIFICATION
#
# Important:
# This is not a genome-wide DE call.
# It is targeted single-gene validation.
# ============================================================

classify <- function(
  mean_effect,
  ci_low,
  ci_high,
  positive_donors,
  negative_donors,
  n
) {

  if (
    mean_effect > 0 &&
    ci_low > 0
  ) {
    return("Consistent positive post-washout JAK3 shift")
  }

  if (
    mean_effect < 0 &&
    ci_high < 0
  ) {
    return("Consistent negative post-washout JAK3 shift")
  }

  if (
    positive_donors == n
  ) {
    return("All donors positive; uncertainty remains")
  }

  if (
    negative_donors == n
  ) {
    return("All donors negative; uncertainty remains")
  }

  return("Heterogeneous / uncertain")
}

summary_results$classification <- mapply(
  classify,
  mean_effect =
    summary_results$mean_paired_effect_lnRPKM,
  ci_low =
    summary_results$CI95_low,
  ci_high =
    summary_results$CI95_high,
  positive_donors =
    summary_results$positive_donors,
  negative_donors =
    summary_results$negative_donors,
  n =
    summary_results$n_donors
)

write.csv(
  summary_results,
  file.path(
    outdir,
    "19a_GSE58310_JAK3_final_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 8. LONG FORMAT FOR FIGURE
# ============================================================

plotdat <- rbind(
  data.frame(
    donor = dat$donor,
    condition = "Day6 RPMI",
    expression = dat$Day6_RPMI
  ),
  data.frame(
    donor = dat$donor,
    condition = "Day6 LPS",
    expression = dat$Day6_LPS
  ),
  data.frame(
    donor = dat$donor,
    condition = "Day6 beta-glucan",
    expression = dat$Day6_BG
  )
)

plotdat$condition <- factor(
  plotdat$condition,
  levels = c(
    "Day6 RPMI",
    "Day6 LPS",
    "Day6 beta-glucan"
  )
)

# ============================================================
# 9. PAIRED DONOR FIGURE
# ============================================================

p <- ggplot(
  plotdat,
  aes(
    x = condition,
    y = expression,
    group = donor
  )
) +
  geom_line(
    alpha = 0.6
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    title = "GSE58310: JAK3 expression after stimulus washout",
    subtitle = "Paired Day6 measurements across four independent donors",
    x = NULL,
    y = "JAK3 expression (ln RPKM)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

ggsave(
  file.path(
    figdir,
    "19a_GSE58310_JAK3_paired_expression.png"
  ),
  p,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "19a_GSE58310_JAK3_paired_expression.pdf"
  ),
  p,
  width = 7,
  height = 5
)

# ============================================================
# 10. WRITE STUDY METADATA
# ============================================================

metadata <- data.frame(
  dataset = "GSE58310",
  species = "Homo sapiens",
  cell_system = "Primary monocytes differentiated toward macrophage state",
  stimulus_duration = "24 h",
  stimuli = "LPS 100 ng/mL; beta-glucan 5 ug/mL",
  washout = "Yes",
  late_timepoint = "Day 6",
  matched_control = "Day6 RPMI",
  biological_replicates = 4,
  assay = "RNA-seq",
  processed_scale = "natural-log RPKM",
  primary_question = "Does JAK3 remain deregulated at Day6 after stimulus removal?",
  stringsAsFactors = FALSE
)

write.csv(
  metadata,
  file.path(
    outdir,
    "19a_GSE58310_study_metadata.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("FINAL CLASSIFICATIONS\n")
cat("========================================\n")

print(
  summary_results[
    ,
    c(
      "contrast",
      "mean_paired_effect_lnRPKM",
      "CI95_low",
      "CI95_high",
      "paired_t_p",
      "wilcoxon_p",
      "positive_donors",
      "negative_donors",
      "classification"
    )
  ],
  row.names = FALSE
)

cat("\n19a COMPLETE\n")

