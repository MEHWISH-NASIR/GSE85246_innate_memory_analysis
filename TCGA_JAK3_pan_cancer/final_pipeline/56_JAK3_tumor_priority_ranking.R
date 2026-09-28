# ============================================================
# 56 — JAK3 TCGA TUMOR PRIORITY RANKING
#
# Priority criteria:
# 1. Strong tumor-vs-normal upregulation
# 2. Upregulation already present in Stage-I tumors
# 3. Persistence after measured immune-content adjustment
#
# The priority score is descriptive, not an inferential test.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(ggplot2))

indir  <- "TCGA_JAK3_pan_cancer/results/release"
outdir <- "TCGA_JAK3_pan_cancer/results/release"
figdir <- "TCGA_JAK3_pan_cancer/figures/final"

dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Read canonical final tables
# ------------------------------------------------------------

master <- read.csv(
  file.path(indir, "JAK3_master_pan_cancer_summary.csv"),
  check.names = FALSE
)

flags <- read.csv(
  file.path(indir, "JAK3_master_key_flags_corrected.csv"),
  check.names = FALSE
)

dat <- merge(
  master,
  flags,
  by = "project",
  all.x = TRUE
)

dat$Cancer <- sub("^TCGA-", "", dat$project)

# ------------------------------------------------------------
# Helper: relative strength among significant positive cancers
#
# Non-significant / non-positive result = 0
# Significant positive effects are rank-scaled from 0.5 to 1.
# ------------------------------------------------------------

strength_score <- function(effect, significant_up) {

  out <- rep(0, length(effect))

  idx <- which(
    significant_up %in% TRUE &
    !is.na(effect) &
    effect > 0
  )

  n <- length(idx)

  if (n == 1) {
    out[idx] <- 1
  }

  if (n > 1) {
    r <- rank(effect[idx], ties.method = "average")

    out[idx] <- 0.5 + 0.5 * (
      (r - 1) / (n - 1)
    )
  }

  out
}

# ------------------------------------------------------------
# Component A — tumor-vs-normal upregulation
# ------------------------------------------------------------

dat$Tumor_UP_strength <- strength_score(
  dat$TN_median_log2TPM1_difference,
  dat$Canonical_TN_significant_UP
)

# ------------------------------------------------------------
# Component B — Stage-I-vs-normal upregulation
# ------------------------------------------------------------

dat$Early_UP_strength <- strength_score(
  dat$StageI_median_log2TPM1_difference,
  dat$Canonical_StageI_significant_UP
)

# ------------------------------------------------------------
# Component C — immune robustness
#
# Four binary persistence checks:
# TN after PTPRC
# TN after 10-marker ImmuneScore
# Stage-I after PTPRC
# Stage-I after 10-marker ImmuneScore
# ------------------------------------------------------------

immune_mat <- cbind(
  as.numeric(dat$TN_positive_persists_after_PTPRC),
  as.numeric(dat$TN_positive_persists_after_ImmuneScore),
  as.numeric(dat$StageI_positive_persists_after_PTPRC),
  as.numeric(dat$StageI_positive_persists_after_ImmuneScore)
)

dat$Immune_robustness <- rowMeans(
  immune_mat,
  na.rm = TRUE
)

dat$Immune_robustness[
  rowSums(!is.na(immune_mat)) == 0
] <- NA_real_

# ------------------------------------------------------------
# Ranking eligibility
#
# Require tumor-normal, Stage-I, and immune-adjustment information.
# ------------------------------------------------------------

dat$Ranking_eligible <-
  !is.na(dat$TN_median_log2TPM1_difference) &
  !is.na(dat$StageI_median_log2TPM1_difference) &
  !is.na(dat$Immune_robustness)

# ------------------------------------------------------------
# Conservative core candidate
#
# Passes ALL requested criteria.
# ------------------------------------------------------------

dat$Core_candidate <-
  dat$Ranking_eligible &
  dat$Canonical_TN_significant_UP %in% TRUE &
  dat$Canonical_StageI_significant_UP %in% TRUE &
  dat$TN_positive_persists_after_PTPRC %in% TRUE &
  dat$TN_positive_persists_after_ImmuneScore %in% TRUE &
  dat$StageI_positive_persists_after_PTPRC %in% TRUE &
  dat$StageI_positive_persists_after_ImmuneScore %in% TRUE

# ------------------------------------------------------------
# Integrated score
#
# Equal weights:
# 1/3 tumor upregulation
# 1/3 early-stage upregulation
# 1/3 immune robustness
# ------------------------------------------------------------

dat$Priority_score <- NA_real_

i <- dat$Ranking_eligible

dat$Priority_score[i] <- 100 * (
  dat$Tumor_UP_strength[i] +
  dat$Early_UP_strength[i] +
  dat$Immune_robustness[i]
) / 3

# ------------------------------------------------------------
# Sort ranking
#
# Strict core candidates first, then integrated score.
# ------------------------------------------------------------

ranking <- dat[dat$Ranking_eligible, ]

ranking <- ranking[
  order(
    -as.numeric(ranking$Core_candidate),
    -ranking$Priority_score,
    -ranking$TN_median_log2TPM1_difference,
    -ranking$StageI_median_log2TPM1_difference
  ),
]

ranking$Priority_rank <- seq_len(nrow(ranking))

# ------------------------------------------------------------
# Final review table
# ------------------------------------------------------------

ranking_out <- ranking[
  ,
  c(
    "Priority_rank",
    "project",
    "Cancer",
    "Priority_score",
    "Core_candidate",

    "TN_median_log2TPM1_difference",
    "TN_FDR",
    "Canonical_TN_significant_UP",
    "Tumor_UP_strength",

    "StageI_median_log2TPM1_difference",
    "StageI_FDR",
    "Canonical_StageI_significant_UP",
    "Early_UP_strength",

    "TN_positive_persists_after_PTPRC",
    "TN_positive_persists_after_ImmuneScore",
    "StageI_positive_persists_after_PTPRC",
    "StageI_positive_persists_after_ImmuneScore",

    "Immune_robustness",

    "LF_Spearman_rho",
    "LF_Spearman_FDR",
    "CIBERSORT_strongest_signals"
  )
]

write.csv(
  ranking_out,
  file.path(outdir, "JAK3_tumor_priority_ranking.csv"),
  row.names = FALSE
)

# ============================================================
# HEATMAP
# ============================================================

hm <- ranking_out

hm$Cancer_label <- ifelse(
  hm$Core_candidate,
  paste0(hm$Cancer, " *"),
  hm$Cancer
)

hm$Cancer_label <- factor(
  hm$Cancer_label,
  levels = rev(hm$Cancer_label)
)

heat <- rbind(

  data.frame(
    Cancer = hm$Cancer_label,
    Feature = "Tumor UP",
    Fill = hm$Tumor_UP_strength,
    Label = sprintf(
      "%.2f",
      hm$TN_median_log2TPM1_difference
    )
  ),

  data.frame(
    Cancer = hm$Cancer_label,
    Feature = "Stage-I UP",
    Fill = hm$Early_UP_strength,
    Label = sprintf(
      "%.2f",
      hm$StageI_median_log2TPM1_difference
    )
  ),

  data.frame(
    Cancer = hm$Cancer_label,
    Feature = "Immune robustness",
    Fill = hm$Immune_robustness,
    Label = sprintf(
      "%.2f",
      hm$Immune_robustness
    )
  ),

  data.frame(
    Cancer = hm$Cancer_label,
    Feature = "Priority score",
    Fill = hm$Priority_score / 100,
    Label = sprintf(
      "%.1f",
      hm$Priority_score
    )
  )
)

heat$Feature <- factor(
  heat$Feature,
  levels = c(
    "Tumor UP",
    "Stage-I UP",
    "Immune robustness",
    "Priority score"
  )
)

p <- ggplot(
  heat,
  aes(
    x = Feature,
    y = Cancer,
    fill = Fill
  )
) +
  geom_tile(
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = Label),
    size = 3.5
  ) +
  scale_fill_gradient(
    limits = c(0, 1),
    name = "Relative\npriority"
  ) +
  labs(
    title = "JAK3 TCGA tumor prioritization",
    subtitle = paste(
      "Tumor upregulation, early-stage upregulation,",
      "and persistence after immune adjustment"
    ),
    x = NULL,
    y = NULL,
    caption = "* passes all conservative core criteria"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

ggsave(
  file.path(
    figdir,
    "Figure_8_JAK3_priority_heatmap.png"
  ),
  p,
  width = 9,
  height = 8,
  dpi = 320
)

ggsave(
  file.path(
    figdir,
    "Figure_8_JAK3_priority_heatmap.pdf"
  ),
  p,
  width = 9,
  height = 8
)

# ============================================================
# Console summary
# ============================================================

cat("\n========================================\n")
cat("JAK3 TCGA TUMOR PRIORITY RANKING\n")
cat("========================================\n\n")

cat(
  "Ranking-eligible cancers:",
  nrow(ranking_out),
  "\n"
)

cat(
  "Strict core candidates:",
  sum(ranking_out$Core_candidate),
  "\n\n"
)

cat("TOP RANKED CANCERS\n\n")

print(
  ranking_out[
    seq_len(min(10, nrow(ranking_out))),
    c(
      "Priority_rank",
      "Cancer",
      "Priority_score",
      "Core_candidate",
      "TN_median_log2TPM1_difference",
      "StageI_median_log2TPM1_difference",
      "Immune_robustness",
      "LF_Spearman_rho"
    )
  ],
  row.names = FALSE
)

cat("\nSaved:\n")
cat("  ", file.path(outdir, "JAK3_tumor_priority_ranking.csv"), "\n")
cat("  ", file.path(figdir, "Figure_8_JAK3_priority_heatmap.png"), "\n")

cat("\n========================================\n")
cat("56 COMPLETE\n")
cat("========================================\n")
