# ============================================================
# 55 — JAK3 FINAL TCGA PAN-CANCER FIGURES
#
# Uses canonical release tables only.
# No raw TCGA data are required.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(
  library(ggplot2)
)

indir <- "TCGA_JAK3_pan_cancer/results/release"
outdir <- "TCGA_JAK3_pan_cancer/figures/final"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clean_project <- function(x) {
  sub("^TCGA-", "", x)
}

save_plot <- function(p, name, width = 9, height = 7) {

  ggsave(
    file.path(outdir, paste0(name, ".png")),
    plot = p,
    width = width,
    height = height,
    dpi = 320
  )

  ggsave(
    file.path(outdir, paste0(name, ".pdf")),
    plot = p,
    width = width,
    height = height
  )
}

theme_final <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )

# ============================================================
# FIGURE 1 — TUMOR VS NORMAL
# ============================================================

tn <- read.csv(
  file.path(indir, "JAK3_complete_tumor_normal.csv"),
  check.names = FALSE
)

tn$Cancer <- clean_project(tn$project)
tn$Significant <- tn$Wilcoxon_FDR < 0.05
tn <- tn[order(tn$Median_log2TPM1_difference), ]
tn$Cancer <- factor(tn$Cancer, levels = tn$Cancer)

p1 <- ggplot(
  tn,
  aes(
    x = Median_log2TPM1_difference,
    y = Cancer,
    shape = Significant
  )
) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_segment(
    aes(
      x = 0,
      xend = Median_log2TPM1_difference,
      yend = Cancer
    ),
    linewidth = 0.6
  ) +
  geom_point(size = 3) +
  labs(
    title = "JAK3 expression: primary tumor versus normal",
    subtitle = "Positive values indicate higher median JAK3 expression in tumor",
    x = "Median log2(TPM + 1) difference",
    y = NULL,
    shape = "FDR < 0.05"
  ) +
  theme_final

save_plot(
  p1,
  "Figure_1_JAK3_tumor_vs_normal",
  width = 8,
  height = 8
)

# ============================================================
# FIGURE 2 — MATCHED TUMOR-NORMAL
# ============================================================

matched <- read.csv(
  file.path(indir, "JAK3_complete_matched_tumor_normal.csv"),
  check.names = FALSE
)

matched$Cancer <- clean_project(matched$project)
matched$Significant <- matched$Paired_Wilcoxon_FDR < 0.05
matched <- matched[order(matched$median_log2_difference), ]
matched$Cancer <- factor(matched$Cancer, levels = matched$Cancer)

p2 <- ggplot(
  matched,
  aes(
    x = median_log2_difference,
    y = Cancer,
    shape = Significant
  )
) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_segment(
    aes(
      x = 0,
      xend = median_log2_difference,
      yend = Cancer
    ),
    linewidth = 0.6
  ) +
  geom_point(size = 3) +
  labs(
    title = "JAK3 matched tumor-normal expression",
    subtitle = "Patient-matched tumor minus normal expression",
    x = "Median paired log2(TPM + 1) difference",
    y = NULL,
    shape = "FDR < 0.05"
  ) +
  theme_final

save_plot(
  p2,
  "Figure_2_JAK3_matched_tumor_normal",
  width = 8,
  height = 8
)

# ============================================================
# FIGURE 3 — STAGE-I VS NORMAL
# ============================================================

early <- read.csv(
  file.path(indir, "JAK3_complete_stageI_vs_normal.csv"),
  check.names = FALSE
)

early$Cancer <- clean_project(early$project)
early$Significant <- early$Wilcoxon_FDR < 0.05
early <- early[order(early$Median_log2TPM1_difference), ]
early$Cancer <- factor(early$Cancer, levels = early$Cancer)

p3 <- ggplot(
  early,
  aes(
    x = Median_log2TPM1_difference,
    y = Cancer,
    shape = Significant
  )
) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_segment(
    aes(
      x = 0,
      xend = Median_log2TPM1_difference,
      yend = Cancer
    ),
    linewidth = 0.6
  ) +
  geom_point(size = 3) +
  labs(
    title = "JAK3 expression: Stage-I tumor versus normal",
    subtitle = "Stage I includes IA/IB/IC where applicable",
    x = "Median log2(TPM + 1) difference",
    y = NULL,
    shape = "FDR < 0.05"
  ) +
  theme_final

save_plot(
  p3,
  "Figure_3_JAK3_stageI_vs_normal",
  width = 8,
  height = 7
)

# ============================================================
# FIGURE 4 — STAGE TREND
# ============================================================

stage <- read.csv(
  file.path(indir, "JAK3_stage_statistics_complete.csv"),
  check.names = FALSE
)

stage$Cancer <- clean_project(stage$project)
stage$Significant <- stage$Stage_trend_FDR < 0.05
stage <- stage[order(stage$Spearman_rho), ]
stage$Cancer <- factor(stage$Cancer, levels = stage$Cancer)

p4 <- ggplot(
  stage,
  aes(
    x = Spearman_rho,
    y = Cancer,
    shape = Significant
  )
) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_segment(
    aes(
      x = 0,
      xend = Spearman_rho,
      yend = Cancer
    ),
    linewidth = 0.6
  ) +
  geom_point(size = 3) +
  labs(
    title = "JAK3 stage-wise expression trends",
    x = "Spearman rho",
    y = NULL,
    shape = "FDR < 0.05"
  ) +
  theme_final

save_plot(
  p4,
  "Figure_4_JAK3_stage_trends",
  width = 8,
  height = 7
)

# ============================================================
# FIGURE 5 — GRADE TREND
# ============================================================

grade <- read.csv(
  file.path(indir, "JAK3_grade_statistics_complete.csv"),
  check.names = FALSE
)

grade$Cancer <- clean_project(grade$project)
grade$Significant <- grade$Grade_trend_FDR < 0.05
grade <- grade[order(grade$Spearman_rho), ]
grade$Cancer <- factor(grade$Cancer, levels = grade$Cancer)

p5 <- ggplot(
  grade,
  aes(
    x = Spearman_rho,
    y = Cancer,
    shape = Significant
  )
) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_segment(
    aes(
      x = 0,
      xend = Spearman_rho,
      yend = Cancer
    ),
    linewidth = 0.6
  ) +
  geom_point(size = 3) +
  labs(
    title = "JAK3 grade-wise expression trends",
    x = "Spearman rho",
    y = NULL,
    shape = "FDR < 0.05"
  ) +
  theme_final

save_plot(
  p5,
  "Figure_5_JAK3_grade_trends",
  width = 8,
  height = 6
)

# ============================================================
# FIGURE 6 — SURVIVAL FOREST
# ============================================================

surv <- read.csv(
  file.path(indir, "JAK3_TCGA_CDR_survival_complete.csv"),
  check.names = FALSE
)

surv <- surv[
  surv$Analysis_status == "Analysed" &
  !is.na(surv$HR) &
  !is.na(surv$CI95_low) &
  !is.na(surv$CI95_high),
]

surv$Cancer <- clean_project(surv$project)
surv$Significant <- surv$Cox_FDR < 0.05

surv <- surv[order(surv$HR), ]
surv$Cancer <- factor(surv$Cancer, levels = surv$Cancer)

p6 <- ggplot(
  surv,
  aes(
    x = HR,
    y = Cancer,
    shape = Significant
  )
) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbarh(
    aes(
      xmin = CI95_low,
      xmax = CI95_high
    ),
    height = 0.18
  ) +
  geom_point(size = 2.7) +
  scale_x_log10() +
  labs(
    title = "JAK3 pan-cancer survival associations",
    subtitle = "Continuous log2(TPM + 1) Cox models using curated TCGA-CDR endpoints",
    x = "Hazard ratio (log scale)",
    y = NULL,
    shape = "FDR < 0.05"
  ) +
  theme_final

save_plot(
  p6,
  "Figure_6_JAK3_survival_forest",
  width = 9,
  height = 9
)

# ============================================================
# FIGURE 7 — LEUKOCYTE FRACTION
# ============================================================

lf <- read.csv(
  file.path(indir, "JAK3_PanImmune_leukocyte_correlations.csv"),
  check.names = FALSE
)

lf$Cancer <- clean_project(lf$project)
lf$Significant <- lf$Spearman_FDR < 0.05

lf <- lf[order(lf$Spearman_rho), ]
lf$Cancer <- factor(lf$Cancer, levels = lf$Cancer)

p7 <- ggplot(
  lf,
  aes(
    x = Spearman_rho,
    y = Cancer,
    shape = Significant
  )
) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_segment(
    aes(
      x = 0,
      xend = Spearman_rho,
      yend = Cancer
    ),
    linewidth = 0.6
  ) +
  geom_point(size = 3) +
  labs(
    title = "Association between JAK3 expression and leukocyte fraction",
    subtitle = "PanImmune leukocyte-fraction estimates",
    x = "Spearman rho",
    y = NULL,
    shape = "FDR < 0.05"
  ) +
  theme_final

save_plot(
  p7,
  "Figure_7_JAK3_leukocyte_fraction",
  width = 8,
  height = 9
)

cat("\n========================================\n")
cat("55 FINAL JAK3 FIGURES COMPLETE\n")
cat("========================================\n\n")

print(
  list.files(
    outdir,
    pattern = "\\.(png|pdf)$",
    full.names = TRUE
  )
)
