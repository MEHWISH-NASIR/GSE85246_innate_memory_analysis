
rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
})

outdir <- "results/18_JAK3_EPHB2_reconciliation/final"
figdir <- "figures/18_JAK3_EPHB2_reconciliation"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD H3K27ac + H3K4me1 RESULTS
# ============================================================

h27 <- read.csv(
  "results/18_JAK3_EPHB2_reconciliation/H3K27ac/18b_PRIMARY_JAK3_EPHB2_H3K27ac_summary.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

h4 <- read.csv(
  "results/18_JAK3_EPHB2_reconciliation/H3K4me1/18b_PRIMARY_JAK3_EPHB2_H3K4me1_summary.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

h27$Mark <- "H3K27ac"
h4$Mark <- "H3K4me1"

res <- rbind(h27, h4)

# ============================================================
# 2. KEEP CORE FIELDS
# ============================================================

summary_table <- res[
  ,
  c(
    "SYMBOL",
    "Mark",
    "contrast",
    "logFC",
    "CI_low",
    "CI_high",
    "P_value",
    "FDR_target2"
  )
]

names(summary_table) <- c(
  "Gene",
  "Mark",
  "Contrast",
  "Effect",
  "CI_low",
  "CI_high",
  "P_value",
  "FDR_target2"
)

# ============================================================
# 3. ADD INTERPRETATION
# ============================================================

summary_table$Interpretation <- ""

summary_table$Interpretation[
  summary_table$Gene == "JAK3" &
  summary_table$Mark == "H3K27ac" &
  summary_table$Contrast == "Day1_effect"
] <- "Positive Day1 LPS effect"

summary_table$Interpretation[
  summary_table$Gene == "JAK3" &
  summary_table$Mark == "H3K27ac" &
  summary_table$Contrast == "Day6_effect"
] <- "Positive Day6 LPS effect"

summary_table$Interpretation[
  summary_table$Gene == "JAK3" &
  summary_table$Mark == "H3K27ac" &
  summary_table$Contrast == "Trajectory_change"
] <- "No evidence of trajectory change; consistent with persistence"

summary_table$Interpretation[
  summary_table$Gene == "EPHB2" &
  summary_table$Mark == "H3K27ac" &
  summary_table$Contrast == "Trajectory_change"
] <- "Directional upward shift, but CI overlaps zero"

summary_table$Interpretation[
  summary_table$Gene == "JAK3" &
  summary_table$Mark == "H3K4me1" &
  summary_table$Contrast == "Trajectory_change"
] <- "Directional attenuation, uncertain"

summary_table$Interpretation[
  summary_table$Gene == "EPHB2" &
  summary_table$Mark == "H3K4me1" &
  summary_table$Contrast == "Trajectory_change"
] <- "Directional attenuation, uncertain"

write.csv(
  summary_table,
  file.path(
    outdir,
    "18d_Tobias_reconciliation_results.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 4. FINAL CLASSIFICATION TABLE
# ============================================================

classification <- data.frame(

  Gene = c(
    "JAK3",
    "EPHB2"
  ),

  H3K27ac = c(
    "Persistent positive LPS-associated signal",
    "Upward post-washout shift, trajectory uncertain"
  ),

  H3K4me1 = c(
    "Positive Day1 signal; weaker/uncertain Day6",
    "Positive Day1 signal; absent/uncertain Day6"
  ),

  Final_classification = c(
    "Persistent H3K27ac response",
    "Candidate post-washout shift; not confirmed washout-emergent chromatin"
  ),

  stringsAsFactors = FALSE
)

write.csv(
  classification,
  file.path(
    outdir,
    "18d_final_JAK3_EPHB2_classification.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 5. EXPLICIT DEFINITIONS
# ============================================================

definitions <- data.frame(

  Term = c(
    "Persistent",
    "Washout-emergent"
  ),

  Definition = c(
    paste(
      "A treatment-associated effect is retained at the post-washout time point.",
      "Persistence is assessed from effect sizes and uncertainty;",
      "it is not inferred from significance labels alone."
    ),

    paste(
      "The post-washout treatment effect becomes more positive than the earlier effect.",
      "Strong evidence requires support from the direct trajectory interaction",
      "(Day6 treatment effect minus Day1 treatment effect),",
      "rather than significance at only one time point."
    )
  ),

  stringsAsFactors = FALSE
)

write.csv(
  definitions,
  file.path(
    outdir,
    "18d_operational_definitions.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 6. FOREST-PLOT STYLE FIGURE
# ============================================================

plot_df <- summary_table[
  summary_table$Contrast %in%
    c(
      "Day1_effect",
      "Day6_effect",
      "Trajectory_change"
    ),
]

plot_df$Contrast <- factor(
  plot_df$Contrast,
  levels = c(
    "Day1_effect",
    "Day6_effect",
    "Trajectory_change"
  ),
  labels = c(
    "Day 1: LPS - RPMI",
    "Day 6: LPS - RPMI",
    "Trajectory change"
  )
)

p <- ggplot(
  plot_df,
  aes(
    x = Effect,
    y = Contrast
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    height = 0.15
  ) +
  geom_point(
    size = 3
  ) +
  facet_grid(
    Gene ~ Mark
  ) +
  labs(
    title = "JAK3 and EPHB2 chromatin trajectory reconciliation",
    subtitle = "Effect sizes with 95% confidence intervals",
    x = "Effect on log2(signal + 1)",
    y = NULL
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(
    figdir,
    "Figure_18_JAK3_EPHB2_trajectory_reconciliation.png"
  ),
  p,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "Figure_18_JAK3_EPHB2_trajectory_reconciliation.pdf"
  ),
  p,
  width = 10,
  height = 6
)

# ============================================================
# 7. WRITE SHORT INTERPRETATION
# ============================================================

txt <- c(
  "# JAK3 / EPHB2 chromatin reconciliation",
  "",
  "## JAK3",
  "",
  "JAK3 shows a positive H3K27ac LPS-versus-RPMI effect at both Day 1 and Day 6.",
  "The direct Day6-minus-Day1 treatment interaction is approximately zero and its confidence interval includes zero.",
  "This pattern is consistent with persistence of the H3K27ac response rather than a changing trajectory.",
  "H3K4me1 is positive at Day 1 but weaker and uncertain at Day 6.",
  "",
  "## EPHB2",
  "",
  "EPHB2 H3K27ac shifts upward from Day 1 to Day 6, but the direct trajectory contrast remains uncertain.",
  "H3K4me1 does not show a corresponding post-washout emergence and instead decreases toward zero.",
  "Therefore EPHB2 should not be described as confirmed washout-emergent chromatin.",
  "",
  "## Interpretation rule",
  "",
  "Trajectory labels are based on estimated treatment effects and the direct Day-by-treatment interaction.",
  "Different significance labels at two time points are not treated as evidence of different trajectories."
)

writeLines(
  txt,
  file.path(
    outdir,
    "18d_Tobias_reconciliation_interpretation.md"
  )
)

cat("\n========================================\n")
cat("FINAL CLASSIFICATION\n")
cat("========================================\n")
print(classification, row.names = FALSE)

cat("\nFiles written to:\n")
cat(outdir, "\n")
cat(figdir, "\n")

cat("\n18d COMPLETE\n")

