
rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

setDTthreads(1)

cat("\n=====================================================\n")
cat("20d — GSE168468 JAK3 TIME-COURSE TRAJECTORY ANALYSIS\n")
cat("4h -> 24h -> Day6 post-washout\n")
cat("=====================================================\n\n")

indir <- "macrophage_transcriptomics_validation/GSE168468/data/processed"
outdir <- "macrophage_transcriptomics_validation/GSE168468/results"
figdir <- "macrophage_transcriptomics_validation/GSE168468/figures"

helper <- paste0(
  "macrophage_transcriptomics_validation/",
  "GSE168468/scripts/mmseq.R"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(helper)) {
  stop("Missing MMSEQ helper: ", helper)
}

source(helper)

JAK3_ID <- "ENSG00000105639"

# ============================================================
# 1. EXACT SAMPLE MANIFEST
#
# 5 donors:
# A1-A5
#
# Three primary timepoints:
# 4 h
# 24 h
# Day 6 after removal of the initial stimulus
#
# Four conditions:
# RPMI
# LPS
# BCG Denmark
# BCG Bulgaria
#
# Day6 + 4h LPS restimulation samples are intentionally
# EXCLUDED from this primary persistence analysis.
# ============================================================

donors <- paste0("A", 1:5)

add_group <- function(
  time,
  condition,
  files
) {

  if (length(files) != 5) {
    stop("Each group must contain exactly five donors.")
  }

  data.frame(
    donor = donors,
    time = time,
    condition = condition,
    file = files,
    stringsAsFactors = FALSE
  )
}

samples <- rbind(

  # ----------------------------------------------------------
  # 4 HOURS
  # ----------------------------------------------------------

  add_group(
    "T4h",
    "RPMI",
    c(
      "GSM5141074_A1_R4h.gene.mmseq.txt.gz",
      "GSM5141075_A2_R4h.gene.mmseq.txt.gz",
      "GSM5141076_A3_R4h.gene.mmseq.txt.gz",
      "GSM5141077_A4_R4h.gene.mmseq.txt.gz",
      "GSM5141078_A5_R4h.gene.mmseq.txt.gz"
    )
  ),

  add_group(
    "T4h",
    "LPS",
    c(
      "GSM5141079_A1_L4h.gene.mmseq.txt.gz",
      "GSM5141080_A2_L4h.gene.mmseq.txt.gz",
      "GSM5141081_A3_L4h.gene.mmseq.txt.gz",
      "GSM5141082_A4_L4h.gene.mmseq.txt.gz",
      "GSM5141083_A5_L4h.gene.mmseq.txt.gz"
    )
  ),

  add_group(
    "T4h",
    "BCG_Denmark",
    c(
      "GSM5141084_A1_D4h.gene.mmseq.txt.gz",
      "GSM5141085_A2_D4h.gene.mmseq.txt.gz",
      "GSM5141086_A3_D4h.gene.mmseq.txt.gz",
      "GSM5141087_A4_D4h.gene.mmseq.txt.gz",
      "GSM5141088_A5_D4h.gene.mmseq.txt.gz"
    )
  ),

  add_group(
    "T4h",
    "BCG_Bulgaria",
    c(
      "GSM5141089_A1_B4h.gene.mmseq.txt.gz",
      "GSM5141090_A2_B4h.gene.mmseq.txt.gz",
      "GSM5141091_A3_B4h.gene.mmseq.txt.gz",
      "GSM5141092_A4_B4h.gene.mmseq.txt.gz",
      "GSM5141093_A5_B4h.gene.mmseq.txt.gz"
    )
  ),

  # ----------------------------------------------------------
  # 24 HOURS
  # ----------------------------------------------------------

  add_group(
    "T24h",
    "RPMI",
    c(
      "GSM5141094_A1_24h_RPMI.gene.mmseq.txt.gz",
      "GSM5141095_A2_d1R.gene.mmseq.txt.gz",
      "GSM5141096_A3_d1R.gene.mmseq.txt.gz",
      "GSM5141097_A4_d1R.gene.mmseq.txt.gz",
      "GSM5141098_A5_d1R.gene.mmseq.txt.gz"
    )
  ),

  add_group(
    "T24h",
    "LPS",
    c(
      "GSM5141099_A1_24h_LPS.gene.mmseq.txt.gz",
      "GSM5141100_A2_d1L.gene.mmseq.txt.gz",
      "GSM5141101_A3_d1L.gene.mmseq.txt.gz",
      "GSM5141102_A4_d1L.gene.mmseq.txt.gz",
      "GSM5141103_A5_d1L.gene.mmseq.txt.gz"
    )
  ),

  add_group(
    "T24h",
    "BCG_Denmark",
    c(
      "GSM5141104_A1_24h_Den.gene.mmseq.txt.gz",
      "GSM5141105_A2_d1D.gene.mmseq.txt.gz",
      "GSM5141106_A3_d1D.gene.mmseq.txt.gz",
      "GSM5141107_A4_d1D.gene.mmseq.txt.gz",
      "GSM5141108_A5_d1D.gene.mmseq.txt.gz"
    )
  ),

  add_group(
    "T24h",
    "BCG_Bulgaria",
    c(
      "GSM5141109_A1_24h_Bul.gene.mmseq.txt.gz",
      "GSM5141110_A2_d1B.gene.mmseq.txt.gz",
      "GSM5141111_A3_d1B.gene.mmseq.txt.gz",
      "GSM5141112_A4_d1B.gene.mmseq.txt.gz",
      "GSM5141113_A5_d1B.gene.mmseq.txt.gz"
    )
  ),

  # ----------------------------------------------------------
  # DAY 6 — POST-WASHOUT / RECOVERY STATE
  # ----------------------------------------------------------

  add_group(
    "D6",
    "RPMI",
    c(
      "GSM5141114_A1_D6_RPMI.gene.mmseq.txt.gz",
      "GSM5141115_A2_d6R.gene.mmseq.txt.gz",
      "GSM5141116_A3_d6R.gene.mmseq.txt.gz",
      "GSM5141117_A4_d6R.gene.mmseq.txt.gz",
      "GSM5141118_A5_d6R.gene.mmseq.txt.gz"
    )
  ),

  add_group(
    "D6",
    "LPS",
    c(
      "GSM5141119_A1_d6L.gene.mmseq.txt.gz",
      "GSM5141120_A2_d6L.gene.mmseq.txt.gz",
      "GSM5141121_A3_d6L.gene.mmseq.txt.gz",
      "GSM5141122_A4_d6L.gene.mmseq.txt.gz",
      "GSM5141123_A5_d6L.gene.mmseq.txt.gz"
    )
  ),

  add_group(
    "D6",
    "BCG_Denmark",
    c(
      "GSM5141124_A1_D6_Den.gene.mmseq.txt.gz",
      "GSM5141125_A2_d6D.gene.mmseq.txt.gz",
      "GSM5141126_A3_d6D.gene.mmseq.txt.gz",
      "GSM5141127_A4_d6D.gene.mmseq.txt.gz",
      "GSM5141128_A5_d6D.gene.mmseq.txt.gz"
    )
  ),

  add_group(
    "D6",
    "BCG_Bulgaria",
    c(
      "GSM5141129_A1_D6_Bul.gene.mmseq.txt.gz",
      "GSM5141130_A2_d6B.gene.mmseq.txt.gz",
      "GSM5141131_A3_d6B.gene.mmseq.txt.gz",
      "GSM5141132_A4_d6B.gene.mmseq.txt.gz",
      "GSM5141133_A5_d6B.gene.mmseq.txt.gz"
    )
  )
)

samples$path <- file.path(
  indir,
  samples$file
)

missing_files <- samples$path[
  !file.exists(samples$path)
]

if (length(missing_files) > 0) {

  stop(
    paste(
      "Missing files:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}

cat("Total samples:", nrow(samples), "\n")

cat(
  "Unique donors:",
  length(unique(samples$donor)),
  "\n"
)

cat(
  "Timepoints:",
  paste(unique(samples$time), collapse = ", "),
  "\n"
)

cat(
  "Conditions:",
  paste(unique(samples$condition), collapse = ", "),
  "\n"
)

if (nrow(samples) != 60) {
  stop("Expected exactly 60 samples.")
}

# ============================================================
# 2. JOINT MMSEQ NORMALIZATION
#
# All 60 primary trajectory samples are normalized together.
#
# This follows the MMSEQ multi-sample scaling workflow.
# ============================================================

samples$sample_name <- paste(
  samples$donor,
  samples$time,
  samples$condition,
  sep = "_"
)

cat("\n========================================\n")
cat("RUNNING JOINT MMSEQ NORMALIZATION\n")
cat("========================================\n")

mm <- readmmseq(
  mmseq_files =
    samples$path,

  sample_names =
    samples$sample_name,

  normalize = TRUE
)

expr <- mm$log_mu

if (ncol(expr) != nrow(samples)) {
  stop("Normalized expression/sample count mismatch.")
}

if (!JAK3_ID %in% rownames(expr)) {
  stop("JAK3 is absent from normalized matrix.")
}

cat(
  "\nNormalized expression matrix:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "samples\n"
)

# ============================================================
# 3. SAVE NORMALIZED JAK3 VALUES
# ============================================================

jak3 <- data.table(
  donor = samples$donor,
  time = samples$time,
  condition = samples$condition,
  sample = samples$sample_name,
  normalized_log_mu =
    as.numeric(expr[JAK3_ID, ])
)

jak3[, time := factor(
  time,
  levels = c(
    "T4h",
    "T24h",
    "D6"
  )
)]

jak3[, condition := factor(
  condition,
  levels = c(
    "RPMI",
    "LPS",
    "BCG_Denmark",
    "BCG_Bulgaria"
  )
)]

fwrite(
  jak3,
  file.path(
    outdir,
    "20d_GSE168468_JAK3_normalized_timecourse_values.csv"
  )
)

cat("\n===== NORMALIZED JAK3 TIME-COURSE VALUES =====\n")

print(jak3)

# ============================================================
# 4. GENE FILTERING
#
# Keep genes observed in at least 25% of the 60 samples.
#
# This removes very low-information features before empirical
# Bayes variance estimation.
# ============================================================

observed <- mm$observed

min_observed <- ceiling(
  0.25 * ncol(observed)
)

keep <- rowSums(
  observed == 1,
  na.rm = TRUE
) >= min_observed

cat("\nGenes before filtering:",
    nrow(expr), "\n")

cat(
  "Minimum observed samples required:",
  min_observed,
  "\n"
)

cat(
  "Genes retained:",
  sum(keep),
  "\n"
)

if (!keep[JAK3_ID]) {
  stop("JAK3 failed observed-sample filter.")
}

expr_f <- expr[
  keep,
  ,
  drop = FALSE
]

# ============================================================
# 5. MODEL DESIGN
#
# group captures each Time x Condition combination.
#
# donor is included as a fixed blocking factor because the
# same five donors contribute measurements to every group.
# ============================================================

group_levels <- c(
  "T4h_RPMI",
  "T4h_LPS",
  "T4h_BCG_Denmark",
  "T4h_BCG_Bulgaria",

  "T24h_RPMI",
  "T24h_LPS",
  "T24h_BCG_Denmark",
  "T24h_BCG_Bulgaria",

  "D6_RPMI",
  "D6_LPS",
  "D6_BCG_Denmark",
  "D6_BCG_Bulgaria"
)

group <- factor(
  paste(
    samples$time,
    samples$condition,
    sep = "_"
  ),
  levels = group_levels
)

donor <- factor(
  samples$donor,
  levels = donors
)

design <- model.matrix(
  ~ 0 + group + donor
)

colnames(design) <- sub(
  "^group",
  "",
  colnames(design)
)

cat("\n===== DESIGN MATRIX DIMENSIONS =====\n")
cat(
  nrow(design),
  "samples x",
  ncol(design),
  "coefficients\n"
)

cat(
  "Design rank:",
  qr(design)$rank,
  "\n"
)

if (qr(design)$rank != ncol(design)) {
  stop("Design matrix is not full rank.")
}

cat("\nDesign coefficients:\n")
print(colnames(design))

# ============================================================
# 6. DEFINE TIME-SPECIFIC EFFECTS AND DIRECT TRAJECTORY TESTS
#
# Main effects:
#   Treated - RPMI at each timepoint
#
# Trajectory:
#   difference-of-differences
#
# Example:
#
# (D6 LPS - D6 RPMI)
# -
# (24h LPS - 24h RPMI)
# ============================================================

contrast_matrix <- makeContrasts(

  # ----------------------------------------------------------
  # LPS treatment effects
  # ----------------------------------------------------------

  LPS_4h =
    T4h_LPS -
    T4h_RPMI,

  LPS_24h =
    T24h_LPS -
    T24h_RPMI,

  LPS_D6 =
    D6_LPS -
    D6_RPMI,

  # ----------------------------------------------------------
  # Denmark treatment effects
  # ----------------------------------------------------------

  Denmark_4h =
    T4h_BCG_Denmark -
    T4h_RPMI,

  Denmark_24h =
    T24h_BCG_Denmark -
    T24h_RPMI,

  Denmark_D6 =
    D6_BCG_Denmark -
    D6_RPMI,

  # ----------------------------------------------------------
  # Bulgaria treatment effects
  # ----------------------------------------------------------

  Bulgaria_4h =
    T4h_BCG_Bulgaria -
    T4h_RPMI,

  Bulgaria_24h =
    T24h_BCG_Bulgaria -
    T24h_RPMI,

  Bulgaria_D6 =
    D6_BCG_Bulgaria -
    D6_RPMI,

  # ----------------------------------------------------------
  # LPS trajectory contrasts
  # ----------------------------------------------------------

  LPS_24h_vs_4h =
    (
      T24h_LPS -
      T24h_RPMI
    ) -
    (
      T4h_LPS -
      T4h_RPMI
    ),

  LPS_D6_vs_24h =
    (
      D6_LPS -
      D6_RPMI
    ) -
    (
      T24h_LPS -
      T24h_RPMI
    ),

  LPS_D6_vs_4h =
    (
      D6_LPS -
      D6_RPMI
    ) -
    (
      T4h_LPS -
      T4h_RPMI
    ),

  # ----------------------------------------------------------
  # Denmark trajectory contrasts
  # ----------------------------------------------------------

  Denmark_24h_vs_4h =
    (
      T24h_BCG_Denmark -
      T24h_RPMI
    ) -
    (
      T4h_BCG_Denmark -
      T4h_RPMI
    ),

  Denmark_D6_vs_24h =
    (
      D6_BCG_Denmark -
      D6_RPMI
    ) -
    (
      T24h_BCG_Denmark -
      T24h_RPMI
    ),

  Denmark_D6_vs_4h =
    (
      D6_BCG_Denmark -
      D6_RPMI
    ) -
    (
      T4h_BCG_Denmark -
      T4h_RPMI
    ),

  # ----------------------------------------------------------
  # Bulgaria trajectory contrasts
  # ----------------------------------------------------------

  Bulgaria_24h_vs_4h =
    (
      T24h_BCG_Bulgaria -
      T24h_RPMI
    ) -
    (
      T4h_BCG_Bulgaria -
      T4h_RPMI
    ),

  Bulgaria_D6_vs_24h =
    (
      D6_BCG_Bulgaria -
      D6_RPMI
    ) -
    (
      T24h_BCG_Bulgaria -
      T24h_RPMI
    ),

  Bulgaria_D6_vs_4h =
    (
      D6_BCG_Bulgaria -
      D6_RPMI
    ) -
    (
      T4h_BCG_Bulgaria -
      T4h_RPMI
    ),

  levels = design
)

cat(
  "\nContrasts defined:",
  ncol(contrast_matrix),
  "\n"
)

# ============================================================
# 7. FIT ALL GENES
#
# Empirical-Bayes moderation is performed across the entire
# retained feature universe BEFORE JAK3 extraction.
# ============================================================

fit0 <- lmFit(
  expr_f,
  design
)

fit1 <- contrasts.fit(
  fit0,
  contrast_matrix
)

fit1 <- eBayes(
  fit1,
  robust = TRUE
)

# ============================================================
# 8. EXTRACT GENOME-WIDE AND JAK3 RESULTS
# ============================================================

all_results <- list()
jak3_results <- list()

for (
  contrast_name
  in colnames(contrast_matrix)
) {

  tt <- topTable(
    fit1,
    coef = contrast_name,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH",
    confint = TRUE
  )

  tt$feature_id <- rownames(tt)
  tt$contrast <- contrast_name

  all_results[[contrast_name]] <- tt

  j <- tt[
    tt$feature_id == JAK3_ID,
    ,
    drop = FALSE
  ]

  if (nrow(j) != 1) {
    stop(
      paste(
        "Expected exactly one JAK3 row for",
        contrast_name
      )
    )
  }

  jak3_results[[contrast_name]] <- j
}

jak3_stats <- rbindlist(
  jak3_results,
  fill = TRUE
)

fwrite(
  jak3_stats,
  file.path(
    outdir,
    "20d_GSE168468_JAK3_trajectory_statistics.csv"
  )
)

cat("\n========================================\n")
cat("JAK3 TIME-SPECIFIC AND TRAJECTORY RESULTS\n")
cat("========================================\n")

print(
  jak3_stats[
    ,
    .(
      contrast,
      logFC,
      CI.L,
      CI.R,
      t,
      P.Value,
      adj.P.Val
    )
  ]
)

# ============================================================
# 9. DONOR-SPECIFIC TREATMENT EFFECTS
# ============================================================

wide <- dcast(
  jak3,
  donor + time ~ condition,
  value.var = "normalized_log_mu"
)

wide[
  ,
  LPS_effect :=
    LPS - RPMI
]

wide[
  ,
  Denmark_effect :=
    BCG_Denmark - RPMI
]

wide[
  ,
  Bulgaria_effect :=
    BCG_Bulgaria - RPMI
]

fwrite(
  wide,
  file.path(
    outdir,
    "20d_GSE168468_JAK3_donor_timepoint_effects.csv"
  )
)

cat("\n===== DONOR-WISE TREATMENT EFFECTS =====\n")
print(wide)

# ============================================================
# 10. DIRECTION CONSISTENCY AT EACH TIMEPOINT
# ============================================================

effect_long <- melt(
  wide,
  id.vars = c(
    "donor",
    "time"
  ),
  measure.vars = c(
    "LPS_effect",
    "Denmark_effect",
    "Bulgaria_effect"
  ),
  variable.name = "stimulus",
  value.name = "effect"
)

effect_long[
  ,
  stimulus := fifelse(
    stimulus == "LPS_effect",
    "LPS",
    fifelse(
      stimulus == "Denmark_effect",
      "BCG_Denmark",
      "BCG_Bulgaria"
    )
  )
]

direction_summary <- effect_long[
  ,
  .(
    mean_effect =
      mean(effect),

    median_effect =
      median(effect),

    positive_donors =
      sum(effect > 0),

    negative_donors =
      sum(effect < 0),

    zero_donors =
      sum(effect == 0),

    direction_consistency =
      max(
        sum(effect > 0),
        sum(effect < 0)
      ) / .N
  ),
  by = .(
    time,
    stimulus
  )
]

fwrite(
  direction_summary,
  file.path(
    outdir,
    "20d_GSE168468_JAK3_direction_consistency.csv"
  )
)

cat("\n===== DIRECTION CONSISTENCY =====\n")
print(direction_summary)

# ============================================================
# 11. FIGURE 1
# DONOR-SPECIFIC TREATMENT EFFECT TRAJECTORIES
# ============================================================

effect_long$time <- factor(
  effect_long$time,
  levels = c(
    "T4h",
    "T24h",
    "D6"
  ),
  labels = c(
    "4 h",
    "24 h",
    "Day 6"
  )
)

effect_long$stimulus <- factor(
  effect_long$stimulus,
  levels = c(
    "LPS",
    "BCG_Denmark",
    "BCG_Bulgaria"
  ),
  labels = c(
    "LPS",
    "BCG Denmark",
    "BCG Bulgaria"
  )
)

p1 <- ggplot(
  effect_long,
  aes(
    x = time,
    y = effect,
    group = donor
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_line(
    alpha = 0.65
  ) +
  geom_point(
    size = 2.8
  ) +
  facet_wrap(
    ~ stimulus,
    nrow = 1
  ) +
  labs(
    title =
      "GSE168468: donor-specific JAK3 treatment-effect trajectories",
    subtitle =
      "Each value is treated minus matched RPMI within the same donor and timepoint",
    x = NULL,
    y =
      "JAK3 treatment effect\n(normalized log expression)"
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "20d_GSE168468_JAK3_donor_effect_trajectories.png"
  ),
  p1,
  width = 10,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "20d_GSE168468_JAK3_donor_effect_trajectories.pdf"
  ),
  p1,
  width = 10,
  height = 5
)

# ============================================================
# 12. FIGURE 2
# MODERATED TREATMENT EFFECTS + 95% CI
# ============================================================

effect_contrasts <- c(
  "LPS_4h",
  "LPS_24h",
  "LPS_D6",
  "Denmark_4h",
  "Denmark_24h",
  "Denmark_D6",
  "Bulgaria_4h",
  "Bulgaria_24h",
  "Bulgaria_D6"
)

forest <- copy(
  jak3_stats[
    contrast %in%
      effect_contrasts
  ]
)

forest[
  ,
  stimulus := fifelse(
    grepl("^LPS", contrast),
    "LPS",
    fifelse(
      grepl("^Denmark", contrast),
      "BCG Denmark",
      "BCG Bulgaria"
    )
  )
]

forest[
  ,
  time := fifelse(
    grepl("_4h$", contrast),
    "4 h",
    fifelse(
      grepl("_24h$", contrast),
      "24 h",
      "Day 6"
    )
  )
]

forest$time <- factor(
  forest$time,
  levels = c(
    "4 h",
    "24 h",
    "Day 6"
  )
)

forest$stimulus <- factor(
  forest$stimulus,
  levels = c(
    "LPS",
    "BCG Denmark",
    "BCG Bulgaria"
  )
)

p2 <- ggplot(
  forest,
  aes(
    x = logFC,
    y = time
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_segment(
    aes(
      x = CI.L,
      xend = CI.R,
      y = time,
      yend = time
    ),
    linewidth = 0.7
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ stimulus,
    nrow = 1
  ) +
  labs(
    title =
      "GSE168468: JAK3 treatment effects across time",
    subtitle =
      "Moderated effect estimates with 95% confidence intervals",
    x =
      "Treatment effect versus matched RPMI",
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "20d_GSE168468_JAK3_effect_sizes_95CI.png"
  ),
  p2,
  width = 10,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "20d_GSE168468_JAK3_effect_sizes_95CI.pdf"
  ),
  p2,
  width = 10,
  height = 5
)

# ============================================================
# 13. TRAJECTORY-ONLY TABLE
# ============================================================

trajectory_stats <- jak3_stats[
  grepl(
    "vs",
    contrast
  )
]

fwrite(
  trajectory_stats,
  file.path(
    outdir,
    "20d_GSE168468_JAK3_direct_trajectory_tests.csv"
  )
)

# ============================================================
# 14. STUDY / MODEL METADATA
# ============================================================

metadata <- data.frame(

  dataset =
    "GSE168468",

  gene =
    "JAK3",

  gene_id =
    JAK3_ID,

  species =
    "Homo sapiens",

  donors =
    5,

  cell_system =
    "Primary human monocytes / differentiated macrophages",

  stimuli =
    "LPS; BCG Denmark; BCG Bulgaria",

  timepoints =
    "4h; 24h; Day6",

  post_washout_state =
    "Day6",

  restimulation_samples_in_primary_model =
    FALSE,

  normalization =
    "MMSEQ readmmseq joint normalization across 60 samples",

  model =
    "limma group-means model with donor fixed blocking factor",

  inference =
    paste(
      "Treatment-vs-RPMI effects plus",
      "direct difference-of-differences trajectory contrasts"
    ),

  empirical_bayes =
    "robust eBayes across filtered gene universe",

  gene_filter =
    paste0(
      "Observed in >= ",
      min_observed,
      " of 60 samples"
    ),

  stringsAsFactors =
    FALSE
)

write.csv(
  metadata,
  file.path(
    outdir,
    "20d_GSE168468_trajectory_analysis_metadata.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 15. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("DAY6 JAK3 EFFECTS\n")
cat("========================================\n")

print(
  jak3_stats[
    contrast %in%
      c(
        "LPS_D6",
        "Denmark_D6",
        "Bulgaria_D6"
      ),
    .(
      contrast,
      logFC,
      CI.L,
      CI.R,
      P.Value,
      adj.P.Val
    )
  ]
)

cat("\n========================================\n")
cat("DIRECT DAY6-vs-24h TRAJECTORY TESTS\n")
cat("========================================\n")

print(
  jak3_stats[
    contrast %in%
      c(
        "LPS_D6_vs_24h",
        "Denmark_D6_vs_24h",
        "Bulgaria_D6_vs_24h"
      ),
    .(
      contrast,
      logFC,
      CI.L,
      CI.R,
      P.Value,
      adj.P.Val
    )
  ]
)

cat("\n========================================\n")
cat("DIRECT DAY6-vs-4h TRAJECTORY TESTS\n")
cat("========================================\n")

print(
  jak3_stats[
    contrast %in%
      c(
        "LPS_D6_vs_4h",
        "Denmark_D6_vs_4h",
        "Bulgaria_D6_vs_4h"
      ),
    .(
      contrast,
      logFC,
      CI.L,
      CI.R,
      P.Value,
      adj.P.Val
    )
  ]
)

cat("\n20d COMPLETE\n")

