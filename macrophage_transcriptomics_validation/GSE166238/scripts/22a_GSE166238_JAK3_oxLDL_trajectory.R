
rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

setDTthreads(1)

cat("\n====================================================\n")
cat("22a — GSE166238 JAK3 oxLDL TRAJECTORY\n")
cat("DAY1 -> DAY6 POST-WASHOUT\n")
cat("====================================================\n\n")

indir <- "macrophage_transcriptomics_validation/GSE166238/data/processed"
outdir <- "macrophage_transcriptomics_validation/GSE166238/results"
figdir <- "macrophage_transcriptomics_validation/GSE166238/figures"

helper <- "macrophage_transcriptomics_validation/GSE168468/scripts/mmseq.R"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(helper)) {
  stop("Missing MMSEQ helper file: ", helper)
}

source(helper)

JAK3_ID <- "ENSG00000105639"

# ============================================================
# 1. EXACT SAMPLE MANIFEST
#
# Primary analysis:
# Day1 RPMI vs oxLDL
# Day6 RPMI vs oxLDL
#
# Day6 + LPS restimulation samples are excluded.
# ============================================================

samples <- data.frame(

  donor = rep(
    c("DA", "DB", "DC"),
    each = 4
  ),

  time = factor(
    rep(
      c("D1", "D1", "D6", "D6"),
      times = 3
    ),
    levels = c("D1", "D6")
  ),

  condition = factor(
    rep(
      c("RPMI", "oxLDL", "RPMI", "oxLDL"),
      times = 3
    ),
    levels = c("RPMI", "oxLDL")
  ),

  file = c(

    # Donor A
    "GSM5066679_DA_rpmi_d1_trimmed.gene.mmseq.txt.gz",
    "GSM5066680_DA_oxLDL_d1_trimmed.gene.mmseq.txt.gz",
    "GSM5066681_DA_rpmi_d6_trimmed.gene.mmseq.txt.gz",
    "GSM5066682_DA_oxLDL_d6_trimmed.gene.mmseq.txt.gz",

    # Donor B
    "GSM5066686_DB_rpmi_d1_trimmed.gene.mmseq.txt.gz",
    "GSM5066687_DB_oxLDL_d1_trimmed.gene.mmseq.txt.gz",
    "GSM5066688_DB_rpmi_d6_trimmed.gene.mmseq.txt.gz",
    "GSM5066689_DB_oxLDL_d6_trimmed.gene.mmseq.txt.gz",

    # Donor C
    "GSM5066693_DC_rpmi_d1_trimmed.gene.mmseq.txt.gz",
    "GSM5066694_DC_oxLDL_d1_trimmed.gene.mmseq.txt.gz",
    "GSM5066695_DC_rpmi_d6_trimmed.gene.mmseq.txt.gz",
    "GSM5066696_DC_oxLDL_d6_trimmed.gene.mmseq.txt.gz"
  ),

  stringsAsFactors = FALSE
)

samples$path <- file.path(
  indir,
  samples$file
)

if (!all(file.exists(samples$path))) {
  stop(
    paste(
      "Missing files:",
      paste(
        samples$path[!file.exists(samples$path)],
        collapse = "\n"
      )
    )
  )
}

samples$sample_name <- paste(
  samples$donor,
  samples$time,
  samples$condition,
  sep = "_"
)

cat("Samples:", nrow(samples), "\n")
cat("Donors:", length(unique(samples$donor)), "\n")

if (nrow(samples) != 12) {
  stop("Expected exactly 12 primary trajectory samples.")
}

# ============================================================
# 2. JOINT MMSEQ NORMALIZATION
# ============================================================

cat("\nRunning joint MMSEQ normalization...\n")

mm <- readmmseq(
  mmseq_files = samples$path,
  sample_names = samples$sample_name,
  normalize = TRUE
)

expr <- mm$log_mu

if (ncol(expr) != 12) {
  stop("Unexpected normalized sample count.")
}

if (!JAK3_ID %in% rownames(expr)) {
  stop("JAK3 not found.")
}

cat(
  "Normalized matrix:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "samples\n"
)

# ============================================================
# 3. NORMALIZED JAK3 VALUES
# ============================================================

jak3 <- data.table(
  donor = samples$donor,
  time = samples$time,
  condition = samples$condition,
  sample = samples$sample_name,
  normalized_log_mu =
    as.numeric(expr[JAK3_ID, ])
)

fwrite(
  jak3,
  file.path(
    outdir,
    "22a_GSE166238_JAK3_normalized_values.csv"
  )
)

cat("\n===== NORMALIZED JAK3 VALUES =====\n")
print(jak3)

# ============================================================
# 4. FILTER LOW-INFORMATION GENES
#
# Require observed signal in >=3 of 12 samples.
# ============================================================

observed <- mm$observed

keep <- rowSums(
  observed == 1,
  na.rm = TRUE
) >= 3

cat("\nGenes before filtering:", nrow(expr), "\n")
cat("Genes retained:", sum(keep), "\n")

if (!keep[JAK3_ID]) {
  stop("JAK3 failed observed-sample filter.")
}

expr_f <- expr[
  keep,
  ,
  drop = FALSE
]

# ============================================================
# 5. GROUP + DONOR DESIGN
# ============================================================

group <- factor(
  paste(
    samples$time,
    samples$condition,
    sep = "_"
  ),
  levels = c(
    "D1_RPMI",
    "D1_oxLDL",
    "D6_RPMI",
    "D6_oxLDL"
  )
)

donor <- factor(
  samples$donor,
  levels = c("DA", "DB", "DC")
)

design <- model.matrix(
  ~ 0 + group + donor
)

colnames(design) <- sub(
  "^group",
  "",
  colnames(design)
)

cat("\n===== DESIGN COEFFICIENTS =====\n")
print(colnames(design))

if (qr(design)$rank != ncol(design)) {
  stop("Design matrix is not full rank.")
}

# ============================================================
# 6. CONTRASTS
#
# Acute:
# Day1 oxLDL vs Day1 RPMI
#
# Persistent:
# Day6 oxLDL vs Day6 RPMI
#
# Trajectory:
# Day6 treatment effect minus Day1 treatment effect
# ============================================================

contrast_matrix <- makeContrasts(

  oxLDL_D1 =
    D1_oxLDL -
    D1_RPMI,

  oxLDL_D6 =
    D6_oxLDL -
    D6_RPMI,

  Trajectory_D6_vs_D1 =
    (
      D6_oxLDL -
      D6_RPMI
    ) -
    (
      D1_oxLDL -
      D1_RPMI
    ),

  levels = design
)

# ============================================================
# 7. LIMMA FIT ACROSS ALL RETAINED GENES
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
# 8. EXTRACT RESULTS
# ============================================================

all_results <- list()
jak3_results <- list()

for (contrast_name in colnames(contrast_matrix)) {

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
    stop("Expected one JAK3 row for ", contrast_name)
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
    "22a_GSE166238_JAK3_trajectory_statistics.csv"
  )
)

cat("\n========================================\n")
cat("JAK3 RESULTS\n")
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
# 9. DONOR-WISE DESCRIPTIVE EFFECTS
# ============================================================

wide <- dcast(
  jak3,
  donor + time ~ condition,
  value.var = "normalized_log_mu"
)

wide[
  ,
  oxLDL_minus_RPMI :=
    oxLDL - RPMI
]

fwrite(
  wide,
  file.path(
    outdir,
    "22a_GSE166238_JAK3_donor_effects.csv"
  )
)

cat("\n===== DONOR-WISE JAK3 EFFECTS =====\n")
print(wide)

# ============================================================
# 10. DIRECTION CONSISTENCY
# ============================================================

direction <- wide[
  ,
  .(
    mean_effect =
      mean(oxLDL_minus_RPMI),

    median_effect =
      median(oxLDL_minus_RPMI),

    positive_donors =
      sum(oxLDL_minus_RPMI > 0),

    negative_donors =
      sum(oxLDL_minus_RPMI < 0)
  ),
  by = time
]

direction[
  ,
  direction_consistency :=
    pmax(
      positive_donors,
      negative_donors
    ) / 3
]

fwrite(
  direction,
  file.path(
    outdir,
    "22a_GSE166238_JAK3_direction_consistency.csv"
  )
)

cat("\n===== DIRECTION CONSISTENCY =====\n")
print(direction)

# ============================================================
# 11. FIGURE 1 — PAIRED EXPRESSION
# ============================================================

plotdat <- copy(jak3)

plotdat$time <- factor(
  plotdat$time,
  levels = c("D1", "D6"),
  labels = c(
    "24 h",
    "Day 6 post-washout"
  )
)

plotdat$condition <- factor(
  plotdat$condition,
  levels = c(
    "RPMI",
    "oxLDL"
  )
)

p1 <- ggplot(
  plotdat,
  aes(
    x = condition,
    y = normalized_log_mu,
    group = donor
  )
) +
  geom_line(
    alpha = 0.65
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ time,
    nrow = 1
  ) +
  labs(
    title =
      "GSE166238: JAK3 response to oxLDL training",
    subtitle =
      "Paired measurements across three independent donors",
    x = NULL,
    y = "MMSEQ-normalized JAK3 log expression"
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "22a_GSE166238_JAK3_paired_expression.png"
  ),
  p1,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "22a_GSE166238_JAK3_paired_expression.pdf"
  ),
  p1,
  width = 8,
  height = 5
)

# ============================================================
# 12. FIGURE 2 — EFFECT SIZES
# ============================================================

forest <- copy(
  jak3_stats[
    contrast %in%
      c(
        "oxLDL_D1",
        "oxLDL_D6"
      )
  ]
)

forest[
  ,
  time_label :=
    fifelse(
      contrast == "oxLDL_D1",
      "24 h",
      "Day 6 post-washout"
    )
]

forest$time_label <- factor(
  forest$time_label,
  levels = c(
    "24 h",
    "Day 6 post-washout"
  )
)

p2 <- ggplot(
  forest,
  aes(
    x = logFC,
    y = time_label
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
      y = time_label,
      yend = time_label
    ),
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    title =
      "GSE166238: JAK3 oxLDL treatment effects",
    subtitle =
      "Moderated effect estimates with 95% confidence intervals",
    x =
      "oxLDL vs RPMI treatment effect",
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "22a_GSE166238_JAK3_effect_sizes.png"
  ),
  p2,
  width = 7,
  height = 4.5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "22a_GSE166238_JAK3_effect_sizes.pdf"
  ),
  p2,
  width = 7,
  height = 4.5
)

# ============================================================
# 13. METADATA
# ============================================================

metadata <- data.frame(
  dataset = "GSE166238",
  species = "Homo sapiens",
  biological_donors = 3,
  stimulus = "oxLDL",
  initial_exposure = "24 h",
  washout = "Yes",
  recovery = "5 days",
  late_state = "Day6",
  expression_source =
    "Deposited gene-level MMSEQ files",
  normalization =
    "MMSEQ readmmseq joint normalization",
  statistical_model =
    "limma group-means model with donor fixed effect",
  primary_tests =
    paste(
      "Day1 oxLDL vs RPMI;",
      "Day6 oxLDL vs RPMI;",
      "Day6-vs-Day1 interaction"
    ),
  restimulation_samples_in_primary_model =
    FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  metadata,
  file.path(
    outdir,
    "22a_GSE166238_analysis_metadata.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("22a COMPLETE\n")
cat("========================================\n")

