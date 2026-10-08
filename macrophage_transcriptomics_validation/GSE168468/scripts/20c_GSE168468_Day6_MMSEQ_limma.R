
rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

setDTthreads(1)

cat("\n========================================\n")
cat("20c — GSE168468 DAY6 MMSEQ-NORMALIZED\n")
cat("PAIRED DONOR LIMMA ANALYSIS\n")
cat("========================================\n\n")

indir <- "macrophage_transcriptomics_validation/GSE168468/data/processed"
outdir <- "macrophage_transcriptomics_validation/GSE168468/results"
figdir <- "macrophage_transcriptomics_validation/GSE168468/figures"
helper <- "macrophage_transcriptomics_validation/GSE168468/scripts/mmseq.R"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(helper)) {
  stop("Missing MMSEQ helper file: ", helper)
}

source(helper)

JAK3_ID <- "ENSG00000105639"

# ============================================================
# 1. SAMPLE MANIFEST
# ============================================================

samples <- data.frame(
  donor = rep(paste0("A", 1:5), 4),

  condition = factor(
    rep(
      c(
        "RPMI",
        "LPS",
        "BCG_Denmark",
        "BCG_Bulgaria"
      ),
      each = 5
    ),
    levels = c(
      "RPMI",
      "LPS",
      "BCG_Denmark",
      "BCG_Bulgaria"
    )
  ),

  file = c(
    # RPMI
    "GSM5141114_A1_D6_RPMI.gene.mmseq.txt.gz",
    "GSM5141115_A2_d6R.gene.mmseq.txt.gz",
    "GSM5141116_A3_d6R.gene.mmseq.txt.gz",
    "GSM5141117_A4_d6R.gene.mmseq.txt.gz",
    "GSM5141118_A5_d6R.gene.mmseq.txt.gz",

    # LPS
    "GSM5141119_A1_d6L.gene.mmseq.txt.gz",
    "GSM5141120_A2_d6L.gene.mmseq.txt.gz",
    "GSM5141121_A3_d6L.gene.mmseq.txt.gz",
    "GSM5141122_A4_d6L.gene.mmseq.txt.gz",
    "GSM5141123_A5_d6L.gene.mmseq.txt.gz",

    # BCG Denmark
    "GSM5141124_A1_D6_Den.gene.mmseq.txt.gz",
    "GSM5141125_A2_d6D.gene.mmseq.txt.gz",
    "GSM5141126_A3_d6D.gene.mmseq.txt.gz",
    "GSM5141127_A4_d6D.gene.mmseq.txt.gz",
    "GSM5141128_A5_d6D.gene.mmseq.txt.gz",

    # BCG Bulgaria
    "GSM5141129_A1_D6_Bul.gene.mmseq.txt.gz",
    "GSM5141130_A2_d6B.gene.mmseq.txt.gz",
    "GSM5141131_A3_d6B.gene.mmseq.txt.gz",
    "GSM5141132_A4_d6B.gene.mmseq.txt.gz",
    "GSM5141133_A5_d6B.gene.mmseq.txt.gz"
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
      paste(samples$path[!file.exists(samples$path)], collapse = "\n")
    )
  )
}

samples$sample_name <- paste(
  samples$donor,
  samples$condition,
  sep = "_"
)

cat("Samples:", nrow(samples), "\n")
cat("Donors:", length(unique(samples$donor)), "\n")

# ============================================================
# 2. OFFICIAL MMSEQ MULTI-SAMPLE NORMALIZATION
# ============================================================

cat("\nRunning readmmseq normalization...\n")

mm <- readmmseq(
  mmseq_files = samples$path,
  sample_names = samples$sample_name,
  normalize = TRUE
)

expr <- mm$log_mu

if (ncol(expr) != nrow(samples)) {
  stop("Normalized matrix sample count mismatch.")
}

if (!JAK3_ID %in% rownames(expr)) {
  stop("JAK3 not found in normalized MMSEQ matrix.")
}

cat(
  "Normalized matrix:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "samples\n"
)

# Save normalized JAK3 values
jak3_norm <- data.frame(
  donor = samples$donor,
  condition = samples$condition,
  normalized_log_mu =
    as.numeric(expr[JAK3_ID, ]),
  stringsAsFactors = FALSE
)

fwrite(
  jak3_norm,
  file.path(
    outdir,
    "20c_GSE168468_Day6_JAK3_normalized_values.csv"
  )
)

cat("\n===== NORMALIZED JAK3 VALUES =====\n")
print(jak3_norm)

# ============================================================
# 3. PAIRED DONOR DESIGN
#
# donor = fixed blocking factor
# condition = treatment effect
# RPMI = reference
# ============================================================

donor <- factor(samples$donor)

condition <- factor(
  samples$condition,
  levels = c(
    "RPMI",
    "LPS",
    "BCG_Denmark",
    "BCG_Bulgaria"
  )
)

design <- model.matrix(
  ~ donor + condition
)

cat("\n===== DESIGN MATRIX =====\n")
print(design)

# ============================================================
# 4. FILTER VERY LOW-INFORMATION GENES
#
# Retain genes observed in >= 5 of 20 samples.
# This is for stable genome-wide variance estimation.
# JAK3 is checked separately and must survive.
# ============================================================

observed_mat <- mm$observed

keep <- rowSums(
  observed_mat == 1,
  na.rm = TRUE
) >= 5

cat("\nGenes before filtering:", nrow(expr), "\n")
cat("Genes retained:", sum(keep), "\n")

if (!keep[JAK3_ID]) {
  stop("JAK3 failed observed-sample filter.")
}

expr_f <- expr[keep, , drop = FALSE]

# ============================================================
# 5. LIMMA FIT ACROSS ALL RETAINED GENES
# ============================================================

fit <- lmFit(
  expr_f,
  design
)

fit <- eBayes(
  fit,
  robust = TRUE
)

cat("\nModel coefficients:\n")
print(colnames(design))

# Expected coefficients:
# conditionLPS
# conditionBCG_Denmark
# conditionBCG_Bulgaria

coef_names <- c(
  "conditionLPS",
  "conditionBCG_Denmark",
  "conditionBCG_Bulgaria"
)

if (!all(coef_names %in% colnames(design))) {
  stop(
    paste(
      "Expected condition coefficients not found:",
      paste(
        setdiff(coef_names, colnames(design)),
        collapse = ", "
      )
    )
  )
}

# ============================================================
# 6. FULL GENOME-WIDE RESULTS + JAK3 EXTRACTION
# ============================================================

all_results <- list()
jak3_results <- list()

contrast_labels <- c(
  conditionLPS =
    "Day6_LPS_vs_RPMI",

  conditionBCG_Denmark =
    "Day6_BCG_Denmark_vs_RPMI",

  conditionBCG_Bulgaria =
    "Day6_BCG_Bulgaria_vs_RPMI"
)

for (coef in coef_names) {

  tt <- topTable(
    fit,
    coef = coef,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH",
    confint = TRUE
  )

  tt$feature_id <- rownames(tt)
  tt$contrast <- contrast_labels[[coef]]

  all_results[[coef]] <- tt

  j <- tt[
    tt$feature_id == JAK3_ID,
    ,
    drop = FALSE
  ]

  if (nrow(j) != 1) {
    stop("Expected one JAK3 row for ", coef)
  }

  jak3_results[[coef]] <- j
}

genomewide <- rbindlist(
  all_results,
  fill = TRUE
)

jak3_stats <- rbindlist(
  jak3_results,
  fill = TRUE
)

fwrite(
  genomewide,
  file.path(
    outdir,
    "20c_GSE168468_Day6_genomewide_limma.csv"
  )
)

fwrite(
  jak3_stats,
  file.path(
    outdir,
    "20c_GSE168468_Day6_JAK3_limma.csv"
  )
)

cat("\n===== JAK3 MODERATED RESULTS =====\n")

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
# 7. DONOR-WISE NORMALIZED EFFECTS
# ============================================================

wide <- dcast(
  as.data.table(jak3_norm),
  donor ~ condition,
  value.var = "normalized_log_mu"
)

wide[, LPS_minus_RPMI :=
       LPS - RPMI]

wide[, Denmark_minus_RPMI :=
       BCG_Denmark - RPMI]

wide[, Bulgaria_minus_RPMI :=
       BCG_Bulgaria - RPMI]

fwrite(
  wide,
  file.path(
    outdir,
    "20c_GSE168468_Day6_JAK3_normalized_paired_effects.csv"
  )
)

cat("\n===== NORMALIZED DONOR-WISE EFFECTS =====\n")
print(wide)

# ============================================================
# 8. DIRECTION CONSISTENCY SUMMARY
# ============================================================

effect_summary <- data.table(
  contrast = c(
    "Day6_LPS_vs_RPMI",
    "Day6_BCG_Denmark_vs_RPMI",
    "Day6_BCG_Bulgaria_vs_RPMI"
  ),

  mean_effect = c(
    mean(wide$LPS_minus_RPMI),
    mean(wide$Denmark_minus_RPMI),
    mean(wide$Bulgaria_minus_RPMI)
  ),

  positive_donors = c(
    sum(wide$LPS_minus_RPMI > 0),
    sum(wide$Denmark_minus_RPMI > 0),
    sum(wide$Bulgaria_minus_RPMI > 0)
  ),

  negative_donors = c(
    sum(wide$LPS_minus_RPMI < 0),
    sum(wide$Denmark_minus_RPMI < 0),
    sum(wide$Bulgaria_minus_RPMI < 0)
  )
)

effect_summary[
  ,
  direction_consistency :=
    pmax(
      positive_donors,
      negative_donors
    ) / 5
]

fwrite(
  effect_summary,
  file.path(
    outdir,
    "20c_GSE168468_Day6_JAK3_direction_summary.csv"
  )
)

cat("\n===== DIRECTION CONSISTENCY =====\n")
print(effect_summary)

# ============================================================
# 9. PAIRED DONOR FIGURE
# ============================================================

plotdat <- copy(
  as.data.table(jak3_norm)
)

plotdat$condition <- factor(
  plotdat$condition,
  levels = c(
    "RPMI",
    "LPS",
    "BCG_Denmark",
    "BCG_Bulgaria"
  ),
  labels = c(
    "RPMI",
    "LPS",
    "BCG Denmark",
    "BCG Bulgaria"
  )
)

p <- ggplot(
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
  labs(
    title =
      "GSE168468: JAK3 expression at Day 6 after stimulus removal",
    subtitle =
      "MMSEQ-normalized expression; paired measurements across five donors",
    x = NULL,
    y = "Normalized JAK3 log expression"
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
    "20c_GSE168468_Day6_JAK3_normalized_paired.png"
  ),
  p,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "20c_GSE168468_Day6_JAK3_normalized_paired.pdf"
  ),
  p,
  width = 7.5,
  height = 5.5
)

# ============================================================
# 10. MODEL METADATA
# ============================================================

metadata <- data.frame(
  dataset = "GSE168468",
  species = "Homo sapiens",
  cell_system =
    "Primary monocytes differentiated to macrophages",
  donors = 5,
  timepoint = "Day 6",
  initial_exposure = "24 h",
  recovery = "5 days after stimulus removal",
  expression_source = "MMSEQ gene-level log_mu",
  normalization =
    "Official MMSEQ readmmseq multi-sample normalization",
  statistical_model =
    "limma: normalized log expression ~ donor + condition",
  variance_moderation =
    "eBayes across genes; robust=TRUE",
  primary_control = "Day6 RPMI",
  stringsAsFactors = FALSE
)

write.csv(
  metadata,
  file.path(
    outdir,
    "20c_GSE168468_Day6_analysis_metadata.csv"
  ),
  row.names = FALSE
)

cat("\n20c COMPLETE\n")

