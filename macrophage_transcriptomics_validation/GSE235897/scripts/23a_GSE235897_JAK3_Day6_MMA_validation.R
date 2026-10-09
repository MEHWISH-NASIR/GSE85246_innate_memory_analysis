
rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

setDTthreads(1)

cat("\n====================================================\n")
cat("23a — GSE235897 JAK3 DAY6 MMA VALIDATION\n")
cat("5-DAY RESTED TRAINED MACROPHAGES VS DAY6 NAIVE\n")
cat("====================================================\n\n")

indir <- "macrophage_transcriptomics_validation/GSE235897/data/processed"
outdir <- "macrophage_transcriptomics_validation/GSE235897/results"
figdir <- "macrophage_transcriptomics_validation/GSE235897/figures"

helper <- "macrophage_transcriptomics_validation/GSE168468/scripts/mmseq.R"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(helper)) {
  stop("Missing MMSEQ helper: ", helper)
}

source(helper)

JAK3_ID <- "ENSG00000105639"

# ============================================================
# 1. PRIMARY SAMPLE MANIFEST
#
# Formal comparison:
# Day6 MMA (24 h training + 5-day rest)
# versus
# Day6 naive macrophages
#
# GEO reports three biological replicates, but does not
# explicitly provide donor IDs linking replicate numbers.
# Therefore formal inference is intentionally UNPAIRED.
# ============================================================

samples <- data.frame(
  sample = c(
    "Naive_1",
    "Naive_2",
    "Naive_3",
    "MMA_D6_1",
    "MMA_D6_2",
    "MMA_D6_3"
  ),

  replicate_index = c(
    "1", "2", "3",
    "1", "2", "3"
  ),

  condition = factor(
    c(
      "Naive",
      "Naive",
      "Naive",
      "MMA_D6",
      "MMA_D6",
      "MMA_D6"
    ),
    levels = c(
      "Naive",
      "MMA_D6"
    )
  ),

  file = c(
    "GSM7511666_human_day6_Naive_1_RNAseq.gene.mmseq.txt.gz",
    "GSM7511667_human_day6_Naive_2_RNAseq.gene.mmseq.txt.gz",
    "GSM7511668_human_day6_Naive_3_RNAseq.gene.mmseq.txt.gz",

    "GSM7511672_human_day6_MMA_1_RNAseq.gene.mmseq.txt.gz",
    "GSM7511673_human_day6_MMA_2_RNAseq.gene.mmseq.txt.gz",
    "GSM7511674_human_day6_MMA_3_RNAseq.gene.mmseq.txt.gz"
  ),

  stringsAsFactors = FALSE
)

samples$path <- file.path(
  indir,
  samples$file
)

missing <- samples$path[
  !file.exists(samples$path)
]

if (length(missing) > 0) {
  stop(
    paste(
      "Missing files:",
      paste(missing, collapse = "\n")
    )
  )
}

cat("Primary samples:", nrow(samples), "\n")
cat("Naive:", sum(samples$condition == "Naive"), "\n")
cat("MMA Day6:", sum(samples$condition == "MMA_D6"), "\n")

# ============================================================
# 2. JOINT MMSEQ NORMALIZATION
# ============================================================

cat("\nRunning joint MMSEQ normalization...\n")

mm <- readmmseq(
  mmseq_files = samples$path,
  sample_names = samples$sample,
  normalize = TRUE
)

expr <- mm$log_mu

if (ncol(expr) != 6) {
  stop("Expected six normalized samples.")
}

if (!JAK3_ID %in% rownames(expr)) {
  stop("JAK3 not found.")
}

cat(
  "\nNormalized matrix:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "samples\n"
)

# ============================================================
# 3. SAVE NORMALIZED JAK3 VALUES
# ============================================================

jak3 <- data.table(
  sample = samples$sample,
  replicate_index = samples$replicate_index,
  condition = samples$condition,
  normalized_log_mu =
    as.numeric(expr[JAK3_ID, ])
)

fwrite(
  jak3,
  file.path(
    outdir,
    "23a_GSE235897_JAK3_normalized_Day6_values.csv"
  )
)

cat("\n===== NORMALIZED JAK3 VALUES =====\n")
print(jak3)

# ============================================================
# 4. FILTER LOW-INFORMATION GENES
#
# Keep features observed in >=3 of 6 samples.
# ============================================================

observed <- mm$observed

keep <- rowSums(
  observed == 1,
  na.rm = TRUE
) >= 3

cat("\nGenes before filtering:", nrow(expr), "\n")
cat("Genes retained:", sum(keep), "\n")

if (!keep[JAK3_ID]) {
  stop("JAK3 failed observed-feature filter.")
}

expr_f <- expr[
  keep,
  ,
  drop = FALSE
]

# ============================================================
# 5. UNPAIRED TWO-GROUP LIMMA MODEL
# ============================================================

condition <- factor(
  samples$condition,
  levels = c(
    "Naive",
    "MMA_D6"
  )
)

design <- model.matrix(
  ~ 0 + condition
)

colnames(design) <- c(
  "Naive",
  "MMA_D6"
)

contrast_matrix <- makeContrasts(
  MMA_D6_vs_Naive =
    MMA_D6 - Naive,
  levels = design
)

cat("\n===== DESIGN =====\n")
print(design)

# ============================================================
# 6. FIT ALL RETAINED GENES
# ============================================================

fit <- lmFit(
  expr_f,
  design
)

fit <- contrasts.fit(
  fit,
  contrast_matrix
)

fit <- eBayes(
  fit,
  robust = TRUE
)

tt <- topTable(
  fit,
  coef = "MMA_D6_vs_Naive",
  number = Inf,
  sort.by = "none",
  adjust.method = "BH",
  confint = TRUE
)

tt$feature_id <- rownames(tt)
tt$contrast <- "MMA_D6_vs_Naive"

fwrite(
  as.data.table(tt),
  file.path(
    outdir,
    "23a_GSE235897_Day6_MMA_vs_Naive_genomewide_limma.csv"
  )
)

# ============================================================
# 7. JAK3 FORMAL RESULT
# ============================================================

j <- tt[
  tt$feature_id == JAK3_ID,
  ,
  drop = FALSE
]

if (nrow(j) != 1) {
  stop("Expected exactly one JAK3 result.")
}

jak3_result <- data.table(
  dataset = "GSE235897",
  gene = "JAK3",
  gene_id = JAK3_ID,
  contrast = "Day6_MMA_vs_Day6_Naive",
  logFC = j$logFC,
  CI_low = j$CI.L,
  CI_high = j$CI.R,
  t = j$t,
  P_value = j$P.Value,
  FDR = j$adj.P.Val
)

fwrite(
  jak3_result,
  file.path(
    outdir,
    "23a_GSE235897_JAK3_Day6_limma_result.csv"
  )
)

cat("\n========================================\n")
cat("JAK3 FORMAL RESULT\n")
cat("========================================\n")
print(jak3_result)

# ============================================================
# 8. DESCRIPTIVE REPLICATE-INDEX DIFFERENCES
#
# IMPORTANT:
# Matching 1->1, 2->2, 3->3 is descriptive only.
# GEO does not explicitly identify these as matched donors.
# ============================================================

wide <- dcast(
  jak3,
  replicate_index ~ condition,
  value.var = "normalized_log_mu"
)

wide[
  ,
  MMA_minus_Naive :=
    MMA_D6 - Naive
]

fwrite(
  wide,
  file.path(
    outdir,
    "23a_GSE235897_JAK3_replicate_index_differences.csv"
  )
)

cat("\n===== DESCRIPTIVE REPLICATE-INDEX DIFFERENCES =====\n")
print(wide)

direction <- data.table(
  positive = sum(
    wide$MMA_minus_Naive > 0
  ),

  negative = sum(
    wide$MMA_minus_Naive < 0
  ),

  mean_difference = mean(
    wide$MMA_minus_Naive
  ),

  median_difference = median(
    wide$MMA_minus_Naive
  )
)

direction[
  ,
  direction_consistency :=
    max(
      positive,
      negative
    ) / 3
]

fwrite(
  direction,
  file.path(
    outdir,
    "23a_GSE235897_JAK3_direction_summary.csv"
  )
)

cat("\n===== DESCRIPTIVE DIRECTION SUMMARY =====\n")
print(direction)

# ============================================================
# 9. FIGURE 1 — ALL SIX REPLICATES
# ============================================================

plotdat <- copy(jak3)

plotdat$condition <- factor(
  plotdat$condition,
  levels = c(
    "Naive",
    "MMA_D6"
  ),
  labels = c(
    "Day6 Naive",
    "Day6 MMA\n(5-day rest)"
  )
)

p1 <- ggplot(
  plotdat,
  aes(
    x = condition,
    y = normalized_log_mu
  )
) +
  geom_jitter(
    width = 0.07,
    height = 0,
    size = 3
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.35
  ) +
  labs(
    title =
      "GSE235897: JAK3 after five-day trained-immunity rest",
    subtitle =
      "Three biological replicates per condition",
    x = NULL,
    y = "MMSEQ-normalized JAK3 log expression"
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "23a_GSE235897_JAK3_Day6_expression.png"
  ),
  p1,
  width = 6.5,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "23a_GSE235897_JAK3_Day6_expression.pdf"
  ),
  p1,
  width = 6.5,
  height = 5
)

# ============================================================
# 10. FIGURE 2 — EFFECT + 95% CI
# ============================================================

forest <- copy(jak3_result)

p2 <- ggplot(
  forest,
  aes(
    x = logFC,
    y = "Day6 MMA vs Naive"
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_segment(
    aes(
      x = CI_low,
      xend = CI_high,
      y = "Day6 MMA vs Naive",
      yend = "Day6 MMA vs Naive"
    ),
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    title =
      "GSE235897: JAK3 Day6 treatment effect",
    subtitle =
      "Limma moderated estimate with 95% confidence interval",
    x =
      "MMA Day6 minus naive Day6",
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "23a_GSE235897_JAK3_Day6_effect_size.png"
  ),
  p2,
  width = 7,
  height = 4,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "23a_GSE235897_JAK3_Day6_effect_size.pdf"
  ),
  p2,
  width = 7,
  height = 4
)

# ============================================================
# 11. ANALYSIS METADATA
# ============================================================

metadata <- data.frame(
  dataset = "GSE235897",
  species = "Homo sapiens",
  cell_system =
    "Human monocyte-derived macrophages",
  biological_replicates = 3,
  training =
    "Al(OH)3 + PHAD + Mannan (MMA)",
  training_duration =
    "24 h",
  rest_period =
    "5 days",
  primary_comparison =
    "Day6 MMA vs Day6 naive",
  expression_source =
    "Deposited gene-level MMSEQ output",
  normalization =
    "MMSEQ readmmseq joint normalization",
  formal_model =
    "Unpaired two-group limma model",
  pairing_assumption =
    "None for formal inference",
  replicate_index_matching =
    "Descriptive only; GEO does not explicitly provide matched donor IDs",
  stringsAsFactors = FALSE
)

write.csv(
  metadata,
  file.path(
    outdir,
    "23a_GSE235897_analysis_metadata.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("23a COMPLETE\n")
cat("========================================\n")

