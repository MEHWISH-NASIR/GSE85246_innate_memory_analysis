
rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

setDTthreads(1)

cat("\n====================================================\n")
cat("26a — GSE141656 JAK3 BETA-GLUCAN VALIDATION\n")
cat("DAY6 POST-REST PAIRED LIMMA ANALYSIS\n")
cat("====================================================\n\n")

indir <- "macrophage_transcriptomics_validation/GSE141656/data"
outdir <- "macrophage_transcriptomics_validation/GSE141656/results"
figdir <- "macrophage_transcriptomics_validation/GSE141656/figures"

helper <- "macrophage_transcriptomics_validation/GSE168468/scripts/mmseq.R"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(helper)) {
  stop("Missing MMSEQ helper: ", helper)
}

source(helper)

JAK3_ID <- "ENSG00000105639"

samples <- data.frame(
  sample = c(
    "A_RPMI",
    "A_BG",
    "B_RPMI",
    "B_BG",
    "C_RPMI",
    "C_BG"
  ),

  donor = factor(
    c(
      "A", "A",
      "B", "B",
      "C", "C"
    )
  ),

  condition = factor(
    c(
      "RPMI", "BG",
      "RPMI", "BG",
      "RPMI", "BG"
    ),
    levels = c("RPMI", "BG")
  ),

  file = c(
    "GSM4210636_A_RPMI_RPMI.gene.mmseq.txt.gz",
    "GSM4210634_A_BG_RPMI.gene.mmseq.txt.gz",
    "GSM4210640_B_RPMI_RPMI.gene.mmseq.txt.gz",
    "GSM4210638_B_BG_RPMI.gene.mmseq.txt.gz",
    "GSM4210644_C_RPMI_RPMI.gene.mmseq.txt.gz",
    "GSM4210642_C_BG_RPMI.gene.mmseq.txt.gz"
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

cat("Samples:", nrow(samples), "\n")
cat("Donors:", length(unique(samples$donor)), "\n")

# ============================================================
# 1. JOINT MMSEQ NORMALIZATION
# ============================================================

cat("\nRunning joint MMSEQ normalization...\n")

mm <- readmmseq(
  mmseq_files = samples$path,
  sample_names = samples$sample,
  normalize = TRUE
)

expr <- mm$log_mu

if (!JAK3_ID %in% rownames(expr)) {
  stop("JAK3 not found.")
}

cat(
  "Expression matrix:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "samples\n"
)

# ============================================================
# 2. NORMALIZED JAK3 VALUES
# ============================================================

jak3 <- data.table(
  sample = samples$sample,
  donor = samples$donor,
  condition = samples$condition,
  normalized_log_mu =
    as.numeric(expr[JAK3_ID, ])
)

fwrite(
  jak3,
  file.path(
    outdir,
    "26a_GSE141656_JAK3_normalized_values.csv"
  )
)

cat("\n===== NORMALIZED JAK3 VALUES =====\n")
print(jak3)

# ============================================================
# 3. FILTER LOW-INFORMATION GENES
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
# 4. PAIRED LIMMA MODEL
# ============================================================

donor <- factor(samples$donor)

condition <- factor(
  samples$condition,
  levels = c(
    "RPMI",
    "BG"
  )
)

design <- model.matrix(
  ~ donor + condition
)

cat("\n===== DESIGN MATRIX =====\n")
print(design)

cat("\nCoefficients:\n")
print(colnames(design))

if (qr(design)$rank != ncol(design)) {
  stop("Design matrix is not full rank.")
}

# ============================================================
# 5. GENOME-WIDE FIT
# ============================================================

fit <- lmFit(
  expr_f,
  design
)

fit <- eBayes(
  fit,
  robust = TRUE
)

coef_name <- grep(
  "^conditionBG$",
  colnames(design),
  value = TRUE
)

if (length(coef_name) != 1) {
  stop("Could not identify BG coefficient.")
}

tt <- topTable(
  fit,
  coef = coef_name,
  number = Inf,
  sort.by = "none",
  adjust.method = "BH",
  confint = TRUE
)

tt$feature_id <- rownames(tt)

fwrite(
  as.data.table(tt),
  file.path(
    outdir,
    "26a_GSE141656_BG_vs_RPMI_genomewide_limma.csv"
  )
)

# ============================================================
# 6. JAK3 FORMAL RESULT
# ============================================================

j <- tt[
  tt$feature_id == JAK3_ID,
  ,
  drop = FALSE
]

if (nrow(j) != 1) {
  stop("Expected exactly one JAK3 row.")
}

jak3_result <- data.table(
  dataset = "GSE141656",
  gene = "JAK3",
  gene_id = JAK3_ID,
  comparison = "Day6_BG_vs_Day6_RPMI",
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
    "26a_GSE141656_JAK3_paired_limma_result.csv"
  )
)

cat("\n========================================\n")
cat("JAK3 FORMAL RESULT\n")
cat("========================================\n")
print(jak3_result)

# ============================================================
# 7. DONOR-WISE EFFECTS
# ============================================================

wide <- dcast(
  jak3,
  donor ~ condition,
  value.var = "normalized_log_mu"
)

wide[
  ,
  BG_minus_RPMI :=
    BG - RPMI
]

fwrite(
  wide,
  file.path(
    outdir,
    "26a_GSE141656_JAK3_donor_effects.csv"
  )
)

cat("\n===== DONOR-WISE JAK3 EFFECTS =====\n")
print(wide)

# ============================================================
# 8. DIRECTION CONSISTENCY
# ============================================================

direction <- data.table(
  mean_effect =
    mean(
      wide$BG_minus_RPMI
    ),

  median_effect =
    median(
      wide$BG_minus_RPMI
    ),

  positive_donors =
    sum(
      wide$BG_minus_RPMI > 0
    ),

  negative_donors =
    sum(
      wide$BG_minus_RPMI < 0
    )
)

direction[
  ,
  direction_consistency :=
    max(
      positive_donors,
      negative_donors
    ) /
    nrow(wide)
]

fwrite(
  direction,
  file.path(
    outdir,
    "26a_GSE141656_JAK3_direction_summary.csv"
  )
)

cat("\n===== DIRECTION CONSISTENCY =====\n")
print(direction)

# ============================================================
# 9. FIGURE 1
# ============================================================

p1 <- ggplot(
  jak3,
  aes(
    x = condition,
    y = normalized_log_mu,
    group = donor
  )
) +
  geom_line(
    alpha = 0.7
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    title =
      "GSE141656: JAK3 after beta-glucan training and five-day rest",
    subtitle =
      "Paired Day6 measurements across donors A, B and C",
    x = NULL,
    y = "MMSEQ-normalized JAK3 log expression"
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "26a_GSE141656_JAK3_paired_expression.png"
  ),
  p1,
  width = 6.5,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "26a_GSE141656_JAK3_paired_expression.pdf"
  ),
  p1,
  width = 6.5,
  height = 5
)

# ============================================================
# 10. FIGURE 2
# ============================================================

p2 <- ggplot(
  jak3_result,
  aes(
    x = logFC,
    y = "Beta-glucan vs RPMI"
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
      y = "Beta-glucan vs RPMI",
      yend = "Beta-glucan vs RPMI"
    ),
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    title =
      "GSE141656: JAK3 Day6 beta-glucan effect",
    subtitle =
      "Paired limma estimate with 95% confidence interval",
    x =
      "Beta-glucan minus RPMI normalized expression",
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "26a_GSE141656_JAK3_effect_size.png"
  ),
  p2,
  width = 7,
  height = 4,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "26a_GSE141656_JAK3_effect_size.pdf"
  ),
  p2,
  width = 7,
  height = 4
)

# ============================================================
# 11. METADATA
# ============================================================

metadata <- data.frame(
  dataset = "GSE141656",
  species = "Homo sapiens",
  system =
    "Human monocyte-derived macrophages",
  training =
    "Beta-glucan",
  training_duration =
    "24 h",
  rest_period =
    "5 days",
  analysis_state =
    "Day6, no Mtb restimulation",
  donors =
    3,
  expression_source =
    "Deposited MMSEQ gene-level output",
  normalization =
    "Joint MMSEQ readmmseq normalization",
  formal_model =
    "~ donor + condition",
  inference =
    "Paired limma moderated analysis",
  stringsAsFactors = FALSE
)

write.csv(
  metadata,
  file.path(
    outdir,
    "26a_GSE141656_analysis_metadata.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("26a COMPLETE\n")
cat("========================================\n")

