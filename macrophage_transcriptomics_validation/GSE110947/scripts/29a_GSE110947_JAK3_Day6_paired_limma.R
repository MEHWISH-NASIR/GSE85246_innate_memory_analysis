
rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

setDTthreads(1)

cat("\n====================================================\n")
cat("29a — GSE110947 JAK3 DAY6 VALIDATION\n")
cat("PAIRED 5-DONOR MMSEQ + LIMMA\n")
cat("====================================================\n\n")

indir  <- "macrophage_transcriptomics_validation/GSE110947/data"
outdir <- "macrophage_transcriptomics_validation/GSE110947/results"
figdir <- "macrophage_transcriptomics_validation/GSE110947/figures"

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
    "HD34_RPMI","HD34_Heme","HD34_BG",
    "HD37_RPMI","HD37_Heme","HD37_BG",
    "HD48_RPMI","HD48_Heme","HD48_BG",
    "HD49_RPMI","HD49_Heme","HD49_BG",
    "HD51_RPMI","HD51_Heme","HD51_BG"
  ),

  donor = factor(
    rep(c("HD34","HD37","HD48","HD49","HD51"), each = 3)
  ),

  condition = factor(
    rep(c("RPMI","Heme","BG"), times = 5),
    levels = c("RPMI","Heme","BG")
  ),

  file = c(
    "GSM3019186_bowtie_HD34_RPMI_d6_RNAseq_13578.gene.mmseq.txt.gz",
    "GSM3019187_bowtie_HD34_Heme_d6_RNAseq_13579.gene.mmseq.txt.gz",
    "GSM3019188_bowtie_HD34_BG_d6_RNAseq_13580.gene.mmseq.txt.gz",

    "GSM3019200_bowtie_HD37_RPMI_d6_RNAseq_13897.gene.mmseq.txt.gz",
    "GSM3019201_bowtie_HD37_Heme_d6_RNAseq_13898.gene.mmseq.txt.gz",
    "GSM3019202_bowtie_HD37_BG_d6_RNAseq_13899.gene.mmseq.txt.gz",

    "GSM3019214_bowtie_HD48_RPMI_d6_RNAseq_13795.gene.mmseq.txt.gz",
    "GSM3019215_bowtie_HD48_Heme_d6_RNAseq_13796.gene.mmseq.txt.gz",
    "GSM3019216_bowtie_HD48_BG_d6_RNAseq_13797.gene.mmseq.txt.gz",

    "GSM3019228_bowtie_HD49_RPMI_d6_RNAseq_13820.gene.mmseq.txt.gz",
    "GSM3019229_bowtie_HD49_Heme_d6_RNAseq_13821.gene.mmseq.txt.gz",
    "GSM3019230_bowtie_HD49_BG_d6_RNAseq_13822.gene.mmseq.txt.gz",

    "GSM3019248_bowtie_HD51_RPMI_d6_RNAseq_13883.gene.mmseq.txt.gz",
    "GSM3019249_bowtie_HD51_Heme_d6_RNAseq_13884.gene.mmseq.txt.gz",
    "GSM3019250_bowtie_HD51_BG_d6_RNAseq_13885.gene.mmseq.txt.gz"
  ),

  stringsAsFactors = FALSE
)

samples$path <- file.path(indir, samples$file)

missing <- samples$path[!file.exists(samples$path)]

if (length(missing) > 0) {
  stop(paste("Missing files:", paste(missing, collapse = "\n")))
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
  normalized_log_mu = as.numeric(expr[JAK3_ID, ])
)

fwrite(
  jak3,
  file.path(outdir, "29a_GSE110947_JAK3_normalized_Day6_values.csv")
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
) >= 5

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
# 4. PAIRED DONOR-ADJUSTED MODEL
# ============================================================

donor <- factor(samples$donor)

condition <- factor(
  samples$condition,
  levels = c("RPMI","Heme","BG")
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

fit <- lmFit(
  expr_f,
  design
)

fit <- eBayes(
  fit,
  robust = TRUE
)

# ============================================================
# 5. FORMAL CONTRASTS
# ============================================================

contrast_info <- data.table(
  contrast = c(
    "Day6_Heme_vs_RPMI",
    "Day6_BG_vs_RPMI"
  ),
  coefficient = c(
    "conditionHeme",
    "conditionBG"
  )
)

jak3_results <- list()

for (i in seq_len(nrow(contrast_info))) {

  cname <- contrast_info$contrast[i]
  coefn <- contrast_info$coefficient[i]

  tt <- topTable(
    fit,
    coef = coefn,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH",
    confint = TRUE
  )

  tt$gene_id <- rownames(tt)

  tt_dt <- as.data.table(tt)

  fwrite(
    tt_dt,
    file.path(
      outdir,
      paste0("29a_GSE110947_", cname, "_genomewide_limma.csv")
    )
  )

  j <- tt_dt[gene_id == JAK3_ID]

  if (nrow(j) != 1) {
    stop("Expected one JAK3 row for ", cname)
  }

  jak3_results[[cname]] <- data.table(
    dataset = "GSE110947",
    gene = "JAK3",
    gene_id = JAK3_ID,
    contrast = cname,
    logFC = j$logFC,
    CI_low = j$CI.L,
    CI_high = j$CI.R,
    t = j$t,
    P_value = j$P.Value,
    FDR = j$adj.P.Val
  )
}

jak3_formal <- rbindlist(jak3_results)

fwrite(
  jak3_formal,
  file.path(outdir, "29a_GSE110947_JAK3_formal_results.csv")
)

cat("\n========================================\n")
cat("JAK3 FORMAL RESULTS\n")
cat("========================================\n")
print(jak3_formal)

# ============================================================
# 6. DONOR-WISE EFFECTS
# ============================================================

wide <- dcast(
  jak3,
  donor ~ condition,
  value.var = "normalized_log_mu"
)

wide[, Heme_minus_RPMI := Heme - RPMI]
wide[, BG_minus_RPMI := BG - RPMI]

fwrite(
  wide,
  file.path(outdir, "29a_GSE110947_JAK3_donor_effects.csv")
)

cat("\n===== DONOR-WISE JAK3 EFFECTS =====\n")
print(wide)

# ============================================================
# 7. DIRECTION CONSISTENCY
# ============================================================

direction <- rbind(
  data.table(
    contrast = "Heme_vs_RPMI",
    mean_effect = mean(wide$Heme_minus_RPMI),
    median_effect = median(wide$Heme_minus_RPMI),
    positive_donors = sum(wide$Heme_minus_RPMI > 0),
    negative_donors = sum(wide$Heme_minus_RPMI < 0),
    direction_consistency =
      max(
        sum(wide$Heme_minus_RPMI > 0),
        sum(wide$Heme_minus_RPMI < 0)
      ) / nrow(wide)
  ),
  data.table(
    contrast = "BG_vs_RPMI",
    mean_effect = mean(wide$BG_minus_RPMI),
    median_effect = median(wide$BG_minus_RPMI),
    positive_donors = sum(wide$BG_minus_RPMI > 0),
    negative_donors = sum(wide$BG_minus_RPMI < 0),
    direction_consistency =
      max(
        sum(wide$BG_minus_RPMI > 0),
        sum(wide$BG_minus_RPMI < 0)
      ) / nrow(wide)
  )
)

fwrite(
  direction,
  file.path(outdir, "29a_GSE110947_JAK3_direction_summary.csv")
)

cat("\n===== DIRECTION CONSISTENCY =====\n")
print(direction)

# ============================================================
# 8. FIGURE — PAIRED EXPRESSION
# ============================================================

plotdat <- copy(jak3)

plotdat$condition <- factor(
  plotdat$condition,
  levels = c("RPMI","Heme","BG")
)

p1 <- ggplot(
  plotdat,
  aes(
    x = condition,
    y = normalized_log_mu,
    group = donor
  )
) +
  geom_line(alpha = 0.65) +
  geom_point(size = 3) +
  labs(
    title = "GSE110947: Day6 JAK3 expression",
    subtitle = "Paired resting macrophage measurements across five healthy donors",
    x = NULL,
    y = "MMSEQ-normalized JAK3 log expression"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(figdir, "29a_GSE110947_JAK3_Day6_paired_expression.png"),
  p1,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figdir, "29a_GSE110947_JAK3_Day6_paired_expression.pdf"),
  p1,
  width = 7,
  height = 5
)

# ============================================================
# 9. FIGURE — EFFECT ESTIMATES
# ============================================================

p2 <- ggplot(
  jak3_formal,
  aes(
    x = logFC,
    y = contrast
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
      yend = contrast
    ),
    linewidth = 0.8
  ) +
  geom_point(size = 3) +
  labs(
    title = "GSE110947: JAK3 Day6 treatment effects",
    subtitle = "Paired donor-adjusted limma estimates with 95% confidence intervals",
    x = "Treatment minus RPMI normalized expression",
    y = NULL
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(figdir, "29a_GSE110947_JAK3_Day6_effect_sizes.png"),
  p2,
  width = 7,
  height = 4.5,
  dpi = 300
)

ggsave(
  file.path(figdir, "29a_GSE110947_JAK3_Day6_effect_sizes.pdf"),
  p2,
  width = 7,
  height = 4.5
)

# ============================================================
# 10. METADATA
# ============================================================

metadata <- data.frame(
  dataset = "GSE110947",
  species = "Homo sapiens",
  system = "Primary human monocyte-derived macrophages",
  donors = 5,
  initial_conditions = "RPMI; Heme; Beta-glucan",
  initial_exposure = "24 h",
  rest_period = "5 days",
  analysis_state = "Day6 resting macrophages, before secondary LPS",
  expression = "Deposited MMSEQ gene-level output",
  normalization = "Joint readmmseq normalization",
  statistical_model = "~ donor + condition",
  inference = "Paired limma moderated analysis",
  stringsAsFactors = FALSE
)

write.csv(
  metadata,
  file.path(outdir, "29a_GSE110947_analysis_metadata.csv"),
  row.names = FALSE
)

cat("\n========================================\n")
cat("29a COMPLETE\n")
cat("========================================\n")

