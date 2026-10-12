
rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

setDTthreads(1)

cat("\n====================================================\n")
cat("27a — GSE212282 JAK3 DAY6 MULTI-STIMULUS ANALYSIS\n")
cat("PAIRED DONOR MMSEQ + LIMMA\n")
cat("====================================================\n\n")

indir  <- "macrophage_transcriptomics_validation/GSE212282/data"
outdir <- "macrophage_transcriptomics_validation/GSE212282/results"
figdir <- "macrophage_transcriptomics_validation/GSE212282/figures"

helper <- "macrophage_transcriptomics_validation/GSE168468/scripts/mmseq.R"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(helper)) {
  stop("Missing MMSEQ helper: ", helper)
}

source(helper)

JAK3_ID <- "ENSG00000105639"

# ============================================================
# 1. SAMPLE METADATA
# ============================================================

samples <- data.frame(
  sample = c(
    "A_RPMI","A_BG","A_DMI","A_Ita","A_BG_DMI",
    "B_RPMI","B_BG","B_DMI","B_Ita","B_BG_DMI",
    "C_RPMI","C_BG","C_DMI","C_Ita","C_BG_DMI"
  ),

  donor = factor(rep(c("A","B","C"), each = 5)),

  condition = factor(
    rep(
      c("RPMI","BG","DMI","Itaconate","BG_DMI"),
      times = 3
    ),
    levels = c(
      "RPMI",
      "BG",
      "DMI",
      "Itaconate",
      "BG_DMI"
    )
  ),

  file = c(
    "GSM6523215_A_d6_RPMI.gene.mmseq.txt.gz",
    "GSM6523216_A_d6_BG.gene.mmseq.txt.gz",
    "GSM6523217_A_d6_DMI.gene.mmseq.txt.gz",
    "GSM6523218_A_d6_Ita.gene.mmseq.txt.gz",
    "GSM6523219_A_d6_BG_DMI.gene.mmseq.txt.gz",

    "GSM6523240_B_d6_RPMI.gene.mmseq.txt.gz",
    "GSM6523241_B_d6_BG.gene.mmseq.txt.gz",
    "GSM6523242_B_d6_DMI.gene.mmseq.txt.gz",
    "GSM6523243_B_d6_Ita.gene.mmseq.txt.gz",
    "GSM6523244_B_d6_BG_DMI.gene.mmseq.txt.gz",

    "GSM6523265_C_d6_RPMI.gene.mmseq.txt.gz",
    "GSM6523266_C_d6_BG.gene.mmseq.txt.gz",
    "GSM6523267_C_d6_DMI.gene.mmseq.txt.gz",
    "GSM6523268_C_d6_Ita.gene.mmseq.txt.gz",
    "GSM6523269_C_d6_BG_DMI.gene.mmseq.txt.gz"
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
cat("Conditions:", paste(levels(samples$condition), collapse = ", "), "\n")

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

if (!JAK3_ID %in% rownames(expr)) {
  stop("JAK3 not found.")
}

cat(
  "\nExpression matrix:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "samples\n"
)

# ============================================================
# 3. NORMALIZED JAK3 VALUES
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
    "27a_GSE212282_JAK3_normalized_Day6_values.csv"
  )
)

cat("\n===== NORMALIZED JAK3 VALUES =====\n")
print(jak3)

# ============================================================
# 4. FILTER GENES
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
# 5. DONOR-ADJUSTED MULTI-CONDITION MODEL
# ============================================================

donor <- factor(samples$donor)

condition <- factor(
  samples$condition,
  levels = c(
    "RPMI",
    "BG",
    "DMI",
    "Itaconate",
    "BG_DMI"
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

fit <- lmFit(
  expr_f,
  design
)

fit <- eBayes(
  fit,
  robust = TRUE
)

# ============================================================
# 6. FOUR PRIMARY CONTRASTS
#
# Because RPMI is reference level:
#
# conditionBG
# conditionDMI
# conditionItaconate
# conditionBG_DMI
# ============================================================

contrast_info <- data.table(
  contrast = c(
    "Day6_BG_vs_RPMI",
    "Day6_DMI_vs_RPMI",
    "Day6_Itaconate_vs_RPMI",
    "Day6_BG_DMI_vs_RPMI"
  ),

  coefficient = c(
    "conditionBG",
    "conditionDMI",
    "conditionItaconate",
    "conditionBG_DMI"
  )
)

all_results <- list()
jak3_results <- list()

for (i in seq_len(nrow(contrast_info))) {

  cname <- contrast_info$contrast[i]
  coefn <- contrast_info$coefficient[i]

  if (!coefn %in% colnames(design)) {
    stop("Coefficient missing: ", coefn)
  }

  tt <- topTable(
    fit,
    coef = coefn,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH",
    confint = TRUE
  )

  tt$gene_id <- rownames(tt)
  tt$contrast <- cname

  tt_dt <- as.data.table(tt)

  all_results[[cname]] <- tt_dt

  j <- tt_dt[
    gene_id == JAK3_ID
  ]

  if (nrow(j) != 1) {
    stop("Expected one JAK3 row for ", cname)
  }

  jak3_results[[cname]] <- data.table(
    dataset = "GSE212282",
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

  fwrite(
    tt_dt,
    file.path(
      outdir,
      paste0(
        "27a_GSE212282_",
        cname,
        "_genomewide_limma.csv"
      )
    )
  )
}

jak3_formal <- rbindlist(
  jak3_results
)

fwrite(
  jak3_formal,
  file.path(
    outdir,
    "27a_GSE212282_JAK3_formal_multicontrast_results.csv"
  )
)

cat("\n========================================\n")
cat("JAK3 FORMAL RESULTS\n")
cat("========================================\n")
print(jak3_formal)

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

wide[
  ,
  DMI_minus_RPMI :=
    DMI - RPMI
]

wide[
  ,
  Itaconate_minus_RPMI :=
    Itaconate - RPMI
]

wide[
  ,
  BG_DMI_minus_RPMI :=
    BG_DMI - RPMI
]

fwrite(
  wide,
  file.path(
    outdir,
    "27a_GSE212282_JAK3_donor_effects.csv"
  )
)

cat("\n===== DONOR-WISE JAK3 EFFECTS =====\n")
print(wide)

# ============================================================
# 8. DIRECTION CONSISTENCY
# ============================================================

effect_cols <- c(
  "BG_minus_RPMI",
  "DMI_minus_RPMI",
  "Itaconate_minus_RPMI",
  "BG_DMI_minus_RPMI"
)

pretty_names <- c(
  "BG_vs_RPMI",
  "DMI_vs_RPMI",
  "Itaconate_vs_RPMI",
  "BG_DMI_vs_RPMI"
)

direction_list <- list()

for (i in seq_along(effect_cols)) {

  x <- wide[[effect_cols[i]]]

  direction_list[[i]] <- data.table(
    contrast = pretty_names[i],
    mean_effect = mean(x),
    median_effect = median(x),
    positive_donors = sum(x > 0),
    negative_donors = sum(x < 0),
    zero_donors = sum(x == 0),
    direction_consistency =
      max(
        sum(x > 0),
        sum(x < 0)
      ) / length(x)
  )
}

direction <- rbindlist(
  direction_list
)

fwrite(
  direction,
  file.path(
    outdir,
    "27a_GSE212282_JAK3_direction_summary.csv"
  )
)

cat("\n===== DIRECTION CONSISTENCY =====\n")
print(direction)

# ============================================================
# 9. COMBINE FORMAL + DIRECTION SUMMARY
# ============================================================

formal2 <- copy(jak3_formal)

formal2[
  ,
  short_contrast := fifelse(
    contrast == "Day6_BG_vs_RPMI",
    "BG_vs_RPMI",
    fifelse(
      contrast == "Day6_DMI_vs_RPMI",
      "DMI_vs_RPMI",
      fifelse(
        contrast == "Day6_Itaconate_vs_RPMI",
        "Itaconate_vs_RPMI",
        "BG_DMI_vs_RPMI"
      )
    )
  )
]

final_summary <- merge(
  formal2,
  direction,
  by.x = "short_contrast",
  by.y = "contrast",
  all.x = TRUE
)

fwrite(
  final_summary,
  file.path(
    outdir,
    "27a_GSE212282_JAK3_final_summary.csv"
  )
)

# ============================================================
# 10. FIGURE — PAIRED VALUES
# ============================================================

plotdat <- copy(jak3)

plotdat$condition <- factor(
  plotdat$condition,
  levels = c(
    "RPMI",
    "BG",
    "DMI",
    "Itaconate",
    "BG_DMI"
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
  labs(
    title =
      "GSE212282: Day6 JAK3 expression across training conditions",
    subtitle =
      "Paired measurements across donors A, B and C after five-day differentiation/rest",
    x = NULL,
    y = "MMSEQ-normalized JAK3 log expression"
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "27a_GSE212282_JAK3_Day6_paired_expression.png"
  ),
  p1,
  width = 8,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "27a_GSE212282_JAK3_Day6_paired_expression.pdf"
  ),
  p1,
  width = 8,
  height = 5.5
)

# ============================================================
# 11. FIGURE — EFFECT ESTIMATES
# ============================================================

plot_effect <- copy(jak3_formal)

plot_effect[
  ,
  label := factor(
    contrast,
    levels = rev(c(
      "Day6_BG_vs_RPMI",
      "Day6_DMI_vs_RPMI",
      "Day6_Itaconate_vs_RPMI",
      "Day6_BG_DMI_vs_RPMI"
    ))
  )
]

p2 <- ggplot(
  plot_effect,
  aes(
    x = logFC,
    y = label
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
      yend = label
    ),
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    title =
      "GSE212282: Day6 JAK3 treatment effects",
    subtitle =
      "Donor-adjusted limma estimates with 95% confidence intervals",
    x =
      "Treatment minus RPMI normalized expression",
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "27a_GSE212282_JAK3_Day6_effect_sizes.png"
  ),
  p2,
  width = 8,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "27a_GSE212282_JAK3_Day6_effect_sizes.pdf"
  ),
  p2,
  width = 8,
  height = 5.5
)

# ============================================================
# 12. METADATA
# ============================================================

metadata <- data.frame(
  dataset = "GSE212282",
  species = "Homo sapiens",
  system =
    "Primary human monocytes differentiated to macrophages",
  donors =
    3,
  primary_timepoint =
    "Day6",
  rest_period =
    "5 days after initial exposure",
  conditions =
    "RPMI; BG; DMI; Itaconate; BG+DMI",
  expression =
    "Deposited MMSEQ gene-level output",
  normalization =
    "Joint readmmseq normalization",
  statistical_model =
    "~ donor + condition",
  inference =
    "Donor-adjusted limma moderated analysis",
  multiple_testing =
    "BH FDR separately within each genome-wide contrast",
  stringsAsFactors = FALSE
)

write.csv(
  metadata,
  file.path(
    outdir,
    "27a_GSE212282_analysis_metadata.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("27a COMPLETE\n")
cat("========================================\n")

