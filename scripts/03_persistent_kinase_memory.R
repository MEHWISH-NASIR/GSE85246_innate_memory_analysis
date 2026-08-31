# ============================================================
# 03_persistent_kinase_memory.R
#
# GSE85246 / GSE85243
#
# DAY-6 POST-WASHOUT / MEMORY-ASSOCIATED RNA STATE
#
# Biological question:
#
#   Which transcriptional differences remain at Day 6
#   after the initial LPS exposure?
#
# Primary contrast:
#
#   Day6 LPS-conditioned macrophages
#               -
#   Day6 RPMI macrophages
#
# Design:
#   5 matched LPS/RPMI pairs
#
# Statistical framework:
#   - MMSEQ normalized expression estimates
#   - natural log -> log2 conversion
#   - paired design
#   - limma empirical-Bayes model
#
# Kinase definition:
#   canonical KinHub/OpenKinome reference from Step 01
#
# IMPORTANT:
#   Day-1 significance is NOT required here.
#
#   Step 03 independently defines the Day-6
#   memory-associated kinase candidates.
#
#   Day-1 and Day-6 evidence will be integrated later
#   in Step 05.
# ============================================================


rm(list = ls())


# ============================================================
# 1. REQUIRED PACKAGES / INPUTS
# ============================================================

required_packages <- c(
  "limma",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  
  stop(
    paste(
      "Missing required package(s):",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}


if (!file.exists("scripts/mmseq.R")) {
  stop("Missing scripts/mmseq.R")
}


if (!file.exists(
  "results/01_setup/01_setup_object.rds"
)) {
  
  stop(
    paste0(
      "Step 01 output missing.\n",
      "Run scripts/01_setup_and_data.R first."
    )
  )
}


if (!file.exists(
  "results/02_initial_LPS/02_initial_LPS_response.rds"
)) {
  
  stop(
    paste0(
      "Step 02 output missing.\n",
      "Run scripts/02_initial_LPS_response.R first."
    )
  )
}


source(
  "scripts/mmseq.R"
)


# ============================================================
# 2. OUTPUT DIRECTORIES
# ============================================================

outdir <- "results/03_memory_kinases"
figdir <- "figures/03_memory_kinases"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. LOAD CLEAN PROJECT SETUP
# ============================================================

setup <- readRDS(
  "results/01_setup/01_setup_object.rds"
)

samples <- setup$samples

kinase_symbols <- setup$kinase_symbols


cat("\n========================================\n")
cat("03 — DAY-6 MEMORY-ASSOCIATED RNA STATE\n")
cat("========================================\n\n")

cat(
  "Canonical KinHub/OpenKinome symbols:",
  length(kinase_symbols),
  "\n"
)


# ============================================================
# 4. DAY-6 STATE SAMPLES
# ============================================================

day6_manifest <- samples[
  samples$analysis_set ==
    "Day6_state",
  ,
  drop = FALSE
]


if (nrow(day6_manifest) != 10) {
  
  stop(
    paste(
      "Expected 10 Day-6 state samples; found",
      nrow(day6_manifest)
    )
  )
}


cat("\nDay-6 samples:\n\n")

print(
  day6_manifest[
    ,
    c(
      "sample",
      "GSM",
      "condition",
      "state_pair"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 5. DEFINITIVE ANALYSIS ORDER
# ============================================================

analysis_samples <- c(
  
  # LPS-conditioned
  "LPS_d6_6097",
  "LPS_d6_6719",
  "LPS_d6_8747",
  "LPS_d6_9244",
  "LPS_d6_10649",
  
  # RPMI controls
  "RPMI_d6_6095",
  "RPMI_d6_6718",
  "RPMI_d6_8746",
  "RPMI_d6_9243",
  "RPMI_d6_10647"
)


idx <- match(
  analysis_samples,
  day6_manifest$sample
)


if (anyNA(idx)) {
  
  stop(
    "One or more required Day-6 samples are absent."
  )
}


day6_manifest <- day6_manifest[
  idx,
  ,
  drop = FALSE
]


stopifnot(
  identical(
    day6_manifest$sample,
    analysis_samples
  )
)


if (!all(
  file.exists(
    day6_manifest$file
  )
)) {
  
  stop(
    "One or more Day-6 MMSEQ files are missing."
  )
}


# ============================================================
# 6. READ + NORMALIZE DAY-6 MMSEQ DATA
# ============================================================

cat("\n========================================\n")
cat("READING + NORMALIZING DAY-6 MMSEQ DATA\n")
cat("========================================\n\n")


d6 <- readmmseq(
  
  mmseq_files =
    day6_manifest$file,
  
  sample_names =
    day6_manifest$sample,
  
  normalize = TRUE
)


# ============================================================
# 7. CHECK MMSEQ OBJECT
# ============================================================

required_objects <- c(
  "log_mu",
  "unique_hits"
)


missing_objects <- setdiff(
  required_objects,
  names(d6)
)


if (length(missing_objects) > 0) {
  
  stop(
    paste(
      "readmmseq output missing:",
      paste(
        missing_objects,
        collapse = ", "
      )
    )
  )
}


cat(
  "\nGenes in MMSEQ object:",
  nrow(d6$log_mu),
  "\n"
)


cat(
  "Samples in MMSEQ object:",
  ncol(d6$log_mu),
  "\n"
)


# ============================================================
# 8. NATURAL LOG -> LOG2 EXPRESSION
# ============================================================

expr_d6 <- d6$log_mu / log(2)


colnames(expr_d6) <-
  day6_manifest$sample


# ============================================================
# 9. EXPRESSION FILTER
#
# Require >=1 unique hit in at least 2 samples.
# ============================================================

keep_d6 <- rowSums(
  d6$unique_hits >= 1,
  na.rm = TRUE
) >= 2


cat(
  "\nDay-6 genes before filtering:",
  nrow(expr_d6),
  "\n"
)


cat(
  "Day-6 genes retained:",
  sum(keep_d6),
  "\n"
)


cat(
  "Day-6 genes removed:",
  sum(!keep_d6),
  "\n"
)


expr_d6_f <- expr_d6[
  keep_d6,
  ,
  drop = FALSE
]


# ============================================================
# 10. REMOVE ZERO-VARIANCE GENES
# ============================================================

gene_variance <- apply(
  expr_d6_f,
  1,
  var,
  na.rm = TRUE
)


keep_variance <-
  is.finite(gene_variance) &
  gene_variance > 0


expr_d6_f <- expr_d6_f[
  keep_variance,
  ,
  drop = FALSE
]


cat(
  "Genes after zero-variance removal:",
  nrow(expr_d6_f),
  "\n"
)


# ============================================================
# 11. PAIRED DAY-6 DESIGN
#
# Order:
#
# LPS 609
# LPS 671
# LPS 874
# LPS 924
# LPS donor70
#
# RPMI 609
# RPMI 671
# RPMI 874
# RPMI 924
# RPMI donor70
# ============================================================

pair <- factor(
  c(
    "pair609",
    "pair671",
    "pair874",
    "pair924",
    "pair70",
    
    "pair609",
    "pair671",
    "pair874",
    "pair924",
    "pair70"
  )
)


condition <- factor(
  c(
    rep("LPS", 5),
    rep("RPMI", 5)
  ),
  levels = c(
    "RPMI",
    "LPS"
  )
)


design <- model.matrix(
  ~ pair + condition
)


rownames(design) <-
  day6_manifest$sample


cat("\n========================================\n")
cat("DAY-6 PAIRED DESIGN MATRIX\n")
cat("========================================\n\n")


print(design)


if (qr(design)$rank != ncol(design)) {
  
  stop(
    "Day-6 design matrix is not full rank."
  )
}


# ============================================================
# 12. LIMMA MODEL
#
# coefficient conditionLPS:
#
#   Day6 LPS - Day6 RPMI
#
# Positive:
#   persistent increase after prior LPS exposure
#
# Negative:
#   persistent decrease after prior LPS exposure
# ============================================================

fit <- limma::lmFit(
  expr_d6_f,
  design
)


fit <- limma::eBayes(
  fit,
  trend = TRUE,
  robust = TRUE
)


results_all <- limma::topTable(
  
  fit,
  
  coef = "conditionLPS",
  
  number = Inf,
  
  sort.by = "P",
  
  adjust.method = "BH"
)


results_all$ENSEMBL <-
  rownames(results_all)


names(results_all)[
  names(results_all) == "logFC"
] <- "log2FC_Day6_LPS_vs_RPMI"


names(results_all)[
  names(results_all) == "adj.P.Val"
] <- "FDR_genome"


results_all <- results_all[
  ,
  c(
    "ENSEMBL",
    "log2FC_Day6_LPS_vs_RPMI",
    "AveExpr",
    "t",
    "P.Value",
    "FDR_genome",
    "B"
  )
]


# ============================================================
# 13. ENSEMBL -> HGNC ANNOTATION
# ============================================================

cat("\n========================================\n")
cat("ANNOTATING DAY-6 GENES\n")
cat("========================================\n\n")


ensembl_ids <-
  results_all$ENSEMBL


symbols <- AnnotationDbi::mapIds(
  
  org.Hs.eg.db::org.Hs.eg.db,
  
  keys = ensembl_ids,
  
  keytype = "ENSEMBL",
  
  column = "SYMBOL",
  
  multiVals = "first"
)


gene_names <- AnnotationDbi::mapIds(
  
  org.Hs.eg.db::org.Hs.eg.db,
  
  keys = ensembl_ids,
  
  keytype = "ENSEMBL",
  
  column = "GENENAME",
  
  multiVals = "first"
)


results_all$SYMBOL <-
  unname(
    symbols[
      results_all$ENSEMBL
    ]
  )


results_all$GENENAME <-
  unname(
    gene_names[
      results_all$ENSEMBL
    ]
  )


cat(
  "Genes with HGNC symbols:",
  sum(
    !is.na(
      results_all$SYMBOL
    )
  ),
  "\n"
)


# ============================================================
# 14. DIRECTION / SIGNIFICANCE FLAGS
# ============================================================

results_all$Direction <- ifelse(
  
  results_all$log2FC_Day6_LPS_vs_RPMI > 0,
  
  "UP",
  
  ifelse(
    
    results_all$log2FC_Day6_LPS_vs_RPMI < 0,
    
    "DOWN",
    
    "UNCHANGED"
  )
)


results_all$Nominal_P_lt_0.05 <-
  results_all$P.Value < 0.05


results_all$Genome_FDR_lt_0.05 <-
  results_all$FDR_genome < 0.05


results_all$Abs_log2FC_ge_1 <-
  abs(
    results_all$log2FC_Day6_LPS_vs_RPMI
  ) >= 1


# ============================================================
# 15. PER-PAIR DIRECTIONAL CONSISTENCY
#
# Calculate LPS - RPMI within each of the five pairs.
# This is SUPPORTING evidence, not a second significance test.
# ============================================================

pair_differences <- cbind(
  
  pair609 =
    expr_d6_f[, "LPS_d6_6097"] -
    expr_d6_f[, "RPMI_d6_6095"],
  
  pair671 =
    expr_d6_f[, "LPS_d6_6719"] -
    expr_d6_f[, "RPMI_d6_6718"],
  
  pair874 =
    expr_d6_f[, "LPS_d6_8747"] -
    expr_d6_f[, "RPMI_d6_8746"],
  
  pair924 =
    expr_d6_f[, "LPS_d6_9244"] -
    expr_d6_f[, "RPMI_d6_9243"],
  
  pair70 =
    expr_d6_f[, "LPS_d6_10649"] -
    expr_d6_f[, "RPMI_d6_10647"]
)


n_pairs_up <- rowSums(
  pair_differences > 0,
  na.rm = TRUE
)


n_pairs_down <- rowSums(
  pair_differences < 0,
  na.rm = TRUE
)


results_all$D6_pairs_up <-
  n_pairs_up[
    results_all$ENSEMBL
  ]


results_all$D6_pairs_down <-
  n_pairs_down[
    results_all$ENSEMBL
  ]


results_all$D6_direction_consistent_pairs <-
  ifelse(
    
    results_all$Direction == "UP",
    
    results_all$D6_pairs_up,
    
    ifelse(
      
      results_all$Direction == "DOWN",
      
      results_all$D6_pairs_down,
      
      0
    )
  )


# ============================================================
# 16. SAVE COMPLETE DAY-6 GENOME-WIDE RESULT
# ============================================================

write.csv(
  results_all,
  file.path(
    outdir,
    "03_Day6_LPS_vs_RPMI_all_genes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. CANONICAL KINASE RESULTS
# ============================================================

results_all$Is_KinHub_Kinase <-
  !is.na(results_all$SYMBOL) &
  results_all$SYMBOL %in%
  kinase_symbols


kinase_results <- results_all[
  results_all$Is_KinHub_Kinase,
  ,
  drop = FALSE
]


kinase_results$FDR_kinase <- p.adjust(
  kinase_results$P.Value,
  method = "BH"
)


kinase_results$Kinase_FDR_lt_0.05 <-
  kinase_results$FDR_kinase < 0.05


kinase_results$Kinase_FDR_lt_0.10 <-
  kinase_results$FDR_kinase < 0.10


kinase_results <- kinase_results[
  order(
    kinase_results$P.Value
  ),
  ,
  drop = FALSE
]


cat(
  "\nCanonical KinHub/OpenKinome genes represented:",
  nrow(kinase_results),
  "\n"
)


write.csv(
  kinase_results,
  file.path(
    outdir,
    "03_Day6_KinHub_all_kinases.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. PRIMARY MEMORY-ASSOCIATED KINASE SET
#
# Primary statistical definition:
#
#       Kinase-level FDR < 0.05
#
# This follows the AP-1-style logic:
# memory-associated kinases are selected from the memory
# contrast itself.
#
# Day-1 response is NOT required here.
# ============================================================

memory_kinases_FDR05 <- kinase_results[
  kinase_results$FDR_kinase < 0.05,
  ,
  drop = FALSE
]


memory_kinases_FDR10 <- kinase_results[
  kinase_results$FDR_kinase < 0.10,
  ,
  drop = FALSE
]


memory_kinases_nominal <- kinase_results[
  kinase_results$P.Value < 0.05,
  ,
  drop = FALSE
]


write.csv(
  memory_kinases_FDR05,
  file.path(
    outdir,
    "03_memory_kinases_FDR05.csv"
  ),
  row.names = FALSE
)


write.csv(
  memory_kinases_FDR10,
  file.path(
    outdir,
    "03_memory_kinases_FDR10.csv"
  ),
  row.names = FALSE
)


write.csv(
  memory_kinases_nominal,
  file.path(
    outdir,
    "03_memory_kinases_nominal.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. HIGH-CONSISTENCY SUPPORTING SUBSET
#
# This is NOT the primary definition.
#
# It asks whether the direction of the Day-6 effect is seen
# in at least 4 of the 5 matched pairs.
# ============================================================

memory_consistent <- memory_kinases_FDR05[
  memory_kinases_FDR05$
    D6_direction_consistent_pairs >= 4,
  ,
  drop = FALSE
]


write.csv(
  memory_consistent,
  file.path(
    outdir,
    "03_memory_kinases_FDR05_direction_consistent.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. GENOME-WIDE SUMMARY
# ============================================================

genome_nominal <- results_all[
  results_all$P.Value < 0.05,
  ,
  drop = FALSE
]


genome_fdr05 <- results_all[
  results_all$FDR_genome < 0.05,
  ,
  drop = FALSE
]


genome_fdr05_fc1 <- results_all[
  results_all$FDR_genome < 0.05 &
    abs(
      results_all$log2FC_Day6_LPS_vs_RPMI
    ) >= 1,
  ,
  drop = FALSE
]


# ============================================================
# 21. TOP 25 DAY-6 KINASES
# ============================================================

top_n <- min(
  25,
  nrow(kinase_results)
)


top_kinases <- kinase_results[
  seq_len(top_n),
  c(
    "ENSEMBL",
    "SYMBOL",
    "GENENAME",
    "log2FC_Day6_LPS_vs_RPMI",
    "P.Value",
    "FDR_kinase",
    "FDR_genome",
    "Direction",
    "D6_direction_consistent_pairs"
  ),
  drop = FALSE
]


write.csv(
  top_kinases,
  file.path(
    outdir,
    "03_Day6_top25_kinases.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. DAY-6 KINASE VOLCANO
# ============================================================

plot_p <- pmax(
  kinase_results$P.Value,
  .Machine$double.xmin
)


png(
  file.path(
    figdir,
    "03_Day6_KinHub_kinase_volcano.png"
  ),
  width = 1800,
  height = 1400,
  res = 180
)


plot(
  kinase_results$log2FC_Day6_LPS_vs_RPMI,
  -log10(plot_p),
  pch = 16,
  cex = 0.75,
  xlab = "Day-6 LPS-conditioned vs RPMI log2 fold change",
  ylab = "-log10(P value)",
  main = "Day-6 memory-associated canonical kinase signals"
)


abline(
  v = 0,
  lty = 2
)


abline(
  h = -log10(0.05),
  lty = 2
)


label_n <- min(
  10,
  nrow(kinase_results)
)


label_idx <- seq_len(
  label_n
)


text(
  kinase_results$log2FC_Day6_LPS_vs_RPMI[
    label_idx
  ],
  -log10(
    plot_p[
      label_idx
    ]
  ),
  labels =
    kinase_results$SYMBOL[
      label_idx
    ],
  pos = 3,
  cex = 0.65
)


dev.off()


# ============================================================
# 23. SAVE ANALYSIS OBJECT
# ============================================================

saveRDS(
  list(
    
    manifest =
      day6_manifest,
    
    expression_log2 =
      expr_d6_f,
    
    design =
      design,
    
    pair_differences =
      pair_differences,
    
    all_genes =
      results_all,
    
    kinase_results =
      kinase_results,
    
    memory_kinases_FDR05 =
      memory_kinases_FDR05,
    
    memory_kinases_FDR10 =
      memory_kinases_FDR10,
    
    memory_consistent =
      memory_consistent
    
  ),
  file.path(
    outdir,
    "03_persistent_kinase_memory.rds"
  )
)


# ============================================================
# 24. TEXT SUMMARY
# ============================================================

summary_lines <- c(
  
  "GSE85246 — Step 03 Day-6 Memory-Associated RNA State",
  
  "",
  
  paste(
    "Day-6 samples:",
    nrow(day6_manifest)
  ),
  
  paste(
    "Matched pairs:",
    5
  ),
  
  paste(
    "Genes tested:",
    nrow(results_all)
  ),
  
  paste(
    "Canonical KinHub/OpenKinome kinases tested:",
    nrow(kinase_results)
  ),
  
  "",
  
  paste(
    "Genome-wide nominal P < 0.05:",
    nrow(genome_nominal)
  ),
  
  paste(
    "Genome-wide FDR < 0.05:",
    nrow(genome_fdr05)
  ),
  
  paste(
    "Genome-wide FDR < 0.05 and |log2FC| >= 1:",
    nrow(genome_fdr05_fc1)
  ),
  
  "",
  
  paste(
    "Kinases nominal P < 0.05:",
    nrow(memory_kinases_nominal)
  ),
  
  paste(
    "Kinase-level FDR < 0.10:",
    nrow(memory_kinases_FDR10)
  ),
  
  paste(
    "Kinase-level FDR < 0.05:",
    nrow(memory_kinases_FDR05)
  ),
  
  paste(
    "FDR < 0.05 and direction consistent in >=4/5 pairs:",
    nrow(memory_consistent)
  ),
  
  "",
  
  "Primary memory-kinase definition:",
  "Kinase-level FDR < 0.05 in Day6 LPS vs RPMI.",
  
  "",
  
  paste(
    "Day-1 significance is not required",
    "for Step-03 candidate selection."
  )
)


writeLines(
  summary_lines,
  file.path(
    outdir,
    "03_memory_kinase_summary.txt"
  )
)


# ============================================================
# 25. FINAL CONSOLE REPORT
# ============================================================

cat("\n========================================\n")
cat("DAY-6 MEMORY RNA RESULTS\n")
cat("========================================\n\n")


cat(
  "Genes tested:",
  nrow(results_all),
  "\n"
)


cat(
  "Genome-wide nominal P < 0.05:",
  nrow(genome_nominal),
  "\n"
)


cat(
  "Genome-wide FDR < 0.05:",
  nrow(genome_fdr05),
  "\n"
)


cat(
  "Genome-wide FDR < 0.05 & |log2FC| >= 1:",
  nrow(genome_fdr05_fc1),
  "\n"
)


cat("\n----------------------------------------\n")
cat("CANONICAL MEMORY-ASSOCIATED KINASES\n")
cat("----------------------------------------\n\n")


cat(
  "KinHub/OpenKinome genes tested:",
  nrow(kinase_results),
  "\n"
)


cat(
  "Nominal P < 0.05:",
  nrow(memory_kinases_nominal),
  "\n"
)


cat(
  "Kinase-level FDR < 0.10:",
  nrow(memory_kinases_FDR10),
  "\n"
)


cat(
  "Kinase-level FDR < 0.05:",
  nrow(memory_kinases_FDR05),
  "\n"
)


cat(
  "FDR < 0.05 + >=4/5 directional consistency:",
  nrow(memory_consistent),
  "\n"
)


cat("\nTop Day-6 kinase signals:\n\n")


print(
  top_kinases,
  row.names = FALSE,
  digits = 4
)


if (nrow(memory_kinases_FDR05) > 0) {
  
  cat(
    "\n========================================\n"
  )
  
  cat(
    "PRIMARY DAY-6 MEMORY KINASE CANDIDATES\n"
  )
  
  cat(
    "========================================\n\n"
  )
  
  
  print(
    memory_kinases_FDR05[
      ,
      c(
        "SYMBOL",
        "log2FC_Day6_LPS_vs_RPMI",
        "P.Value",
        "FDR_kinase",
        "FDR_genome",
        "Direction",
        "D6_direction_consistent_pairs"
      )
    ],
    row.names = FALSE,
    digits = 4
  )
  
} else {
  
  cat(
    "\nNo kinase passes kinase-level FDR < 0.05.\n"
  )
  
  cat(
    "Do NOT substitute nominal hits as final candidates yet.\n"
  )
}


cat("\n========================================\n")
cat("03 MEMORY KINASE ANALYSIS COMPLETE\n")
cat("========================================\n")


cat(
  "\nNext canonical step:\n"
)


cat(
  "04_BG_interference.R\n"
)


cat("\n========================================\n")