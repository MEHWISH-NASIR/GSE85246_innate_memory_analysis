# ============================================================
# 02_initial_LPS_response.R
#
# GSE85246 / GSE85243
#
# INITIAL LPS TRANSCRIPTIONAL RESPONSE
#
# Biological question:
#
#   What transcriptional changes are associated with the
#   initial 24-h LPS exposure relative to RPMI control?
#
# Primary contrast:
#
#   Day 1 LPS - Day 1 RPMI
#
# Design:
#
#   Pair 1:
#     LPS_d1_6094  vs  RPMI_d1_6092
#
#   Pair 2:
#     LPS_d1_6716  vs  RPMI_d1_6715
#
# Statistical framework:
#   - MMSEQ normalized expression estimates
#   - natural log -> log2 conversion
#   - paired donor design
#   - limma empirical-Bayes modelling
#
# Kinase annotation:
#   Canonical KinHub/OpenKinome list established in Step 01
#
# IMPORTANT:
#   This script analyses INITIAL response only.
#
#   It does NOT define persistent/memory kinases.
#   It does NOT analyse Day 6.
#   It does NOT analyse beta-glucan rescue.
#   It does NOT analyse restimulation.
#
# Next:
#   03_persistent_kinase_memory.R
# ============================================================


rm(list = ls())


# ============================================================
# 1. REQUIRED PACKAGES / FILES
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
    paste0(
      "Missing required R package(s): ",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}


if (!file.exists("scripts/mmseq.R")) {
  
  stop(
    "Missing scripts/mmseq.R"
  )
}


if (!file.exists(
  "results/01_setup/01_setup_object.rds"
)) {
  
  stop(
    paste0(
      "Step 01 output is missing.\n",
      "Run scripts/01_setup_and_data.R first."
    )
  )
}


# Official/local MMSEQ helper functions
source(
  "scripts/mmseq.R"
)


# ============================================================
# 2. OUTPUT DIRECTORIES
# ============================================================

outdir <- "results/02_initial_LPS"
figdir <- "figures/02_initial_LPS"

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
# 3. LOAD STEP-01 SETUP
# ============================================================

setup <- readRDS(
  "results/01_setup/01_setup_object.rds"
)

samples <- setup$samples

kinase_symbols <- setup$kinase_symbols


cat("\n========================================\n")
cat("02 — INITIAL LPS RESPONSE\n")
cat("========================================\n\n")


cat(
  "Canonical KinHub/OpenKinome symbols:",
  length(kinase_symbols),
  "\n"
)


# ============================================================
# 4. EXTRACT DAY-1 SAMPLE INFORMATION
# ============================================================

day1_manifest <- samples[
  samples$analysis_set ==
    "Day1_state",
  ,
  drop = FALSE
]


if (nrow(day1_manifest) != 4) {
  
  stop(
    paste0(
      "Expected 4 Day-1 state samples, found ",
      nrow(day1_manifest),
      "."
    )
  )
}


cat("\nDay-1 samples from Step 01:\n\n")

print(
  day1_manifest[
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
# 5. DEFINITIVE DAY-1 ANALYSIS ORDER
#
# Keep this exact order because the paired design below
# corresponds to these four samples.
# ============================================================

analysis_samples <- c(
  "LPS_d1_6094",
  "LPS_d1_6716",
  "RPMI_d1_6092",
  "RPMI_d1_6715"
)


idx <- match(
  analysis_samples,
  day1_manifest$sample
)


if (anyNA(idx)) {
  
  stop(
    "One or more required Day-1 samples are absent from manifest."
  )
}


day1_manifest <- day1_manifest[
  idx,
  ,
  drop = FALSE
]


stopifnot(
  identical(
    day1_manifest$sample,
    analysis_samples
  )
)


if (!all(
  file.exists(
    day1_manifest$file
  )
)) {
  
  stop(
    "One or more Day-1 MMSEQ files are missing."
  )
}


# ============================================================
# 6. READ + NORMALIZE DAY-1 MMSEQ DATA
#
# readmmseq() performs the MMSEQ expression scaling used
# by the existing project helper.
#
# Day 1 is normalized independently because this script
# asks only about the Day-1 LPS response.
# ============================================================

cat("\n========================================\n")
cat("READING + NORMALIZING DAY-1 MMSEQ DATA\n")
cat("========================================\n\n")


d1 <- readmmseq(
  
  mmseq_files =
    day1_manifest$file,
  
  sample_names =
    day1_manifest$sample,
  
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
  names(d1)
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
  nrow(d1$log_mu),
  "\n"
)


cat(
  "Samples in MMSEQ object:",
  ncol(d1$log_mu),
  "\n"
)


# ============================================================
# 8. CONVERT MMSEQ NATURAL-LOG EXPRESSION TO LOG2
# ============================================================

expr_d1 <- d1$log_mu / log(2)


colnames(expr_d1) <-
  day1_manifest$sample


# ============================================================
# 9. EXPRESSION FILTER
#
# Require at least one unique hit in at least 2/4 samples.
#
# This removes genes with essentially no sequencing support
# while retaining genes measurable in at least half of the
# Day-1 samples.
# ============================================================

keep_d1 <- rowSums(
  d1$unique_hits >= 1,
  na.rm = TRUE
) >= 2


cat(
  "\nDay-1 genes before filtering:",
  nrow(expr_d1),
  "\n"
)


cat(
  "Day-1 genes retained:",
  sum(keep_d1),
  "\n"
)


cat(
  "Day-1 genes removed:",
  sum(!keep_d1),
  "\n"
)


expr_d1_f <- expr_d1[
  keep_d1,
  ,
  drop = FALSE
]


# Remove any zero-variance genes
gene_variance <- apply(
  expr_d1_f,
  1,
  var,
  na.rm = TRUE
)


keep_variance <-
  is.finite(gene_variance) &
  gene_variance > 0


expr_d1_f <- expr_d1_f[
  keep_variance,
  ,
  drop = FALSE
]


cat(
  "Genes retained after zero-variance removal:",
  nrow(expr_d1_f),
  "\n"
)


# ============================================================
# 10. PAIRED EXPERIMENTAL DESIGN
#
# Sample order:
#
#   LPS_d1_6094
#   LPS_d1_6716
#   RPMI_d1_6092
#   RPMI_d1_6715
#
# Corresponding pairs:
#
#   6094 ↔ 6092
#   6716 ↔ 6715
# ============================================================

pair <- factor(
  c(
    "pair1",
    "pair2",
    "pair1",
    "pair2"
  )
)


condition <- factor(
  c(
    "LPS",
    "LPS",
    "RPMI",
    "RPMI"
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
  day1_manifest$sample


cat("\n========================================\n")
cat("DAY-1 PAIRED DESIGN MATRIX\n")
cat("========================================\n\n")


print(design)


if (qr(design)$rank != ncol(design)) {
  
  stop(
    "Day-1 design matrix is not full rank."
  )
}


# ============================================================
# 11. DIFFERENTIAL EXPRESSION
#
# Coefficient:
#
#   conditionLPS = LPS - RPMI
#
# Positive log2FC:
#   higher expression after initial LPS exposure
#
# Negative log2FC:
#   lower expression after initial LPS exposure
# ============================================================

fit <- limma::lmFit(
  expr_d1_f,
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


# ============================================================
# 12. CLEAN RESULT COLUMN NAMES
# ============================================================

names(results_all)[
  names(results_all) == "logFC"
] <- "log2FC_Day1_LPS_vs_RPMI"


names(results_all)[
  names(results_all) == "adj.P.Val"
] <- "FDR_genome"


results_all <- results_all[
  ,
  c(
    "ENSEMBL",
    "log2FC_Day1_LPS_vs_RPMI",
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
cat("ANNOTATING GENES\n")
cat("========================================\n\n")


ensembl_ids <-
  results_all$ENSEMBL


gene_symbols <- AnnotationDbi::mapIds(
  
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
    gene_symbols[
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
  "Genes with HGNC symbol:",
  sum(
    !is.na(
      results_all$SYMBOL
    )
  ),
  "\n"
)


# ============================================================
# 14. ADD DIRECTION + SIGNIFICANCE FLAGS
# ============================================================

results_all$Direction <- ifelse(
  results_all$log2FC_Day1_LPS_vs_RPMI > 0,
  "UP",
  ifelse(
    results_all$log2FC_Day1_LPS_vs_RPMI < 0,
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
    results_all$log2FC_Day1_LPS_vs_RPMI
  ) >= 1


# ============================================================
# 15. SAVE COMPLETE GENOME-WIDE RESULT
# ============================================================

write.csv(
  results_all,
  file.path(
    outdir,
    "02_Day1_LPS_vs_RPMI_all_genes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. EXTRACT CANONICAL KINHUB / OPENKINOME GENES
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


# Kinase-level multiple-testing correction
#
# This is separate from the genome-wide FDR.
# It answers:
#
#   Among the measured canonical kinases, which have
#   evidence of an initial LPS response?
#
kinase_results$FDR_kinase <- p.adjust(
  kinase_results$P.Value,
  method = "BH"
)


kinase_results$Kinase_FDR_lt_0.05 <-
  kinase_results$FDR_kinase < 0.05


kinase_results$Kinase_FDR_lt_0.10 <-
  kinase_results$FDR_kinase < 0.10


# Sort by raw P value
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
    "02_Day1_LPS_vs_RPMI_KinHub_kinases.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. INITIAL LPS RESPONSE SUMMARIES
# ============================================================

all_nominal <- results_all[
  results_all$P.Value < 0.05,
  ,
  drop = FALSE
]


all_fdr05 <- results_all[
  results_all$FDR_genome < 0.05,
  ,
  drop = FALSE
]


all_fdr05_fc1 <- results_all[
  results_all$FDR_genome < 0.05 &
    abs(
      results_all$log2FC_Day1_LPS_vs_RPMI
    ) >= 1,
  ,
  drop = FALSE
]


kinase_nominal <- kinase_results[
  kinase_results$P.Value < 0.05,
  ,
  drop = FALSE
]


kinase_fdr10 <- kinase_results[
  kinase_results$FDR_kinase < 0.10,
  ,
  drop = FALSE
]


kinase_fdr05 <- kinase_results[
  kinase_results$FDR_kinase < 0.05,
  ,
  drop = FALSE
]


kinase_genome_fdr05 <- kinase_results[
  kinase_results$FDR_genome < 0.05,
  ,
  drop = FALSE
]


# ============================================================
# 18. SAVE USEFUL INITIAL-RESPONSE SUBSETS
# ============================================================

write.csv(
  kinase_nominal,
  file.path(
    outdir,
    "02_Day1_nominal_kinase_responders.csv"
  ),
  row.names = FALSE
)


write.csv(
  kinase_fdr05,
  file.path(
    outdir,
    "02_Day1_kinase_FDR05_responders.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. TOP KINASE TABLE FOR QUICK REVIEW
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
    "log2FC_Day1_LPS_vs_RPMI",
    "P.Value",
    "FDR_kinase",
    "FDR_genome",
    "Direction"
  ),
  drop = FALSE
]


write.csv(
  top_kinases,
  file.path(
    outdir,
    "02_Day1_top25_kinases.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. SIMPLE KINASE VOLCANO FIGURE
#
# This is descriptive.
# Final integrated kinase figure will be produced later.
# ============================================================

plot_p <- pmax(
  kinase_results$P.Value,
  .Machine$double.xmin
)


png(
  file.path(
    figdir,
    "02_Day1_KinHub_kinase_volcano.png"
  ),
  width = 1800,
  height = 1400,
  res = 180
)


plot(
  kinase_results$log2FC_Day1_LPS_vs_RPMI,
  -log10(plot_p),
  pch = 16,
  cex = 0.75,
  xlab = "Day-1 LPS vs RPMI log2 fold change",
  ylab = "-log10(P value)",
  main = "Initial LPS response — canonical kinases"
)


abline(
  v = 0,
  lty = 2
)


abline(
  h = -log10(0.05),
  lty = 2
)


# Label top ten kinase signals only
label_n <- min(
  10,
  nrow(kinase_results)
)


label_idx <- seq_len(
  label_n
)


text(
  kinase_results$log2FC_Day1_LPS_vs_RPMI[
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
# 21. SAVE ANALYSIS OBJECT
# ============================================================

saveRDS(
  list(
    
    manifest =
      day1_manifest,
    
    expression_log2 =
      expr_d1_f,
    
    design =
      design,
    
    all_genes =
      results_all,
    
    kinase_results =
      kinase_results,
    
    kinase_nominal =
      kinase_nominal,
    
    kinase_fdr10 =
      kinase_fdr10,
    
    kinase_fdr05 =
      kinase_fdr05
    
  ),
  file.path(
    outdir,
    "02_initial_LPS_response.rds"
  )
)


# ============================================================
# 22. SAVE TEXT SUMMARY
# ============================================================

summary_lines <- c(
  
  "GSE85246 — Step 02 Initial LPS Response",
  
  "",
  
  paste(
    "Day-1 samples:",
    nrow(day1_manifest)
  ),
  
  paste(
    "Genes tested:",
    nrow(results_all)
  ),
  
  paste(
    "Canonical KinHub/OpenKinome genes tested:",
    nrow(kinase_results)
  ),
  
  "",
  
  paste(
    "Genome-wide nominal P < 0.05:",
    nrow(all_nominal)
  ),
  
  paste(
    "Genome-wide FDR < 0.05:",
    nrow(all_fdr05)
  ),
  
  paste(
    "Genome-wide FDR < 0.05 and |log2FC| >= 1:",
    nrow(all_fdr05_fc1)
  ),
  
  "",
  
  paste(
    "Kinases nominal P < 0.05:",
    nrow(kinase_nominal)
  ),
  
  paste(
    "Kinases FDR < 0.10:",
    nrow(kinase_fdr10)
  ),
  
  paste(
    "Kinases FDR < 0.05:",
    nrow(kinase_fdr05)
  ),
  
  paste(
    "Kinases also genome-wide FDR < 0.05:",
    nrow(kinase_genome_fdr05)
  ),
  
  "",
  
  paste(
    "Contrast:",
    "Day1 LPS - Day1 RPMI"
  ),
  
  "",
  
  paste(
    "IMPORTANT:",
    "This step characterizes initial response only."
  ),
  
  paste(
    "Initial response is NOT used by itself",
    "to define cellular memory."
  )
)


writeLines(
  summary_lines,
  file.path(
    outdir,
    "02_initial_LPS_response_summary.txt"
  )
)


# ============================================================
# 23. FINAL CONSOLE REPORT
# ============================================================

cat("\n========================================\n")
cat("DAY-1 INITIAL LPS RESULTS\n")
cat("========================================\n\n")


cat(
  "Genes tested:",
  nrow(results_all),
  "\n"
)


cat(
  "Genome-wide nominal P < 0.05:",
  nrow(all_nominal),
  "\n"
)


cat(
  "Genome-wide FDR < 0.05:",
  nrow(all_fdr05),
  "\n"
)


cat(
  "Genome-wide FDR < 0.05 & |log2FC| >= 1:",
  nrow(all_fdr05_fc1),
  "\n"
)


cat("\n----------------------------------------\n")
cat("CANONICAL KINASE RESPONSE\n")
cat("----------------------------------------\n\n")


cat(
  "KinHub/OpenKinome genes tested:",
  nrow(kinase_results),
  "\n"
)


cat(
  "Nominal kinase P < 0.05:",
  nrow(kinase_nominal),
  "\n"
)


cat(
  "Kinase-level FDR < 0.10:",
  nrow(kinase_fdr10),
  "\n"
)


cat(
  "Kinase-level FDR < 0.05:",
  nrow(kinase_fdr05),
  "\n"
)


cat(
  "Also genome-wide FDR < 0.05:",
  nrow(kinase_genome_fdr05),
  "\n"
)


cat("\nTop kinase signals:\n\n")


print(
  top_kinases,
  row.names = FALSE,
  digits = 4
)


cat("\n========================================\n")
cat("02 INITIAL LPS RESPONSE COMPLETE\n")
cat("========================================\n")


cat(
  "\nThis analysis does NOT yet define memory kinases.\n"
)


cat(
  "\nNext canonical step:\n"
)


cat(
  "scripts/03_persistent_kinase_memory.R\n"
)


cat("\n========================================\n")