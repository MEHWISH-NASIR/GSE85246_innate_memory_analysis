# ============================================================
# 12k_allpromoter_ebayes_kinaseBH.R
#
# PURPOSE
# Test whether Tobias's analysis used:
#
#   all measurable promoters
#       -> corrected RPMI_rep2
#       -> log2(signal + 1)
#       -> unpaired limma
#       -> eBayes across ALL promoters
#       -> extract kinase genes
#       -> BH correction within kinase genes only
#
# This differs from 12j, where kinase genes were subset
# BEFORE lmFit/eBayes.
# ============================================================

rm(list = ls())

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Package 'limma' is required.")
}

outdir <- "results/12_kinase_reevaluation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

raw_file <- file.path(
  outdir,
  "12d_day6_promoter_signal_raw.csv"
)

scale_file <- file.path(
  outdir,
  "12e_RPMI_rep2_scaling_factors.csv"
)

kinase_file <- paste0(
  "results/03_memory_kinases/",
  "03_Day6_KinHub_all_kinases.csv"
)

for (f in c(raw_file, scale_file, kinase_file)) {
  if (!file.exists(f)) {
    stop(paste("Missing required file:", f))
  }
}

cat("\n========================================\n")
cat("12k — ALL-PROMOTER eBAYES + KINASE BH\n")
cat("========================================\n\n")


# ============================================================
# 1. LOAD FULL PROMOTER MATRIX
# ============================================================

pm <- read.csv(
  raw_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Promoters loaded:", nrow(pm), "\n")


# ============================================================
# 2. LOAD SCALING FACTORS
# ============================================================

sc <- read.csv(
  scale_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_scale_cols <- c(
  "mark",
  "primary_scaling_factor"
)

if (!all(required_scale_cols %in% names(sc))) {
  stop("Scaling table missing required columns.")
}

cat("\nPrimary peer-median scaling factors:\n")
print(
  sc[, required_scale_cols],
  row.names = FALSE
)


# ============================================================
# 3. APPLY PEER-MEDIAN CORRECTION
#    ONLY TO RPMI_rep2 FOR EACH MARK
# ============================================================

marks <- c(
  "H3K27ac",
  "H3K4me1",
  "H3K4me3"
)

pm_corrected <- pm

for (mark in marks) {

  col <- paste0(mark, "_RPMI_rep2")

  if (!col %in% names(pm_corrected)) {
    stop(paste("Missing signal column:", col))
  }

  sf <- sc$primary_scaling_factor[
    sc$mark == mark
  ]

  if (length(sf) != 1 || !is.finite(sf)) {
    stop(
      paste(
        "Could not resolve scaling factor for",
        mark
      )
    )
  }

  cat(
    mark,
    "RPMI_rep2 scaling factor =",
    sf,
    "\n"
  )

  pm_corrected[[col]] <-
    pm_corrected[[col]] * sf
}


# ============================================================
# 4. KINASE UNIVERSE
# ============================================================

kin <- read.csv(
  kinase_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kinase_symbols <- unique(
  kin$SYMBOL[
    !is.na(kin$SYMBOL) &
      kin$SYMBOL != "" &
      kin$Is_KinHub_Kinase
  ]
)

kinase_idx <-
  !is.na(pm_corrected$SYMBOL) &
  pm_corrected$SYMBOL %in% kinase_symbols

cat("\nKinHub symbols defined:", length(kinase_symbols), "\n")
cat(
  "Kinase promoters measurable:",
  sum(kinase_idx),
  "\n"
)


# ============================================================
# 5. UNPAIRED DESIGN
# ============================================================

condition <- factor(
  c(
    "BG", "BG",
    "LPS", "LPS",
    "RPMI", "RPMI"
  ),
  levels = c(
    "BG",
    "LPS",
    "RPMI"
  )
)

design <- model.matrix(
  ~ 0 + condition
)

colnames(design) <- levels(condition)

contrast_matrix <- limma::makeContrasts(
  LPS_vs_RPMI = LPS - RPMI,
  BG_vs_RPMI  = BG  - RPMI,
  levels = design
)

cat("\nDesign matrix:\n")
print(design)

cat("\nContrasts:\n")
print(contrast_matrix)


# ============================================================
# 6. FIT EACH MARK USING ALL PROMOTERS
# ============================================================

all_kinase_results <- list()

model_summary <- list()

for (mark in marks) {

  cat("\n========================================\n")
  cat("MARK:", mark, "\n")
  cat("========================================\n")

  signal_cols <- paste0(
    mark,
    c(
      "_BG_rep1",
      "_BG_rep2",
      "_LPS_rep1",
      "_LPS_rep2",
      "_RPMI_rep1",
      "_RPMI_rep2"
    )
  )

  missing_cols <- setdiff(
    signal_cols,
    names(pm_corrected)
  )

  if (length(missing_cols) > 0) {
    stop(
      paste(
        "Missing columns:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  mat_raw <- as.matrix(
    pm_corrected[, signal_cols, drop = FALSE]
  )

  storage.mode(mat_raw) <- "numeric"

  if (anyNA(mat_raw)) {
    stop(
      paste(
        "NA signal encountered for",
        mark
      )
    )
  }

  if (any(mat_raw < 0)) {
    stop(
      paste(
        "Negative signal encountered for",
        mark
      )
    )
  }

  # Tobias: log transformed.
  # Keep same transformation as 12j for this diagnostic:
  mat <- log2(mat_raw + 1)

  rownames(mat) <- seq_len(nrow(pm_corrected))

  # ----------------------------------------------------------
  # CRITICAL DIFFERENCE FROM 12J:
  # fit ALL promoters before extracting kinases.
  # ----------------------------------------------------------

  fit <- limma::lmFit(
    mat,
    design
  )

  fit2 <- limma::contrasts.fit(
    fit,
    contrast_matrix
  )

  fit2 <- limma::eBayes(
    fit2
  )

  for (contrast_name in colnames(contrast_matrix)) {

    tt <- limma::topTable(
      fit2,
      coef = contrast_name,
      number = Inf,
      sort.by = "none",
      adjust.method = "none"
    )

    if (nrow(tt) != nrow(pm_corrected)) {
      stop("Unexpected limma result row count.")
    }

    res <- data.frame(
      pm_corrected[
        ,
        c(
          "ENTREZID",
          "SYMBOL",
          "ENSEMBL",
          "chr",
          "start",
          "end",
          "strand",
          "width"
        ),
        drop = FALSE
      ],
      mark = mark,
      contrast = contrast_name,
      logFC = tt$logFC,
      AveExpr = tt$AveExpr,
      t = tt$t,
      P.Value = tt$P.Value,
      B = tt$B,
      stringsAsFactors = FALSE
    )

    # --------------------------------------------------------
    # KINASE RESTRICTION OCCURS ONLY AFTER eBAYES
    # --------------------------------------------------------

    kres <- res[
      kinase_idx,
      ,
      drop = FALSE
    ]

    kres$FDR_kinase <- p.adjust(
      kres$P.Value,
      method = "BH"
    )

    kres$Direction <- ifelse(
      kres$logFC >= 0,
      "UP",
      "DOWN"
    )

    all_kinase_results[[paste(mark, contrast_name, sep = "__")]] <- kres

    cat(
      mark,
      contrast_name,
      ":",
      nrow(kres),
      "kinases tested\n"
    )
  }

  model_summary[[mark]] <- data.frame(
    mark = mark,
    promoters_used_for_lmFit = nrow(mat),
    kinase_promoters_for_BH = sum(kinase_idx),
    transformation = "log2(signal + 1)",
    design = "unpaired",
    eBayes_universe = "all_promoters",
    FDR_universe = "KinHub_kinases",
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 7. SAVE ALL KINASE RESULTS
# ============================================================

kinase_results <- do.call(
  rbind,
  all_kinase_results
)

rownames(kinase_results) <- NULL

write.csv(
  kinase_results,
  file.path(
    outdir,
    "12k_allpromoter_eBayes_kinaseBH.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 8. TARGET GENES
# ============================================================

targets <- c(
  "JAK3",
  "EPHB2",
  "MET",
  "MAP3K8",
  "BMPR1A"
)

target_results <- kinase_results[
  kinase_results$SYMBOL %in% targets,
  ,
  drop = FALSE
]

target_results <- target_results[
  order(
    match(
      target_results$SYMBOL,
      targets
    ),
    target_results$mark,
    target_results$contrast
  ),
  ,
  drop = FALSE
]

write.csv(
  target_results,
  file.path(
    outdir,
    "12k_target_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. TOBIAS BENCHMARK
# ============================================================

benchmark <- data.frame(
  SYMBOL = c(
    "JAK3",
    "JAK3",
    "EPHB2",
    "MET"
  ),
  mark = c(
    "H3K27ac",
    "H3K4me3",
    "H3K27ac",
    "H3K27ac"
  ),
  contrast = rep(
    "LPS_vs_RPMI",
    4
  ),
  Tobias_FDR = c(
    0.0065,
    0.0141,
    0.030,
    0.051
  ),
  stringsAsFactors = FALSE
)

benchmark$Our_P <- NA_real_
benchmark$Our_FDR <- NA_real_
benchmark$Our_logFC <- NA_real_

for (i in seq_len(nrow(benchmark))) {

  hit <-
    target_results$SYMBOL == benchmark$SYMBOL[i] &
    target_results$mark == benchmark$mark[i] &
    target_results$contrast == benchmark$contrast[i]

  if (sum(hit) == 1) {

    benchmark$Our_P[i] <-
      target_results$P.Value[hit]

    benchmark$Our_FDR[i] <-
      target_results$FDR_kinase[hit]

    benchmark$Our_logFC[i] <-
      target_results$logFC[hit]
  }
}

benchmark$FDR_difference <-
  benchmark$Our_FDR -
  benchmark$Tobias_FDR

write.csv(
  benchmark,
  file.path(
    outdir,
    "12k_Tobias_benchmark.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. MODEL SUMMARY
# ============================================================

model_summary <- do.call(
  rbind,
  model_summary
)

write.csv(
  model_summary,
  file.path(
    outdir,
    "12k_model_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("TOBIAS BENCHMARK\n")
cat("========================================\n\n")

print(
  benchmark,
  row.names = FALSE
)

cat("\nTarget results — LPS vs RPMI:\n\n")

print(
  target_results[
    target_results$contrast == "LPS_vs_RPMI",
    c(
      "SYMBOL",
      "mark",
      "logFC",
      "P.Value",
      "FDR_kinase"
    )
  ],
  row.names = FALSE
)

cat("\n========================================\n")
cat("12k COMPLETE\n")
cat("========================================\n")

cat(
  "\nFiles written:\n",
  "  12k_allpromoter_eBayes_kinaseBH.csv\n",
  "  12k_target_results.csv\n",
  "  12k_Tobias_benchmark.csv\n",
  "  12k_model_summary.csv\n"
)
