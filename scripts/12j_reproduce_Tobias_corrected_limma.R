
# ============================================================
# 12j — Reproduce Tobias corrected Day-6 kinase promoter model
#
# Fixed model:
#   - corrected primary scaling
#   - KinHub kinase universe
#   - log2(signal + 1)
#   - unpaired limma
#   - BH within kinase family
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(limma)
})

outdir <- "results/12_kinase_reevaluation"

rds_file <- file.path(
  outdir,
  "12e_day6_promoter_signal_rescaled.rds"
)

kinase_file <-
  "results/03_memory_kinases/03_Day6_KinHub_all_kinases.csv"

if (!file.exists(rds_file)) {
  stop("Missing 12e corrected RDS.")
}

if (!file.exists(kinase_file)) {
  stop("Missing KinHub kinase table.")
}

# ------------------------------------------------------------
# 1. Load corrected promoter signals
# ------------------------------------------------------------

x <- readRDS(rds_file)

pm <- x$promoter_manifest
sm <- x$sample_manifest

signal <- x$signal_matrix_primary

if (is.null(signal)) {
  stop("signal_matrix_primary missing.")
}

if (nrow(signal) != nrow(pm)) {
  stop("Promoter manifest and signal matrix do not align.")
}

cat("\n========================================\n")
cat("12j TOBIAS REPRODUCTION\n")
cat("========================================\n")

cat("\nAll promoter rows:", nrow(signal), "\n")
cat("Samples:", ncol(signal), "\n")

# ------------------------------------------------------------
# 2. Read predefined kinase universe
# ------------------------------------------------------------

kin <- read.csv(
  kinase_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"SYMBOL" %in% names(kin)) {
  stop("Kinase table has no SYMBOL column.")
}

kinase_symbols <- unique(
  trimws(
    as.character(
      kin$SYMBOL[
        !is.na(kin$SYMBOL) &
        kin$SYMBOL != ""
      ]
    )
  )
)

keep <- which(
  pm$SYMBOL %in% kinase_symbols
)

kinase_pm <- pm[
  keep,
  ,
  drop = FALSE
]

kinase_signal <- signal[
  keep,
  ,
  drop = FALSE
]

cat(
  "Kinases requested:",
  length(kinase_symbols),
  "\n"
)

cat(
  "Kinase promoters mapped:",
  nrow(kinase_signal),
  "\n"
)

if (nrow(kinase_signal) < 400) {
  stop("Unexpectedly small kinase universe.")
}

# ------------------------------------------------------------
# 3. Confirm target genes
# ------------------------------------------------------------

targets <- c(
  "JAK3",
  "EPHB2",
  "MET",
  "MAP3K8",
  "BMPR1A"
)

cat("\nTarget check:\n")

print(
  kinase_pm[
    kinase_pm$SYMBOL %in% targets,
    c(
      "ENTREZID",
      "SYMBOL",
      "ENSEMBL",
      "chr",
      "start",
      "end"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

missing_targets <- setdiff(
  targets,
  kinase_pm$SYMBOL
)

if (length(missing_targets) > 0) {
  stop(
    paste(
      "Missing targets:",
      paste(missing_targets, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 4. Tobias: log transform
#
# Exact transform was not specified by Tobias.
# Primary explicit implementation:
# log2(signal + 1)
# ------------------------------------------------------------

mat <- log2(
  kinase_signal + 1
)

if (any(!is.finite(mat))) {
  stop("Non-finite values after log transformation.")
}

# ------------------------------------------------------------
# 5. Analyze each histone mark independently
# ------------------------------------------------------------

marks <- c(
  "H3K27ac",
  "H3K4me1",
  "H3K4me3"
)

all_results <- list()
counter <- 0L

for (mark_name in marks) {

  cat(
    "\n----------------------------------------\n",
    mark_name,
    "\n----------------------------------------\n",
    sep = ""
  )

  idx <- which(
    sm$mark == mark_name
  )

  if (length(idx) != 6) {
    stop(
      paste(
        "Expected 6 samples for",
        mark_name,
        "but found",
        length(idx)
      )
    )
  }

  y <- mat[
    ,
    idx,
    drop = FALSE
  ]

  meta <- sm[
    idx,
    ,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # UNPAIRED DESIGN
  # ----------------------------------------------------------

  group <- factor(
    meta$condition,
    levels = c(
      "RPMI",
      "LPS",
      "BG"
    )
  )

  design <- model.matrix(
    ~0 + group
  )

  colnames(design) <- levels(group)

  cat("\nSamples:\n")

  print(
    data.frame(
      sample_id = meta$sample_id,
      condition = meta$condition,
      replicate = meta$replicate
    ),
    row.names = FALSE
  )

  cat("\nDesign:\n")
  print(design)

  if (qr(design)$rank != ncol(design)) {
    stop(
      paste(
        "Design not full rank:",
        mark_name
      )
    )
  }

  # ----------------------------------------------------------
  # 6. limma
  # ----------------------------------------------------------

  fit <- limma::lmFit(
    y,
    design
  )

  cont <- limma::makeContrasts(
    LPS_vs_RPMI = LPS - RPMI,
    BG_vs_RPMI  = BG  - RPMI,
    levels = design
  )

  fit2 <- limma::contrasts.fit(
    fit,
    cont
  )

  fit2 <- limma::eBayes(
    fit2
  )

  # ----------------------------------------------------------
  # 7. BH within kinase testing family
  # ----------------------------------------------------------

  for (contrast_name in colnames(cont)) {

    tt <- limma::topTable(
      fit2,
      coef = contrast_name,
      number = Inf,
      adjust.method = "BH",
      sort.by = "none"
    )

    if (nrow(tt) != nrow(kinase_pm)) {
      stop(
        paste(
          "Unexpected result count:",
          mark_name,
          contrast_name
        )
      )
    }

    counter <- counter + 1L

    all_results[[counter]] <- data.frame(
      kinase_pm,

      Mark = mark_name,
      Contrast = contrast_name,

      Log_signal_difference =
        tt$logFC,

      AveExpr =
        tt$AveExpr,

      t =
        tt$t,

      P_value =
        tt$P.Value,

      FDR =
        tt$adj.P.Val,

      B =
        tt$B,

      Direction =
        ifelse(
          tt$logFC > 0,
          "UP",
          ifelse(
            tt$logFC < 0,
            "DOWN",
            "FLAT"
          )
        ),

      stringsAsFactors = FALSE
    )
  }
}

# ------------------------------------------------------------
# 8. Combine full kinase results
# ------------------------------------------------------------

res <- do.call(
  rbind,
  all_results
)

write.csv(
  res,
  file.path(
    outdir,
    "12j_all_kinase_corrected_limma.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Target summary
# ------------------------------------------------------------

target_res <- res[
  res$SYMBOL %in% targets,
  ,
  drop = FALSE
]

target_res <- target_res[
  order(
    match(
      target_res$SYMBOL,
      targets
    ),
    match(
      target_res$Mark,
      marks
    ),
    target_res$Contrast
  ),
  ,
  drop = FALSE
]

write.csv(
  target_res,
  file.path(
    outdir,
    "12j_target_results.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Tobias benchmark rows
# ------------------------------------------------------------

benchmark <- target_res[
  (
    target_res$SYMBOL == "JAK3" &
    target_res$Mark %in%
      c("H3K27ac", "H3K4me3")
  ) |
  (
    target_res$SYMBOL == "EPHB2" &
    target_res$Mark == "H3K27ac"
  ) |
  (
    target_res$SYMBOL == "MET" &
    target_res$Mark == "H3K27ac"
  ),
  ,
  drop = FALSE
]

benchmark$Tobias_FDR <- NA_real_

benchmark$Tobias_FDR[
  benchmark$SYMBOL == "JAK3" &
  benchmark$Mark == "H3K27ac" &
  benchmark$Contrast == "LPS_vs_RPMI"
] <- 0.0065

benchmark$Tobias_FDR[
  benchmark$SYMBOL == "JAK3" &
  benchmark$Mark == "H3K4me3" &
  benchmark$Contrast == "LPS_vs_RPMI"
] <- 0.0141

benchmark$Tobias_FDR[
  benchmark$SYMBOL == "EPHB2" &
  benchmark$Mark == "H3K27ac"
] <- 0.030

benchmark$Tobias_FDR[
  benchmark$SYMBOL == "MET" &
  benchmark$Mark == "H3K27ac"
] <- 0.051

write.csv(
  benchmark,
  file.path(
    outdir,
    "12j_Tobias_benchmark_rows.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Console output
# ------------------------------------------------------------

cat("\n========================================\n")
cat("TARGET RESULTS\n")
cat("========================================\n\n")

print(
  target_res[
    ,
    c(
      "SYMBOL",
      "Mark",
      "Contrast",
      "Log_signal_difference",
      "P_value",
      "FDR",
      "Direction"
    )
  ],
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("TOBIAS BENCHMARK ROWS\n")
cat("========================================\n\n")

print(
  benchmark[
    ,
    c(
      "SYMBOL",
      "Mark",
      "Contrast",
      "Log_signal_difference",
      "P_value",
      "FDR",
      "Tobias_FDR",
      "Direction"
    )
  ],
  row.names = FALSE,
  digits = 6
)

cat("\nFiles saved:\n")
cat(
  "  12j_all_kinase_corrected_limma.csv\n"
)
cat(
  "  12j_target_results.csv\n"
)
cat(
  "  12j_Tobias_benchmark_rows.csv\n"
)

cat("\n========================================\n")
cat("12j COMPLETE\n")
cat("========================================\n")
