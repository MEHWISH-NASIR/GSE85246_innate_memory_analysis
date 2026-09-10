# ============================================================
# 12f — FORMAL DAY-6 PROMOTER CHROMATIN TEST
#
# GSE85245
#
# Inputs:
#   Step 12e corrected genome-wide promoter matrices
#
# Marks:
#   H3K27ac
#   H3K4me1
#   H3K4me3
#
# Contrasts:
#   LPS vs RPMI
#   BG  vs RPMI
#
# Statistical model:
#   limma
#   condition + replicate block
#   empirical Bayes moderation
#   BH FDR across ALL tested promoters
#
# PRIMARY:
#   normalized-peer-median scaling
#
# SENSITIVITY:
#   RPMI_rep1-matched scaling
#
# IMPORTANT:
#   BigWig promoter means are continuous signal values.
#   Therefore the limma coefficient is a SIGNAL DIFFERENCE,
#   not an RNA-seq log2 fold-change.
# ============================================================

rm(list = ls())

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Package 'limma' is required.")
}

library(limma)

outdir <- "results/12_kinase_reevaluation"

input_file <- file.path(
  outdir,
  "12e_day6_promoter_signal_rescaled.rds"
)

if (!file.exists(input_file)) {
  stop("Missing Step 12e rescaled promoter matrix.")
}

x <- readRDS(input_file)

promoter_manifest <- x$promoter_manifest
sample_manifest <- x$sample_manifest

matrices <- list(
  Primary =
    x$signal_matrix_primary,
  Sensitivity =
    x$signal_matrix_sensitivity
)

marks <- c(
  "H3K27ac",
  "H3K4me1",
  "H3K4me3"
)

targets <- c(
  "JAK3",
  "EPHB2",
  "MET",
  "MAP3K8",
  "BMPR1A"
)

cat("\n========================================\n")
cat("12f — FORMAL PROMOTER LIMMA\n")
cat("========================================\n\n")

cat(
  "Promoters:",
  nrow(promoter_manifest),
  "\n"
)

cat(
  "Samples:",
  nrow(sample_manifest),
  "\n\n"
)

all_results <- list()
target_results <- list()
summary_results <- list()

counter <- 0L


# ============================================================
# ANALYSIS FUNCTION
# ============================================================

run_one_mark <- function(
    signal_matrix,
    method_name,
    mark_name
) {

  idx <- which(
    sample_manifest$mark == mark_name
  )

  if (length(idx) != 6) {
    stop(
      paste(
        mark_name,
        "does not have exactly 6 samples."
      )
    )
  }

  mat <- signal_matrix[
    ,
    idx,
    drop = FALSE
  ]

  condition <- factor(
    sample_manifest$condition[idx],
    levels = c(
      "RPMI",
      "LPS",
      "BG"
    )
  )

  replicate <- factor(
    sample_manifest$replicate[idx],
    levels = c(
      "rep1",
      "rep2"
    )
  )

  cat("\n----------------------------------------\n")
  cat(method_name, "—", mark_name, "\n")
  cat("----------------------------------------\n\n")

  cat("Condition x replicate:\n")

  print(
    table(
      condition,
      replicate
    )
  )

  # ----------------------------------------------------------
  # Paired/block design
  # ----------------------------------------------------------

  design <- model.matrix(
    ~ 0 + condition + replicate
  )

  colnames(design) <- sub(
    "^condition",
    "",
    colnames(design)
  )

  cat("\nDesign matrix:\n")
  print(design)

  cat(
    "\nDesign rank:",
    qr(design)$rank,
    "/",
    ncol(design),
    "\n"
  )

  if (qr(design)$rank != ncol(design)) {
    stop(
      paste(
        "Design matrix is not full rank for",
        mark_name
      )
    )
  }

  # ----------------------------------------------------------
  # Limma
  # ----------------------------------------------------------

  fit <- limma::lmFit(
    mat,
    design
  )

  contrast_matrix <- limma::makeContrasts(
    LPS_vs_RPMI = LPS - RPMI,
    BG_vs_RPMI = BG - RPMI,
    levels = design
  )

  fit2 <- limma::contrasts.fit(
    fit,
    contrast_matrix
  )

  fit2 <- limma::eBayes(
    fit2
  )

  contrast_names <- colnames(
    contrast_matrix
  )

  output <- list()

  for (contrast_name in contrast_names) {

    tt <- limma::topTable(
      fit2,
      coef = contrast_name,
      number = Inf,
      adjust.method = "BH",
      sort.by = "none"
    )

    if (nrow(tt) != nrow(promoter_manifest)) {
      stop(
        paste(
          "Unexpected number of promoters in",
          mark_name,
          contrast_name
        )
      )
    }

    result <- data.frame(
      promoter_manifest,
      Scaling_method =
        method_name,
      Mark =
        mark_name,
      Contrast =
        contrast_name,
      Signal_difference =
        tt$logFC,
      Ave_signal =
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
            "NO_CHANGE"
          )
        ),
      stringsAsFactors = FALSE
    )

    output[[contrast_name]] <- result
  }

  output
}


# ============================================================
# RUN PRIMARY + SENSITIVITY
# ============================================================

for (method_name in names(matrices)) {

  signal_matrix <- matrices[[method_name]]

  if (
    nrow(signal_matrix) !=
      nrow(promoter_manifest)
  ) {
    stop(
      paste(
        "Promoter number mismatch:",
        method_name
      )
    )
  }

  for (mark_name in marks) {

    res <- run_one_mark(
      signal_matrix,
      method_name,
      mark_name
    )

    for (contrast_name in names(res)) {

      counter <- counter + 1L

      z <- res[[contrast_name]]

      key <- paste(
        method_name,
        mark_name,
        contrast_name,
        sep = "__"
      )

      all_results[[key]] <- z

      # -------------------------------------------------------
      # SAVE COMPLETE GENOME-WIDE RESULT
      # -------------------------------------------------------

      filename <- paste0(
        "12f_",
        tolower(method_name),
        "_",
        mark_name,
        "_",
        contrast_name,
        ".csv"
      )

      write.csv(
        z,
        file.path(
          outdir,
          filename
        ),
        row.names = FALSE
      )

      # -------------------------------------------------------
      # SUMMARY
      # -------------------------------------------------------

      summary_results[[key]] <- data.frame(
        Scaling_method =
          method_name,
        Mark =
          mark_name,
        Contrast =
          contrast_name,
        Promoters_tested =
          nrow(z),
        FDR_lt_0.05 =
          sum(
            z$FDR < 0.05,
            na.rm = TRUE
          ),
        FDR_UP =
          sum(
            z$FDR < 0.05 &
              z$Signal_difference > 0,
            na.rm = TRUE
          ),
        FDR_DOWN =
          sum(
            z$FDR < 0.05 &
              z$Signal_difference < 0,
            na.rm = TRUE
          ),
        stringsAsFactors = FALSE
      )

      # -------------------------------------------------------
      # TARGET GENES
      # -------------------------------------------------------

      tg <- z[
        z$SYMBOL %in% targets,
      ]

      tg$Target_order <- match(
        tg$SYMBOL,
        targets
      )

      tg <- tg[
        order(tg$Target_order),
      ]

      tg$Target_order <- NULL

      target_results[[key]] <- tg
    }
  }
}


# ============================================================
# COMBINE SUMMARIES
# ============================================================

summary_table <- do.call(
  rbind,
  summary_results
)

rownames(summary_table) <- NULL

write.csv(
  summary_table,
  file.path(
    outdir,
    "12f_formal_limma_summary.csv"
  ),
  row.names = FALSE
)

target_table <- do.call(
  rbind,
  target_results
)

rownames(target_table) <- NULL

write.csv(
  target_table,
  file.path(
    outdir,
    "12f_target_formal_chromatin_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PRIMARY RESULTS ONLY
# ============================================================

primary_targets <- target_table[
  target_table$Scaling_method == "Primary",
]

primary_targets <- primary_targets[
  order(
    match(
      primary_targets$SYMBOL,
      targets
    ),
    primary_targets$Mark,
    primary_targets$Contrast
  ),
]

write.csv(
  primary_targets,
  file.path(
    outdir,
    "12f_PRIMARY_target_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PRIMARY vs SENSITIVITY ROBUSTNESS
# ============================================================

robustness <- merge(
  target_table[
    target_table$Scaling_method == "Primary",
    c(
      "SYMBOL",
      "Mark",
      "Contrast",
      "Signal_difference",
      "P_value",
      "FDR",
      "Direction"
    )
  ],
  target_table[
    target_table$Scaling_method == "Sensitivity",
    c(
      "SYMBOL",
      "Mark",
      "Contrast",
      "Signal_difference",
      "P_value",
      "FDR",
      "Direction"
    )
  ],
  by = c(
    "SYMBOL",
    "Mark",
    "Contrast"
  ),
  suffixes = c(
    "_Primary",
    "_Sensitivity"
  )
)

robustness$Same_direction <-
  robustness$Direction_Primary ==
  robustness$Direction_Sensitivity

robustness$Significant_primary <-
  robustness$FDR_Primary < 0.05

robustness$Significant_sensitivity <-
  robustness$FDR_Sensitivity < 0.05

robustness$Robust_FDR05 <-
  robustness$Significant_primary &
  robustness$Significant_sensitivity &
  robustness$Same_direction

write.csv(
  robustness,
  file.path(
    outdir,
    "12f_target_scaling_robustness.csv"
  ),
  row.names = FALSE
)


# ============================================================
# CONSOLE OUTPUT
# ============================================================

cat("\n========================================\n")
cat("GENOME-WIDE DIFFERENTIAL SUMMARY\n")
cat("========================================\n\n")

print(
  summary_table,
  row.names = FALSE
)

cat("\n========================================\n")
cat("PRIMARY TARGET RESULTS\n")
cat("========================================\n\n")

print(
  primary_targets[
    ,
    c(
      "SYMBOL",
      "Mark",
      "Contrast",
      "Signal_difference",
      "P_value",
      "FDR",
      "Direction"
    )
  ],
  row.names = FALSE
)

cat("\n========================================\n")
cat("SCALING ROBUSTNESS\n")
cat("========================================\n\n")

print(
  robustness[
    ,
    c(
      "SYMBOL",
      "Mark",
      "Contrast",
      "FDR_Primary",
      "FDR_Sensitivity",
      "Direction_Primary",
      "Direction_Sensitivity",
      "Robust_FDR05"
    )
  ],
  row.names = FALSE
)

# ------------------------------------------------------------
# Tobias benchmark subset
# ------------------------------------------------------------

benchmark <- primary_targets[
  (
    primary_targets$SYMBOL == "JAK3" &
    primary_targets$Mark %in%
      c("H3K27ac", "H3K4me3") &
    primary_targets$Contrast ==
      "LPS_vs_RPMI"
  ) |
  (
    primary_targets$SYMBOL == "MET" &
    primary_targets$Mark ==
      "H3K27ac" &
    primary_targets$Contrast ==
      "LPS_vs_RPMI"
  ) |
  (
    primary_targets$SYMBOL == "EPHB2" &
    primary_targets$Mark ==
      "H3K27ac" &
    primary_targets$Contrast ==
      "LPS_vs_RPMI"
  ),
]

cat("\n========================================\n")
cat("TOBIAS BENCHMARK CHECK\n")
cat("========================================\n\n")

print(
  benchmark[
    ,
    c(
      "SYMBOL",
      "Mark",
      "Contrast",
      "Signal_difference",
      "P_value",
      "FDR",
      "Direction"
    )
  ],
  row.names = FALSE
)

cat("\nReference values reported by Tobias:\n")
cat(" JAK3 H3K27ac FDR ~ 0.0065\n")
cat(" JAK3 H3K4me3 FDR ~ 0.0141\n")
cat(" MET H3K27ac FDR ~ 0.051\n")
cat(" EPHB2 H3K27ac FDR ~ 0.030\n")
cat("\nThese are comparison benchmarks, not forced targets.\n")

cat("\nSaved:\n")
cat(
  " ",
  file.path(
    outdir,
    "12f_formal_limma_summary.csv"
  ),
  "\n"
)
cat(
  " ",
  file.path(
    outdir,
    "12f_PRIMARY_target_results.csv"
  ),
  "\n"
)
cat(
  " ",
  file.path(
    outdir,
    "12f_target_scaling_robustness.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("12f COMPLETE\n")
cat("========================================\n")
