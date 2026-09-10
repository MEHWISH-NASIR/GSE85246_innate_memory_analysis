# ============================================================
# 12g — TOBIAS MODEL REPRODUCTION / SENSITIVITY
#
# Confirmed by Tobias:
#   1. promoter signals were log transformed
#   2. statistical testing was limited to kinase genes
#   3. design was unpaired
#
# Purpose:
#   Reproduce the formal Day-6 promoter chromatin analysis
#   while testing plausible pseudocounts for log transform.
#
# IMPORTANT:
#   All transformation results are retained.
#   Tobias's reported FDRs are used only as reproduction
#   benchmarks, NOT to determine significance.
# ============================================================

rm(list = ls())

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("limma is required.")
}

library(limma)

outdir <- "results/12_kinase_reevaluation"

matrix_file <- file.path(
  outdir,
  "12e_day6_promoter_signal_rescaled.rds"
)

kinase_file <-
  "results/03_memory_kinases/03_Day6_KinHub_all_kinases.csv"

if (!file.exists(matrix_file)) {
  stop("Missing Step 12e corrected promoter matrix.")
}

if (!file.exists(kinase_file)) {
  stop("Missing KinHub kinase table.")
}

x <- readRDS(matrix_file)

promoter_manifest <- x$promoter_manifest
sample_manifest <- x$sample_manifest

# Primary correction from Step 12e:
# NotNormalized RPMI_rep2 matched to normalized peers.
signal <- x$signal_matrix_primary

kin <- read.csv(
  kinase_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("12g — TOBIAS MODEL REPRODUCTION\n")
cat("========================================\n\n")

cat("Total promoters:", nrow(signal), "\n")
cat("KinHub rows:", nrow(kin), "\n")

# ============================================================
# 1. DEFINE KINASE UNIVERSE
# ============================================================

kinase_symbols <- unique(
  kin$SYMBOL[
    !is.na(kin$SYMBOL) &
    kin$SYMBOL != ""
  ]
)

keep <- which(
  promoter_manifest$SYMBOL %in%
    kinase_symbols
)

signal_kinase <- signal[
  keep,
  ,
  drop = FALSE
]

kinase_manifest <- promoter_manifest[
  keep,
  ,
  drop = FALSE
]

cat(
  "Kinase promoters successfully mapped:",
  nrow(signal_kinase),
  "\n"
)

cat(
  "KinHub symbols not mapped to hg19 promoter set:",
  length(
    setdiff(
      kinase_symbols,
      kinase_manifest$SYMBOL
    )
  ),
  "\n"
)

if (nrow(signal_kinase) < 400) {
  stop(
    "Unexpectedly small kinase universe; stop before testing."
  )
}

# ============================================================
# 2. TARGET CHECK
# ============================================================

targets <- c(
  "JAK3",
  "MET",
  "MAP3K8",
  "BMPR1A",
  "EPHB2"
)

cat("\nTarget kinase check:\n")

print(
  kinase_manifest[
    kinase_manifest$SYMBOL %in% targets,
    c(
      "ENTREZID",
      "SYMBOL",
      "ENSEMBL",
      "chr",
      "start",
      "end"
    )
  ],
  row.names = FALSE
)

missing_targets <- setdiff(
  targets,
  kinase_manifest$SYMBOL
)

if (length(missing_targets) > 0) {
  warning(
    paste(
      "Targets absent from kinase universe:",
      paste(missing_targets, collapse = ", ")
    )
  )
}

# ============================================================
# 3. TRANSFORMS TO TEST
# ============================================================

pseudocounts <- c(
  1,
  0.1,
  0.01,
  0.001
)

marks <- c(
  "H3K27ac",
  "H3K4me1",
  "H3K4me3"
)

all_results <- list()
target_results <- list()
summary_results <- list()

# ============================================================
# 4. RUN UNPAIRED LIMMA
# ============================================================

for (pc in pseudocounts) {

  transform_name <- paste0(
    "log2_plus_",
    format(
      pc,
      scientific = FALSE,
      trim = TRUE
    )
  )

  cat("\n\n========================================\n")
  cat("TRANSFORM:", transform_name, "\n")
  cat("========================================\n")

  transformed <- log2(
    signal_kinase + pc
  )

  if (any(!is.finite(transformed))) {
    stop(
      paste(
        "Non-finite transformed values:",
        transform_name
      )
    )
  }

  for (mark_name in marks) {

    idx <- which(
      sample_manifest$mark ==
        mark_name
    )

    if (length(idx) != 6) {
      stop(
        paste(
          "Expected 6 samples for",
          mark_name
        )
      )
    }

    mat <- transformed[
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

    # --------------------------------------------------------
    # UNPAIRED MODEL:
    # no replicate / donor blocking term
    # --------------------------------------------------------

    design <- model.matrix(
      ~ 0 + condition
    )

    colnames(design) <- sub(
      "^condition",
      "",
      colnames(design)
    )

    if (qr(design)$rank != ncol(design)) {
      stop("Design is not full rank.")
    }

    fit <- limma::lmFit(
      mat,
      design
    )

    cont <- limma::makeContrasts(
      LPS_vs_RPMI =
        LPS - RPMI,
      BG_vs_RPMI =
        BG - RPMI,
      levels = design
    )

    fit2 <- limma::contrasts.fit(
      fit,
      cont
    )

    # Default empirical Bayes:
    # no trend / robust options unless Tobias specifies them.
    fit2 <- limma::eBayes(
      fit2
    )

    for (
      contrast_name in
      colnames(cont)
    ) {

      tt <- limma::topTable(
        fit2,
        coef = contrast_name,
        number = Inf,
        adjust.method = "BH",
        sort.by = "none"
      )

      if (
        nrow(tt) !=
          nrow(kinase_manifest)
      ) {
        stop(
          "Unexpected topTable row number."
        )
      }

      z <- data.frame(
        kinase_manifest,
        Transform =
          transform_name,
        Pseudocount =
          pc,
        Mark =
          mark_name,
        Contrast =
          contrast_name,
        Log_signal_difference =
          tt$logFC,
        Ave_log_signal =
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

      key <- paste(
        transform_name,
        mark_name,
        contrast_name,
        sep = "__"
      )

      all_results[[key]] <- z

      target_results[[key]] <- z[
        z$SYMBOL %in% targets,
      ]

      summary_results[[key]] <-
        data.frame(
          Transform =
            transform_name,
          Pseudocount =
            pc,
          Mark =
            mark_name,
          Contrast =
            contrast_name,
          Kinases_tested =
            nrow(z),
          FDR_lt_0.05 =
            sum(
              z$FDR < 0.05,
              na.rm = TRUE
            ),
          stringsAsFactors = FALSE
        )
    }
  }
}

# ============================================================
# 5. COMBINE AND SAVE
# ============================================================

all_table <- do.call(
  rbind,
  all_results
)

target_table <- do.call(
  rbind,
  target_results
)

summary_table <- do.call(
  rbind,
  summary_results
)

rownames(all_table) <- NULL
rownames(target_table) <- NULL
rownames(summary_table) <- NULL

write.csv(
  all_table,
  file.path(
    outdir,
    "12g_all_kinase_formal_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  target_table,
  file.path(
    outdir,
    "12g_target_transform_sensitivity.csv"
  ),
  row.names = FALSE
)

write.csv(
  summary_table,
  file.path(
    outdir,
    "12g_transform_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 6. TOBIAS BENCHMARKS
# ============================================================

# Reported after normalization:
#
# JAK3:
#   H3K27ac LPS vs RPMI FDR ~ 0.0065
#   H3K4me3 LPS vs RPMI FDR ~ 0.0141
#
# MET:
#   H3K27ac LPS vs RPMI FDR ~ 0.051
#
# EPHB2:
#   H3K27ac same-direction FDR ~ 0.030
#
# MAP3K8:
#   not significant; reported range ~0.22-0.88
#
# BMPR1A:
#   not significant; reported range ~0.62-0.96

benchmark <- data.frame(
  SYMBOL = c(
    "JAK3",
    "JAK3",
    "MET",
    "EPHB2"
  ),
  Mark = c(
    "H3K27ac",
    "H3K4me3",
    "H3K27ac",
    "H3K27ac"
  ),
  Contrast = rep(
    "LPS_vs_RPMI",
    4
  ),
  Tobias_FDR = c(
    0.0065,
    0.0141,
    0.051,
    0.030
  ),
  stringsAsFactors = FALSE
)

comparison <- merge(
  target_table,
  benchmark,
  by = c(
    "SYMBOL",
    "Mark",
    "Contrast"
  )
)

comparison$FDR_ratio_to_Tobias <-
  comparison$FDR /
  comparison$Tobias_FDR

comparison$Abs_log10_FDR_distance <-
  abs(
    log10(
      pmax(
        comparison$FDR,
        1e-300
      )
    ) -
    log10(
      comparison$Tobias_FDR
    )
  )

write.csv(
  comparison,
  file.path(
    outdir,
    "12g_Tobias_benchmark_comparison.csv"
  ),
  row.names = FALSE
)

# Overall benchmark distance per transform
benchmark_score <- aggregate(
  Abs_log10_FDR_distance ~
    Transform + Pseudocount,
  data = comparison,
  FUN = mean
)

benchmark_score <- benchmark_score[
  order(
    benchmark_score$
      Abs_log10_FDR_distance
  ),
]

write.csv(
  benchmark_score,
  file.path(
    outdir,
    "12g_Tobias_benchmark_scores.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 7. CONSOLE REPORT
# ============================================================

cat("\n\n========================================\n")
cat("TRANSFORM SUMMARY\n")
cat("========================================\n\n")

print(
  summary_table,
  row.names = FALSE
)

cat("\n========================================\n")
cat("KEY LPS-vs-RPMI TARGET RESULTS\n")
cat("========================================\n\n")

key_targets <- target_table[
  target_table$Contrast ==
    "LPS_vs_RPMI" &
  target_table$SYMBOL %in%
    targets,
]

print(
  key_targets[
    ,
    c(
      "Transform",
      "SYMBOL",
      "Mark",
      "Log_signal_difference",
      "P_value",
      "FDR",
      "Direction"
    )
  ],
  row.names = FALSE
)

cat("\n========================================\n")
cat("TOBIAS BENCHMARK COMPARISON\n")
cat("========================================\n\n")

print(
  comparison[
    ,
    c(
      "Transform",
      "SYMBOL",
      "Mark",
      "FDR",
      "Tobias_FDR",
      "FDR_ratio_to_Tobias",
      "Abs_log10_FDR_distance"
    )
  ],
  row.names = FALSE
)

cat("\n========================================\n")
cat("OVERALL BENCHMARK DISTANCE\n")
cat("========================================\n\n")

print(
  benchmark_score,
  row.names = FALSE
)

cat("\nLower distance = closer numerical reproduction.\n")
cat(
  "Do NOT choose a transform from JAK3 alone;",
  "consider the complete pattern.\n"
)

cat("\n========================================\n")
cat("MAP3K8 / BMPR1A NON-SIGNIFICANCE CHECK\n")
cat("========================================\n\n")

check_negative <- target_table[
  target_table$Contrast ==
    "LPS_vs_RPMI" &
  target_table$SYMBOL %in%
    c(
      "MAP3K8",
      "BMPR1A"
    ),
]

print(
  check_negative[
    ,
    c(
      "Transform",
      "SYMBOL",
      "Mark",
      "P_value",
      "FDR",
      "Direction"
    )
  ],
  row.names = FALSE
)

cat("\nSaved:\n")
cat(
  "  results/12_kinase_reevaluation/",
  "12g_target_transform_sensitivity.csv\n",
  sep = ""
)
cat(
  "  results/12_kinase_reevaluation/",
  "12g_Tobias_benchmark_comparison.csv\n",
  sep = ""
)
cat(
  "  results/12_kinase_reevaluation/",
  "12g_Tobias_benchmark_scores.csv\n",
  sep = ""
)

cat("\n========================================\n")
cat("12g COMPLETE\n")
cat("========================================\n")
