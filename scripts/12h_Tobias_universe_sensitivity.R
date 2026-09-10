# ============================================================
# 12h — TOBIAS MODEL: UNIVERSE / SCALING / LOG SENSITIVITY
#
# Tobias confirmed:
#   - log-transformed promoter signals
#   - kinase genes only
#   - unpaired limma design
#
# We test:
#   2 scaling methods
#   4 log pseudocounts
#   3 kinase universes
#   2 contrasts (LPS-vs-RPMI, BG-vs-RPMI)
#
# IMPORTANT:
#   KinHub497 = primary/literal kinase universe.
#   Day6_24 and HIGH_10 are sensitivity analyses only.
# ============================================================

rm(list = ls())

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("limma required.")
}

library(limma)

outdir <- "results/12_kinase_reevaluation"

rds_file <- file.path(
  outdir,
  "12e_day6_promoter_signal_rescaled.rds"
)

file_497 <-
  "results/03_memory_kinases/03_Day6_KinHub_all_kinases.csv"

file_24 <-
  "results/03_memory_kinases/03_memory_kinases_FDR05.csv"

file_10 <-
  "results/05_integrated_trajectory/05_HIGH_priority_for_chromatin.csv"

stopifnot(
  file.exists(rds_file),
  file.exists(file_497),
  file.exists(file_24),
  file.exists(file_10)
)

x <- readRDS(rds_file)

pm <- x$promoter_manifest
sm <- x$sample_manifest

matrices <- list(
  PeerMedian =
    x$signal_matrix_primary,
  RPMIrep1Matched =
    x$signal_matrix_sensitivity
)

read_symbols <- function(f) {

  z <- read.csv(
    f,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  unique(
    z$SYMBOL[
      !is.na(z$SYMBOL) &
      z$SYMBOL != ""
    ]
  )
}

universes <- list(
  KinHub497 =
    read_symbols(file_497),
  Day6_24 =
    read_symbols(file_24),
  HIGH_10 =
    read_symbols(file_10)
)

targets <- c(
  "JAK3",
  "MET",
  "MAP3K8",
  "BMPR1A",
  "EPHB2"
)

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

cat("\n========================================\n")
cat("12h — TOBIAS UNIVERSE SENSITIVITY\n")
cat("========================================\n\n")

cat("Available promoter rows:", nrow(pm), "\n\n")

for (u in names(universes)) {

  mapped <- intersect(
    universes[[u]],
    pm$SYMBOL
  )

  cat(
    u,
    ": requested",
    length(universes[[u]]),
    "| mapped",
    length(mapped),
    "\n"
  )

  cat(
    "  targets:",
    paste(
      intersect(targets, mapped),
      collapse = ", "
    ),
    "\n"
  )
}


# ============================================================
# RUN ALL COMBINATIONS
# ============================================================

all_results <- list()
target_results <- list()
summary_results <- list()

counter <- 0L

for (scaling_name in names(matrices)) {

  full_signal <- matrices[[scaling_name]]

  for (universe_name in names(universes)) {

    symbols <- universes[[universe_name]]

    keep <- which(
      pm$SYMBOL %in% symbols
    )

    manifest <- pm[
      keep,
      ,
      drop = FALSE
    ]

    signal0 <- full_signal[
      keep,
      ,
      drop = FALSE
    ]

    if (nrow(signal0) < 5) {
      stop(
        paste(
          "Too few genes:",
          universe_name
        )
      )
    }

    for (pc in pseudocounts) {

      transform_name <- paste0(
        "log2_plus_",
        format(
          pc,
          scientific = FALSE,
          trim = TRUE
        )
      )

      signal <- log2(
        signal0 + pc
      )

      if (any(!is.finite(signal))) {
        stop("Non-finite transformed values.")
      }

      for (mark_name in marks) {

        idx <- which(
          sm$mark == mark_name
        )

        if (length(idx) != 6) {
          stop(
            paste(
              "Expected 6 samples for",
              mark_name
            )
          )
        }

        y <- signal[
          ,
          idx,
          drop = FALSE
        ]

        condition <- factor(
          sm$condition[idx],
          levels = c(
            "RPMI",
            "LPS",
            "BG"
          )
        )

        # Tobias: unpaired design
        design <- model.matrix(
          ~ 0 + condition
        )

        colnames(design) <- sub(
          "^condition",
          "",
          colnames(design)
        )

        fit <- limma::lmFit(
          y,
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

          z <- data.frame(
            manifest,
            Scaling =
              scaling_name,
            Universe =
              universe_name,
            Universe_N =
              nrow(manifest),
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
            P_value =
              tt$P.Value,
            FDR =
              tt$adj.P.Val,
            t =
              tt$t,
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

          counter <- counter + 1L

          key <- paste(
            scaling_name,
            universe_name,
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
              Scaling =
                scaling_name,
              Universe =
                universe_name,
              Universe_N =
                nrow(manifest),
              Transform =
                transform_name,
              Mark =
                mark_name,
              Contrast =
                contrast_name,
              FDR05 =
                sum(
                  z$FDR < 0.05,
                  na.rm = TRUE
                ),
              stringsAsFactors = FALSE
            )
        }
      }
    }
  }
}


# ============================================================
# COMBINE
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
    "12h_all_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  target_table,
  file.path(
    outdir,
    "12h_target_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  summary_table,
  file.path(
    outdir,
    "12h_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# JAK3 BENCHMARK
#
# Tobias:
# H3K27ac ~0.0065
# H3K4me3 ~0.0141
#
# We deliberately allow BOTH stimulus contrasts because
# Tobias's email did not explicitly specify which contrast
# these two quoted FDRs correspond to.
# ============================================================

jak <- target_table[
  target_table$SYMBOL == "JAK3" &
  target_table$Mark %in%
    c("H3K27ac", "H3K4me3"),
]

jak$Tobias_FDR <- ifelse(
  jak$Mark == "H3K27ac",
  0.0065,
  0.0141
)

jak$Log10_distance <- abs(
  log10(
    pmax(
      jak$FDR,
      1e-300
    )
  ) -
  log10(
    jak$Tobias_FDR
  )
)

write.csv(
  jak,
  file.path(
    outdir,
    "12h_JAK3_benchmark_all_models.csv"
  ),
  row.names = FALSE
)

best_jak <- jak[
  order(
    jak$Mark,
    jak$Log10_distance
  ),
]

best_jak <- do.call(
  rbind,
  lapply(
    split(
      best_jak,
      best_jak$Mark
    ),
    function(z) head(z, 12)
  )
)

write.csv(
  best_jak,
  file.path(
    outdir,
    "12h_JAK3_closest_models.csv"
  ),
  row.names = FALSE
)


# ============================================================
# MET benchmark
# ============================================================

met <- target_table[
  target_table$SYMBOL == "MET" &
  target_table$Mark == "H3K27ac",
]

if (nrow(met) > 0) {

  met$Tobias_FDR <- 0.051

  met$Log10_distance <- abs(
    log10(
      pmax(
        met$FDR,
        1e-300
      )
    ) -
    log10(0.051)
  )

  met <- met[
    order(
      met$Log10_distance
    ),
  ]

  write.csv(
    met,
    file.path(
      outdir,
      "12h_MET_benchmark.csv"
    ),
    row.names = FALSE
  )
}


# ============================================================
# CONSOLE OUTPUT
# ============================================================

cat("\n========================================\n")
cat("UNIVERSE MAPPING\n")
cat("========================================\n\n")

for (u in names(universes)) {

  mapped <- intersect(
    universes[[u]],
    pm$SYMBOL
  )

  cat(
    u,
    "=",
    length(mapped),
    "\n"
  )

  cat(
    " targets:",
    paste(
      intersect(
        targets,
        mapped
      ),
      collapse = ", "
    ),
    "\n\n"
  )
}

cat("\n========================================\n")
cat("CLOSEST JAK3 MODELS\n")
cat("========================================\n\n")

print(
  best_jak[
    ,
    c(
      "SYMBOL",
      "Mark",
      "Scaling",
      "Universe",
      "Universe_N",
      "Transform",
      "Contrast",
      "P_value",
      "FDR",
      "Tobias_FDR",
      "Direction",
      "Log10_distance"
    )
  ],
  row.names = FALSE
)

if (exists("met") && nrow(met) > 0) {

  cat("\n========================================\n")
  cat("CLOSEST MET H3K27ac MODELS\n")
  cat("========================================\n\n")

  print(
    head(
      met[
        ,
        c(
          "Scaling",
          "Universe",
          "Universe_N",
          "Transform",
          "Contrast",
          "P_value",
          "FDR",
          "Tobias_FDR",
          "Direction",
          "Log10_distance"
        )
      ],
      20
    ),
    row.names = FALSE
  )
}


cat("\n========================================\n")
cat("MAP3K8 / BMPR1A CHECK\n")
cat("========================================\n\n")

neg <- target_table[
  target_table$SYMBOL %in%
    c(
      "MAP3K8",
      "BMPR1A"
    ),
]

neg_best <- neg[
  order(
    neg$FDR
  ),
]

print(
  head(
    neg_best[
      ,
      c(
        "SYMBOL",
        "Mark",
        "Scaling",
        "Universe",
        "Transform",
        "Contrast",
        "P_value",
        "FDR",
        "Direction"
      )
    ],
    30
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("EPHB2 DIRECTION CHECK — KINHUB497 ONLY\n")
cat("========================================\n\n")

eph <- target_table[
  target_table$SYMBOL == "EPHB2" &
  target_table$Universe == "KinHub497" &
  target_table$Mark == "H3K27ac",
]

print(
  eph[
    ,
    c(
      "Scaling",
      "Transform",
      "Contrast",
      "Log_signal_difference",
      "P_value",
      "FDR",
      "Direction"
    )
  ],
  row.names = FALSE
)

cat("\n========================================\n")
cat("12h COMPLETE\n")
cat("========================================\n")
