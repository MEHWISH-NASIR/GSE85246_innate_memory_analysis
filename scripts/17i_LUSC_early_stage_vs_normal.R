# ============================================================
# 17i — TCGA-LUSC early-stage tumor vs normal analysis
#
# Primary early-stage definition:
#   Stage I tumors only
#
# Secondary sensitivity:
#   Stage I + Stage II tumors
#
# Genes:
#   JAK3
#   EPHB2
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

expr_file <-
  "results/17_TCGA_LUSC/17a_LUSC_JAK3_EPHB2_expression.csv"

clin_file <-
  "results/17_TCGA_LUSC/17c_LUSC_clinical_metadata.csv"

prev_file <-
  "results/17_TCGA_LUSC/17b_LUSC_prevalence_summary.csv"

outdir <- "results/17_TCGA_LUSC"

expr <- read.csv(
  expr_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

clin <- read.csv(
  clin_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

prev <- read.csv(
  prev_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# Patient-level expression
# ------------------------------------------------------------

expr <- expr[
  expr$extraction_status == "OK" &
  !is.na(expr$TPM),
]

patient <- aggregate(
  TPM ~ gene_name + patient_id + sample_type,
  data = expr,
  FUN = median
)

patient$log2TPM1 <- log2(patient$TPM + 1)

# ------------------------------------------------------------
# Normalize stage
# ------------------------------------------------------------

normalize_stage <- function(x) {

  x <- toupper(trimws(as.character(x)))

  out <- rep(
    NA_character_,
    length(x)
  )

  out[grepl("^STAGE IV", x)] <- "Stage IV"
  out[grepl("^STAGE III", x)] <- "Stage III"

  out[
    grepl("^STAGE II", x) &
    is.na(out)
  ] <- "Stage II"

  out[
    grepl("^STAGE I", x) &
    is.na(out)
  ] <- "Stage I"

  out
}

clin$Stage_group <-
  normalize_stage(
    clin$pathologic_stage
  )

# ------------------------------------------------------------
# Attach stage to tumors
# ------------------------------------------------------------

tumor <- patient[
  patient$sample_type == "Primary Tumor",
]

tumor <- merge(
  tumor,
  clin[
    ,
    c(
      "patient_id",
      "Stage_group"
    )
  ],
  by = "patient_id",
  all.x = TRUE
)

normal <- patient[
  patient$sample_type == "Solid Tissue Normal",
]

# ------------------------------------------------------------
# Analysis
# ------------------------------------------------------------

genes <- c(
  "JAK3",
  "EPHB2"
)

primary_results <- list()
sensitivity_results <- list()

for (gene in genes) {

  tgene <- tumor[
    tumor$gene_name == gene,
  ]

  ngene <- normal[
    normal$gene_name == gene,
  ]

  # ==========================================================
  # PRIMARY: Stage I vs normal
  # ==========================================================

  early1 <- tgene[
    !is.na(tgene$Stage_group) &
    tgene$Stage_group == "Stage I",
  ]

  wt1 <- wilcox.test(
    early1$log2TPM1,
    ngene$log2TPM1,
    exact = FALSE
  )

  diff1 <-
    median(
      early1$log2TPM1,
      na.rm = TRUE
    ) -
    median(
      ngene$log2TPM1,
      na.rm = TRUE
    )

  direction1 <- ifelse(
    diff1 > 0,
    "UP",
    "DOWN"
  )

  p <- prev[
    prev$Gene == gene,
  ]

  threshold <-
    p$Prevalence_threshold_TPM

  if (p$Direction == "UP") {

    aberrant1 <- sum(
      early1$TPM > threshold,
      na.rm = TRUE
    )

  } else {

    aberrant1 <- sum(
      early1$TPM < threshold,
      na.rm = TRUE
    )
  }

  primary_results[[gene]] <- data.frame(
    Gene = gene,

    Comparison =
      "Stage I tumor vs normal",

    Early_stage_definition =
      "Stage I",

    Early_tumor_n =
      nrow(early1),

    Normal_n =
      nrow(ngene),

    Early_tumor_median_TPM =
      median(
        early1$TPM,
        na.rm = TRUE
      ),

    Normal_median_TPM =
      median(
        ngene$TPM,
        na.rm = TRUE
      ),

    Median_log2TPM1_difference =
      diff1,

    Direction =
      direction1,

    Wilcoxon_P =
      wt1$p.value,

    Aberrant_early_n =
      aberrant1,

    Aberrant_early_percent =
      ifelse(
        nrow(early1) > 0,
        100 * aberrant1 / nrow(early1),
        NA_real_
      ),

    stringsAsFactors = FALSE
  )

  # ==========================================================
  # SENSITIVITY: Stage I + II vs normal
  # ==========================================================

  early12 <- tgene[
    tgene$Stage_group %in%
      c(
        "Stage I",
        "Stage II"
      ),
  ]

  wt12 <- wilcox.test(
    early12$log2TPM1,
    ngene$log2TPM1,
    exact = FALSE
  )

  diff12 <-
    median(
      early12$log2TPM1,
      na.rm = TRUE
    ) -
    median(
      ngene$log2TPM1,
      na.rm = TRUE
    )

  direction12 <- ifelse(
    diff12 > 0,
    "UP",
    "DOWN"
  )

  if (p$Direction == "UP") {

    aberrant12 <- sum(
      early12$TPM > threshold,
      na.rm = TRUE
    )

  } else {

    aberrant12 <- sum(
      early12$TPM < threshold,
      na.rm = TRUE
    )
  }

  sensitivity_results[[gene]] <- data.frame(
    Gene = gene,

    Comparison =
      "Stage I+II tumor vs normal",

    Early_stage_definition =
      "Stage I + Stage II",

    Early_tumor_n =
      nrow(early12),

    Normal_n =
      nrow(ngene),

    Early_tumor_median_TPM =
      median(
        early12$TPM,
        na.rm = TRUE
      ),

    Normal_median_TPM =
      median(
        ngene$TPM,
        na.rm = TRUE
      ),

    Median_log2TPM1_difference =
      diff12,

    Direction =
      direction12,

    Wilcoxon_P =
      wt12$p.value,

    Aberrant_early_n =
      aberrant12,

    Aberrant_early_percent =
      ifelse(
        nrow(early12) > 0,
        100 * aberrant12 / nrow(early12),
        NA_real_
      ),

    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Combine + FDR across genes
# ------------------------------------------------------------

primary <- do.call(
  rbind,
  primary_results
)

sensitivity <- do.call(
  rbind,
  sensitivity_results
)

primary$Wilcoxon_FDR <- p.adjust(
  primary$Wilcoxon_P,
  method = "BH"
)

sensitivity$Wilcoxon_FDR <- p.adjust(
  sensitivity$Wilcoxon_P,
  method = "BH"
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  primary,
  file.path(
    outdir,
    "17i_LUSC_stageI_vs_normal.csv"
  ),
  row.names = FALSE
)

write.csv(
  sensitivity,
  file.path(
    outdir,
    "17i_LUSC_stageI_II_vs_normal_sensitivity.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("LUSC EARLY-STAGE PRIMARY ANALYSIS\n")
cat("STAGE I VS NORMAL\n")
cat("========================================\n\n")

print(
  primary,
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("LUSC EARLY-STAGE SENSITIVITY\n")
cat("STAGE I + II VS NORMAL\n")
cat("========================================\n\n")

print(
  sensitivity,
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("17i COMPLETE\n")
cat("========================================\n")
