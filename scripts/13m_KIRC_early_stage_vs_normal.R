# ============================================================
# 13m — TCGA-KIRC early-stage tumor vs normal analysis
#
# Primary:
#   Stage I tumors vs normal
#
# Sensitivity:
#   Stage I + Stage II tumors vs normal
#
# Genes:
#   JAK3
#   EPHB2
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

expr_file <-
  "results/13_TCGA_priority_genes/13g_KIRC_JAK3_EPHB2_expression.csv"

clin_file <-
  "results/13_TCGA_priority_genes/13h_KIRC_clinical_metadata.csv"

prev_file <-
  "results/13_TCGA_priority_genes/13i_KIRC_directional_prevalence_summary.csv"

outdir <-
  "results/13_TCGA_priority_genes"

# ------------------------------------------------------------
# File checks
# ------------------------------------------------------------

if (!file.exists(expr_file)) {
  stop("Missing KIRC combined expression file.")
}

if (!file.exists(clin_file)) {
  stop("Missing KIRC clinical metadata.")
}

if (!file.exists(prev_file)) {
  stop("Missing KIRC directional prevalence summary.")
}

# ------------------------------------------------------------
# Read
# ------------------------------------------------------------

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
# Valid expression
# ------------------------------------------------------------

expr <- expr[
  expr$extraction_status == "OK" &
  !is.na(expr$TPM),
]

# ------------------------------------------------------------
# Patient-level expression
# Collapse duplicate aliquots using median
# ------------------------------------------------------------

patient <- aggregate(
  TPM ~ Gene + patient_id + sample_type,
  data = expr,
  FUN = median
)

patient$log2TPM1 <-
  log2(patient$TPM + 1)

# ------------------------------------------------------------
# Normalize stage
# ------------------------------------------------------------

normalize_stage <- function(x) {

  x <- toupper(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA_character_,
    length(x)
  )

  # Order is important
  out[
    grepl("^STAGE IV", x)
  ] <- "Stage IV"

  out[
    grepl("^STAGE III", x)
  ] <- "Stage III"

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
# Tumor and normal
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

genes <- c(
  "JAK3",
  "EPHB2"
)

primary_list <- list()
sensitivity_list <- list()

# ------------------------------------------------------------
# Analysis
# ------------------------------------------------------------

for (gene in genes) {

  tg <- tumor[
    tumor$Gene == gene,
  ]

  ng <- normal[
    normal$Gene == gene,
  ]

  threshold_info <- prev[
    prev$Gene == gene,
  ]

  if (nrow(threshold_info) != 1) {
    stop(
      paste(
        "Expected one prevalence row for",
        gene
      )
    )
  }

  threshold <-
    threshold_info$Prevalence_threshold_TPM

  prevalence_direction <-
    threshold_info$Direction

  # ==========================================================
  # PRIMARY: Stage I only
  # ==========================================================

  early1 <- tg[
    !is.na(tg$Stage_group) &
    tg$Stage_group == "Stage I",
  ]

  if (
    nrow(early1) >= 2 &&
    nrow(ng) >= 2
  ) {

    wt1 <- wilcox.test(
      early1$log2TPM1,
      ng$log2TPM1,
      exact = FALSE
    )

    p1 <- wt1$p.value

  } else {

    p1 <- NA_real_
  }

  early1_med <-
    if (nrow(early1) > 0)
      median(early1$TPM, na.rm = TRUE)
    else
      NA_real_

  normal_med <-
    if (nrow(ng) > 0)
      median(ng$TPM, na.rm = TRUE)
    else
      NA_real_

  diff1 <-
    if (
      nrow(early1) > 0 &&
      nrow(ng) > 0
    ) {
      median(
        early1$log2TPM1,
        na.rm = TRUE
      ) -
      median(
        ng$log2TPM1,
        na.rm = TRUE
      )
    } else {
      NA_real_
    }

  direction1 <- ifelse(
    is.na(diff1),
    NA_character_,
    ifelse(
      diff1 > 0,
      "UP",
      "DOWN"
    )
  )

  # Use the overall tumor-normal directional definition
  # for prevalence, preserving the existing KIRC framework.
  if (prevalence_direction == "UP") {

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

  primary_list[[gene]] <- data.frame(
    Gene = gene,

    Comparison =
      "Stage I tumor vs normal",

    Early_stage_definition =
      "Stage I",

    Early_tumor_n =
      nrow(early1),

    Normal_n =
      nrow(ng),

    Early_tumor_median_TPM =
      early1_med,

    Normal_median_TPM =
      normal_med,

    Median_log2TPM1_difference =
      diff1,

    Direction =
      direction1,

    Overall_prevalence_direction =
      prevalence_direction,

    Prevalence_threshold_TPM =
      threshold,

    Wilcoxon_P =
      p1,

    Aberrant_early_n =
      aberrant1,

    Aberrant_early_percent =
      ifelse(
        nrow(early1) > 0,
        100 *
          aberrant1 /
          nrow(early1),
        NA_real_
      ),

    stringsAsFactors = FALSE
  )

  # ==========================================================
  # SENSITIVITY: Stage I + II
  # ==========================================================

  early12 <- tg[
    !is.na(tg$Stage_group) &
    tg$Stage_group %in%
      c(
        "Stage I",
        "Stage II"
      ),
  ]

  if (
    nrow(early12) >= 2 &&
    nrow(ng) >= 2
  ) {

    wt12 <- wilcox.test(
      early12$log2TPM1,
      ng$log2TPM1,
      exact = FALSE
    )

    p12 <- wt12$p.value

  } else {

    p12 <- NA_real_
  }

  diff12 <-
    if (
      nrow(early12) > 0 &&
      nrow(ng) > 0
    ) {
      median(
        early12$log2TPM1,
        na.rm = TRUE
      ) -
      median(
        ng$log2TPM1,
        na.rm = TRUE
      )
    } else {
      NA_real_
    }

  direction12 <- ifelse(
    is.na(diff12),
    NA_character_,
    ifelse(
      diff12 > 0,
      "UP",
      "DOWN"
    )
  )

  if (prevalence_direction == "UP") {

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

  sensitivity_list[[gene]] <- data.frame(
    Gene = gene,

    Comparison =
      "Stage I+II tumor vs normal",

    Early_stage_definition =
      "Stage I + Stage II",

    Early_tumor_n =
      nrow(early12),

    Normal_n =
      nrow(ng),

    Early_tumor_median_TPM =
      if (nrow(early12) > 0)
        median(
          early12$TPM,
          na.rm = TRUE
        )
      else
        NA_real_,

    Normal_median_TPM =
      normal_med,

    Median_log2TPM1_difference =
      diff12,

    Direction =
      direction12,

    Overall_prevalence_direction =
      prevalence_direction,

    Prevalence_threshold_TPM =
      threshold,

    Wilcoxon_P =
      p12,

    Aberrant_early_n =
      aberrant12,

    Aberrant_early_percent =
      ifelse(
        nrow(early12) > 0,
        100 *
          aberrant12 /
          nrow(early12),
        NA_real_
      ),

    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------

primary <- do.call(
  rbind,
  primary_list
)

sensitivity <- do.call(
  rbind,
  sensitivity_list
)

# BH correction across the two genes
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
    "13m_KIRC_stageI_vs_normal.csv"
  ),
  row.names = FALSE
)

write.csv(
  sensitivity,
  file.path(
    outdir,
    "13m_KIRC_stageI_II_vs_normal_sensitivity.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("KIRC EARLY-STAGE PRIMARY ANALYSIS\n")
cat("STAGE I VS NORMAL\n")
cat("========================================\n\n")

print(
  primary,
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("KIRC EARLY-STAGE SENSITIVITY\n")
cat("STAGE I + II VS NORMAL\n")
cat("========================================\n\n")

print(
  sensitivity,
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("13m COMPLETE\n")
cat("========================================\n")
