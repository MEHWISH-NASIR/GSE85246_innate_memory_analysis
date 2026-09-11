# ============================================================
# 14d — TCGA-BLCA
# JAK3 + EPHB2 stage, grade and prevalence
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

expr_file <-
  "results/14_TCGA_BLCA/14a_BLCA_JAK3_EPHB2_expression.csv"

clin_file <-
  "results/14_TCGA_BLCA/14c_BLCA_clinical_metadata.csv"

prev_file <-
  "results/14_TCGA_BLCA/14b_BLCA_prevalence_summary.csv"

outdir <- "results/14_TCGA_BLCA"

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
# Patient-level tumor expression
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

patient$log2TPM1 <-
  log2(patient$TPM + 1)

tumor <- patient[
  patient$sample_type == "Primary Tumor",
]

cat("\nPatient-level BLCA tumors:\n")
print(table(tumor$gene_name))

# ------------------------------------------------------------
# Normalize pathological stage
# ------------------------------------------------------------

normalize_stage <- function(x) {

  x <- toupper(trimws(as.character(x)))

  out <- rep(
    NA_character_,
    length(x)
  )

  out[
    grepl("^STAGE 0", x)
  ] <- "Stage 0"

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

clin$Stage_number <- match(
  clin$Stage_group,
  c(
    "Stage 0",
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  )
) - 1

# ------------------------------------------------------------
# Normalize BLCA grade
# ------------------------------------------------------------

g <- toupper(
  trimws(
    as.character(
      clin$tumor_grade
    )
  )
)

clin$Grade_group <- ifelse(
  g == "LOW GRADE",
  "Low Grade",
  ifelse(
    g == "HIGH GRADE",
    "High Grade",
    NA_character_
  )
)

# ------------------------------------------------------------
# Merge
# ------------------------------------------------------------

dat <- merge(
  tumor,
  clin,
  by = "patient_id",
  all.x = TRUE
)

# ------------------------------------------------------------
# Assign normal-derived prevalence status
# ------------------------------------------------------------

dat$Aberrant_expression <- NA

for (gene in c("JAK3", "EPHB2")) {

  p <- prev[
    prev$Gene == gene,
  ]

  idx <- dat$gene_name == gene

  if (p$Direction == "UP") {

    dat$Aberrant_expression[idx] <-
      dat$TPM[idx] >
      p$Prevalence_threshold_TPM

  } else {

    dat$Aberrant_expression[idx] <-
      dat$TPM[idx] <
      p$Prevalence_threshold_TPM
  }
}

write.csv(
  dat,
  file.path(
    outdir,
    "14d_BLCA_patient_expression_clinical.csv"
  ),
  row.names = FALSE
)

genes <- c("JAK3", "EPHB2")

# ============================================================
# STAGE
# Formal trend uses Stage I-IV only.
# Stage 0 is retained in descriptive table.
# ============================================================

stage_stats <- list()
stage_summary <- list()

for (gene in genes) {

  z <- dat[
    dat$gene_name == gene &
    !is.na(dat$Stage_group),
  ]

  # Formal analysis: Stage I-IV only
  ztest <- z[
    z$Stage_group %in%
      c(
        "Stage I",
        "Stage II",
        "Stage III",
        "Stage IV"
      ),
  ]

  ztest$Stage_order <- match(
    ztest$Stage_group,
    c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )
  )

  kw <- kruskal.test(
    log2TPM1 ~ Stage_group,
    data = ztest
  )

  sp <- suppressWarnings(
    cor.test(
      ztest$Stage_order,
      ztest$log2TPM1,
      method = "spearman",
      exact = FALSE
    )
  )

  stage_stats[[gene]] <- data.frame(
    Gene = gene,

    Stage_patient_n =
      nrow(ztest),

    Kruskal_Wallis_P =
      kw$p.value,

    Spearman_rho =
      unname(sp$estimate),

    Stage_trend_P =
      sp$p.value,

    stringsAsFactors = FALSE
  )

  for (
    stage in c(
      "Stage 0",
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )
  ) {

    zz <- z[
      z$Stage_group == stage,
    ]

    stage_summary[[
      paste(gene, stage, sep = "_")
    ]] <- data.frame(
      Gene = gene,
      Stage = stage,

      N =
        nrow(zz),

      Median_TPM =
        if (nrow(zz) > 0)
          median(zz$TPM)
        else NA_real_,

      Mean_TPM =
        if (nrow(zz) > 0)
          mean(zz$TPM)
        else NA_real_,

      Aberrant_n =
        if (nrow(zz) > 0)
          sum(
            zz$Aberrant_expression,
            na.rm = TRUE
          )
        else NA_integer_,

      Aberrant_percent =
        if (nrow(zz) > 0)
          100 *
          mean(
            zz$Aberrant_expression,
            na.rm = TRUE
          )
        else NA_real_,

      stringsAsFactors = FALSE
    )
  }
}

stage_stats <- do.call(
  rbind,
  stage_stats
)

stage_stats$Kruskal_Wallis_FDR <-
  p.adjust(
    stage_stats$Kruskal_Wallis_P,
    method = "BH"
  )

stage_stats$Stage_trend_FDR <-
  p.adjust(
    stage_stats$Stage_trend_P,
    method = "BH"
  )

stage_summary <- do.call(
  rbind,
  stage_summary
)

# ============================================================
# GRADE
# BLCA has Low Grade vs High Grade
# ============================================================

grade_stats <- list()
grade_summary <- list()

for (gene in genes) {

  z <- dat[
    dat$gene_name == gene &
    !is.na(dat$Grade_group),
  ]

  low <- z$log2TPM1[
    z$Grade_group == "Low Grade"
  ]

  high <- z$log2TPM1[
    z$Grade_group == "High Grade"
  ]

  wt <- wilcox.test(
    high,
    low,
    exact = FALSE
  )

  grade_stats[[gene]] <- data.frame(
    Gene = gene,

    Grade_patient_n =
      nrow(z),

    Low_grade_n =
      sum(
        z$Grade_group == "Low Grade"
      ),

    High_grade_n =
      sum(
        z$Grade_group == "High Grade"
      ),

    Low_grade_median_TPM =
      median(
        z$TPM[
          z$Grade_group == "Low Grade"
        ]
      ),

    High_grade_median_TPM =
      median(
        z$TPM[
          z$Grade_group == "High Grade"
        ]
      ),

    Median_log2_difference_High_minus_Low =
      median(high) -
      median(low),

    Wilcoxon_P =
      wt$p.value,

    stringsAsFactors = FALSE
  )

  for (
    grade in c(
      "Low Grade",
      "High Grade"
    )
  ) {

    zz <- z[
      z$Grade_group == grade,
    ]

    grade_summary[[
      paste(gene, grade, sep = "_")
    ]] <- data.frame(
      Gene = gene,
      Grade = grade,

      N =
        nrow(zz),

      Median_TPM =
        median(
          zz$TPM
        ),

      Mean_TPM =
        mean(
          zz$TPM
        ),

      Aberrant_n =
        sum(
          zz$Aberrant_expression,
          na.rm = TRUE
        ),

      Aberrant_percent =
        100 *
        mean(
          zz$Aberrant_expression,
          na.rm = TRUE
        ),

      stringsAsFactors = FALSE
    )
  }
}

grade_stats <- do.call(
  rbind,
  grade_stats
)

grade_stats$Wilcoxon_FDR <-
  p.adjust(
    grade_stats$Wilcoxon_P,
    method = "BH"
  )

grade_summary <- do.call(
  rbind,
  grade_summary
)

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

write.csv(
  stage_stats,
  file.path(
    outdir,
    "14d_BLCA_stage_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  stage_summary,
  file.path(
    outdir,
    "14d_BLCA_stage_expression_prevalence.csv"
  ),
  row.names = FALSE
)

write.csv(
  grade_stats,
  file.path(
    outdir,
    "14d_BLCA_grade_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  grade_summary,
  file.path(
    outdir,
    "14d_BLCA_grade_expression_prevalence.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("BLCA STAGE ASSOCIATION\n")
cat("========================================\n\n")

print(
  stage_stats,
  row.names = FALSE,
  digits = 5
)

cat("\nSTAGE-SPECIFIC EXPRESSION / PREVALENCE\n\n")

print(
  stage_summary,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("BLCA GRADE ASSOCIATION\n")
cat("========================================\n\n")

print(
  grade_stats,
  row.names = FALSE,
  digits = 5
)

cat("\nGRADE-SPECIFIC EXPRESSION / PREVALENCE\n\n")

print(
  grade_summary,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("14d COMPLETE\n")
cat("========================================\n")
