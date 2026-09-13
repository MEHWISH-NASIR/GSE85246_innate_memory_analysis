# ============================================================
# 15d — TCGA-BRCA
# JAK3 + EPHB2 stage association and prevalence
# Grade is not assessable from GDC tumor_grade field
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

expr_file <-
  "results/15_TCGA_BRCA/15a_BRCA_JAK3_EPHB2_expression.csv"

clin_file <-
  "results/15_TCGA_BRCA/15c_BRCA_clinical_metadata.csv"

prev_file <-
  "results/15_TCGA_BRCA/15b_BRCA_prevalence_summary.csv"

outdir <- "results/15_TCGA_BRCA"

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

cat("\nPatient-level BRCA tumors:\n")
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

  # descriptive only
  out[
    grepl("^STAGE 0", x)
  ] <- "Stage 0"

  # order matters: IV before III before II before I
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

  # Stage X remains NA
  out
}

clin$Stage_group <-
  normalize_stage(
    clin$pathologic_stage
  )

# ------------------------------------------------------------
# Grade availability audit
# ------------------------------------------------------------

grade_valid <- !is.na(clin$tumor_grade) &
  clin$tumor_grade != "" &
  !toupper(clin$tumor_grade) %in%
    c(
      "NOT REPORTED",
      "UNKNOWN",
      "GX"
    )

grade_audit <- data.frame(
  Cancer = "TCGA-BRCA",
  Total_clinical_cases = nrow(clin),
  Usable_grade_cases = sum(grade_valid),
  Grade_analysis_status =
    ifelse(
      sum(grade_valid) >= 10,
      "Potentially assessable",
      "Not assessable"
    ),
  Reason =
    "GDC tumor_grade field is essentially unpopulated for TCGA-BRCA",
  stringsAsFactors = FALSE
)

write.csv(
  grade_audit,
  file.path(
    outdir,
    "15d_BRCA_grade_availability.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Merge clinical + expression
# ------------------------------------------------------------

dat <- merge(
  tumor,
  clin,
  by = "patient_id",
  all.x = TRUE
)

# ------------------------------------------------------------
# Assign normal-derived aberrant-expression status
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
    "15d_BRCA_patient_expression_clinical.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Stage association
# Formal testing uses Stage I-IV only
# Stage 0 retained descriptively
# ------------------------------------------------------------

genes <- c("JAK3", "EPHB2")

stage_stats <- list()
stage_summary <- list()

for (gene in genes) {

  z <- dat[
    dat$gene_name == gene &
    !is.na(dat$Stage_group),
  ]

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

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  stage_stats,
  file.path(
    outdir,
    "15d_BRCA_stage_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  stage_summary,
  file.path(
    outdir,
    "15d_BRCA_stage_expression_prevalence.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("BRCA STAGE ASSOCIATION\n")
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
cat("BRCA GRADE AVAILABILITY\n")
cat("========================================\n\n")

print(
  grade_audit,
  row.names = FALSE
)

cat("\n========================================\n")
cat("15d COMPLETE\n")
cat("========================================\n")
