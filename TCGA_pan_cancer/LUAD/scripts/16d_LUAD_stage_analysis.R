# ============================================================
# 16d — TCGA-LUAD
# JAK3 + EPHB2 stage association and prevalence
#
# Grade is not assessable:
# GDC tumor_grade is missing for all TCGA-LUAD cases.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

expr_file <-
  "results/16_TCGA_LUAD/16a_LUAD_JAK3_EPHB2_expression.csv"

clin_file <-
  "results/16_TCGA_LUAD/16c_LUAD_clinical_metadata.csv"

prev_file <-
  "results/16_TCGA_LUAD/16b_LUAD_prevalence_summary.csv"

outdir <- "results/16_TCGA_LUAD"

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
# Collapse multiple tumor aliquots by median
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

cat("\nPatient-level LUAD tumors:\n")
print(
  table(
    tumor$gene_name
  )
)

# ------------------------------------------------------------
# Normalize pathologic stage
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

  # Order matters
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
# Grade availability
# ------------------------------------------------------------

grade_valid <-
  !is.na(clin$tumor_grade) &
  clin$tumor_grade != "" &
  !toupper(clin$tumor_grade) %in%
    c(
      "NOT REPORTED",
      "UNKNOWN",
      "GX"
    )

grade_audit <- data.frame(
  Cancer = "TCGA-LUAD",

  Total_clinical_cases =
    nrow(clin),

  Usable_grade_cases =
    sum(grade_valid),

  Grade_analysis_status =
    ifelse(
      sum(grade_valid) >= 10,
      "Potentially assessable",
      "Not assessable"
    ),

  Reason =
    "GDC tumor_grade field is unpopulated for TCGA-LUAD",

  stringsAsFactors = FALSE
)

write.csv(
  grade_audit,
  file.path(
    outdir,
    "16d_LUAD_grade_availability.csv"
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
# Normal-derived aberrant expression
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
    "16d_LUAD_patient_expression_clinical.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Stage analysis
# Formal testing: Stage I-IV
# ------------------------------------------------------------

genes <- c(
  "JAK3",
  "EPHB2"
)

stage_stats <- list()
stage_summary <- list()

for (gene in genes) {

  z <- dat[
    dat$gene_name == gene &
    dat$Stage_group %in%
      c(
        "Stage I",
        "Stage II",
        "Stage III",
        "Stage IV"
      ),
  ]

  z$Stage_order <- match(
    z$Stage_group,
    c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )
  )

  # ------------------------------------------
  # Overall stage heterogeneity
  # ------------------------------------------

  kw <- kruskal.test(
    log2TPM1 ~ Stage_group,
    data = z
  )

  # ------------------------------------------
  # Ordered monotonic stage trend
  # ------------------------------------------

  sp <- suppressWarnings(
    cor.test(
      z$Stage_order,
      z$log2TPM1,
      method = "spearman",
      exact = FALSE
    )
  )

  stage_stats[[gene]] <- data.frame(
    Gene = gene,

    Stage_patient_n =
      nrow(z),

    Kruskal_Wallis_P =
      kw$p.value,

    Spearman_rho =
      unname(sp$estimate),

    Stage_trend_P =
      sp$p.value,

    stringsAsFactors = FALSE
  )

  # ------------------------------------------
  # Stage-specific descriptive values
  # ------------------------------------------

  for (
    stage in c(
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
      paste(
        gene,
        stage,
        sep = "_"
      )
    ]] <- data.frame(
      Gene = gene,
      Stage = stage,

      N =
        nrow(zz),

      Median_TPM =
        if (nrow(zz) > 0)
          median(
            zz$TPM,
            na.rm = TRUE
          )
        else
          NA_real_,

      Mean_TPM =
        if (nrow(zz) > 0)
          mean(
            zz$TPM,
            na.rm = TRUE
          )
        else
          NA_real_,

      Aberrant_n =
        if (nrow(zz) > 0)
          sum(
            zz$Aberrant_expression,
            na.rm = TRUE
          )
        else
          NA_integer_,

      Aberrant_percent =
        if (nrow(zz) > 0)
          100 *
          mean(
            zz$Aberrant_expression,
            na.rm = TRUE
          )
        else
          NA_real_,

      stringsAsFactors = FALSE
    )
  }
}

stage_stats <- do.call(
  rbind,
  stage_stats
)

stage_summary <- do.call(
  rbind,
  stage_summary
)

# ------------------------------------------------------------
# Multiple testing correction across the two genes
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  stage_stats,
  file.path(
    outdir,
    "16d_LUAD_stage_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  stage_summary,
  file.path(
    outdir,
    "16d_LUAD_stage_expression_prevalence.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("LUAD STAGE ASSOCIATION\n")
cat("========================================\n\n")

print(
  stage_stats,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("LUAD STAGE-SPECIFIC EXPRESSION / PREVALENCE\n")
cat("========================================\n\n")

print(
  stage_summary,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("LUAD GRADE AVAILABILITY\n")
cat("========================================\n\n")

print(
  grade_audit,
  row.names = FALSE
)

cat("\n========================================\n")
cat("16d COMPLETE\n")
cat("========================================\n")
