# ============================================================
# 13i — TCGA-KIRC
# JAK3 + EPHB2
# Stage, grade and directional prevalence analysis
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

expr_file <-
  "results/13_TCGA_priority_genes/13g_KIRC_JAK3_EPHB2_expression.csv"

clin_file <-
  "results/13_TCGA_priority_genes/13h_KIRC_clinical_metadata.csv"

outdir <- "results/13_TCGA_priority_genes"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(expr_file)) stop("Missing expression file")
if (!file.exists(clin_file)) stop("Missing clinical file")

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

# ------------------------------------------------------------
# Keep valid expression
# ------------------------------------------------------------

expr <- expr[
  expr$extraction_status == "OK" &
  !is.na(expr$TPM),
]

# ------------------------------------------------------------
# Patient-level expression
# Collapse duplicate aliquots using median
# ------------------------------------------------------------

patient_expr <- aggregate(
  TPM ~ Gene + patient_id + sample_type,
  data = expr,
  FUN = median
)

patient_expr$log2TPM1 <-
  log2(patient_expr$TPM + 1)

tumor <- patient_expr[
  patient_expr$sample_type == "Primary Tumor",
]

normal <- patient_expr[
  patient_expr$sample_type == "Solid Tissue Normal",
]

cat("\nPatient-level tumors:\n")
print(table(tumor$Gene))

cat("\nPatient-level normals:\n")
print(table(normal$Gene))

# ------------------------------------------------------------
# Clean pathologic stage
# IA / IB -> I
# IIA / IIB -> II
# IIIA / IIIB / IIIC -> III
# IVA / IVB -> IV
# ------------------------------------------------------------

normalize_stage <- function(x) {

  x <- toupper(trimws(as.character(x)))

  out <- rep(
    NA_character_,
    length(x)
  )

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
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  )
)

# ------------------------------------------------------------
# Clean tumor grade
# Keep only G1-G4
# GX = unknown and excluded
# ------------------------------------------------------------

g <- toupper(
  trimws(
    as.character(
      clin$tumor_grade
    )
  )
)

clin$Grade_group <-
  ifelse(
    g %in% c("G1", "G2", "G3", "G4"),
    g,
    NA_character_
  )

clin$Grade_number <- match(
  clin$Grade_group,
  c("G1", "G2", "G3", "G4")
)

# ------------------------------------------------------------
# Merge clinical data with tumor expression
# ------------------------------------------------------------

tumor_clin <- merge(
  tumor,
  clin,
  by = "patient_id",
  all.x = TRUE
)

write.csv(
  tumor_clin,
  file.path(
    outdir,
    "13i_KIRC_patient_expression_clinical.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Determine expression direction and normal thresholds
# ------------------------------------------------------------

genes <- c(
  "JAK3",
  "EPHB2"
)

threshold_list <- list()

for (gene in genes) {

  t <- tumor$TPM[
    tumor$Gene == gene
  ]

  n <- normal$TPM[
    normal$Gene == gene
  ]

  tumor_median <-
    median(t, na.rm = TRUE)

  normal_median <-
    median(n, na.rm = TRUE)

  direction <-
    if (
      tumor_median > normal_median
    ) {
      "UP"
    } else {
      "DOWN"
    }

  p05 <- as.numeric(
    quantile(
      n,
      0.05,
      na.rm = TRUE
    )
  )

  p95 <- as.numeric(
    quantile(
      n,
      0.95,
      na.rm = TRUE
    )
  )

  if (direction == "UP") {

    abnormal_n <-
      sum(
        t > p95,
        na.rm = TRUE
      )

    threshold_used <- p95

    definition <-
      "Tumor TPM > normal 95th percentile"

  } else {

    abnormal_n <-
      sum(
        t < p05,
        na.rm = TRUE
      )

    threshold_used <- p05

    definition <-
      "Tumor TPM < normal 5th percentile"
  }

  threshold_list[[gene]] <-
    data.frame(
      Gene = gene,

      Direction = direction,

      Tumor_patient_n =
        length(t),

      Normal_patient_n =
        length(n),

      Tumor_median_TPM =
        tumor_median,

      Normal_median_TPM =
        normal_median,

      Normal_P05_TPM =
        p05,

      Normal_P95_TPM =
        p95,

      Prevalence_threshold_TPM =
        threshold_used,

      Prevalence_definition =
        definition,

      Aberrant_tumor_n =
        abnormal_n,

      Aberrant_tumor_percent =
        100 * abnormal_n /
        length(t),

      stringsAsFactors = FALSE
    )
}

thresholds <-
  do.call(
    rbind,
    threshold_list
  )

write.csv(
  thresholds,
  file.path(
    outdir,
    "13i_KIRC_directional_prevalence_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Assign aberrant-expression status to each tumor
# ------------------------------------------------------------

tumor_clin$Aberrant_expression <-
  NA

for (gene in genes) {

  z <- thresholds[
    thresholds$Gene == gene,
  ]

  idx <-
    tumor_clin$Gene == gene

  if (z$Direction == "UP") {

    tumor_clin$Aberrant_expression[idx] <-
      tumor_clin$TPM[idx] >
      z$Normal_P95_TPM

  } else {

    tumor_clin$Aberrant_expression[idx] <-
      tumor_clin$TPM[idx] <
      z$Normal_P05_TPM
  }
}

# ------------------------------------------------------------
# Stage association
# Kruskal-Wallis = any difference across stages
# Spearman = monotonic trend I -> IV
# ------------------------------------------------------------

stage_stats <- list()
stage_summary <- list()

for (gene in genes) {

  z <- tumor_clin[
    tumor_clin$Gene == gene &
    !is.na(tumor_clin$Stage_number),
  ]

  kw <- kruskal.test(
    log2TPM1 ~ Stage_group,
    data = z
  )

  sp <- suppressWarnings(
    cor.test(
      z$Stage_number,
      z$log2TPM1,
      method = "spearman",
      exact = FALSE
    )
  )

  stage_stats[[gene]] <-
    data.frame(
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
    ]] <-
      data.frame(
        Gene = gene,
        Stage = stage,
        N = nrow(zz),

        Median_TPM =
          if (nrow(zz) > 0)
            median(
              zz$TPM,
              na.rm = TRUE
            )
          else NA_real_,

        Mean_TPM =
          if (nrow(zz) > 0)
            mean(
              zz$TPM,
              na.rm = TRUE
            )
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

stage_stats <-
  do.call(
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

stage_summary <-
  do.call(
    rbind,
    stage_summary
  )

# ------------------------------------------------------------
# Grade association
# ------------------------------------------------------------

grade_stats <- list()
grade_summary <- list()

for (gene in genes) {

  z <- tumor_clin[
    tumor_clin$Gene == gene &
    !is.na(tumor_clin$Grade_number),
  ]

  kw <- kruskal.test(
    log2TPM1 ~ Grade_group,
    data = z
  )

  sp <- suppressWarnings(
    cor.test(
      z$Grade_number,
      z$log2TPM1,
      method = "spearman",
      exact = FALSE
    )
  )

  grade_stats[[gene]] <-
    data.frame(
      Gene = gene,

      Grade_patient_n =
        nrow(z),

      Kruskal_Wallis_P =
        kw$p.value,

      Spearman_rho =
        unname(sp$estimate),

      Grade_trend_P =
        sp$p.value,

      stringsAsFactors = FALSE
    )

  for (
    grade in c(
      "G1",
      "G2",
      "G3",
      "G4"
    )
  ) {

    zz <- z[
      z$Grade_group == grade,
    ]

    grade_summary[[
      paste(
        gene,
        grade,
        sep = "_"
      )
    ]] <-
      data.frame(
        Gene = gene,
        Grade = grade,
        N = nrow(zz),

        Median_TPM =
          if (nrow(zz) > 0)
            median(
              zz$TPM,
              na.rm = TRUE
            )
          else NA_real_,

        Mean_TPM =
          if (nrow(zz) > 0)
            mean(
              zz$TPM,
              na.rm = TRUE
            )
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

grade_stats <-
  do.call(
    rbind,
    grade_stats
  )

grade_stats$Kruskal_Wallis_FDR <-
  p.adjust(
    grade_stats$Kruskal_Wallis_P,
    method = "BH"
  )

grade_stats$Grade_trend_FDR <-
  p.adjust(
    grade_stats$Grade_trend_P,
    method = "BH"
  )

grade_summary <-
  do.call(
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
    "13i_KIRC_stage_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  stage_summary,
  file.path(
    outdir,
    "13i_KIRC_stage_expression_prevalence.csv"
  ),
  row.names = FALSE
)

write.csv(
  grade_stats,
  file.path(
    outdir,
    "13i_KIRC_grade_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  grade_summary,
  file.path(
    outdir,
    "13i_KIRC_grade_expression_prevalence.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("DIRECTIONAL PREVALENCE\n")
cat("========================================\n\n")

print(
  thresholds,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("STAGE ASSOCIATION\n")
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
cat("GRADE ASSOCIATION\n")
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
cat("13i COMPLETE\n")
cat("========================================\n")
