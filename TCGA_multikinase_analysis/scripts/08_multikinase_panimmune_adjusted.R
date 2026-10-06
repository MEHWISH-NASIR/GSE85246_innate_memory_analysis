
library(data.table)
setDTthreads(1)

root <- "TCGA_multikinase_analysis"

outdir <- file.path(root, "results/immune")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("1. Loading multikinase expression/clinical data...\n")

d <- fread(
  file.path(
    root,
    "results/expression",
    "03_multikinase_expression_clinical_complete.csv"
  ),
  nThread = 1
)

# ============================================================
# Helpers
# ============================================================

safe_cor <- function(x, y) {

  z <- complete.cases(x, y)

  if (sum(z) < 20) {
    return(c(rho = NA_real_, p = NA_real_))
  }

  out <- tryCatch(
    cor.test(
      x[z],
      y[z],
      method = "spearman",
      exact = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(c(rho = NA_real_, p = NA_real_))
  }

  c(
    rho = unname(out$estimate),
    p = out$p.value
  )
}

clean_stage <- function(x) {

  x <- toupper(trimws(as.character(x)))

  out <- rep(NA_real_, length(x))

  out[grepl("^STAGE IV", x)] <- 4
  out[is.na(out) & grepl("^STAGE III", x)] <- 3
  out[is.na(out) & grepl("^STAGE II", x)] <- 2
  out[is.na(out) & grepl("^STAGE I", x)] <- 1

  out
}

clean_grade <- function(project, grade) {

  g <- toupper(trimws(as.character(grade)))

  out <- rep(NA_real_, length(g))

  out[
    project == "TCGA-BLCA" &
    g == "LOW GRADE"
  ] <- 1

  out[
    project == "TCGA-BLCA" &
    g == "HIGH GRADE"
  ] <- 2

  out[
    project != "TCGA-BLCA" &
    g == "G1"
  ] <- 1

  out[
    project != "TCGA-BLCA" &
    g == "G2"
  ] <- 2

  out[
    project != "TCGA-BLCA" &
    g == "G3"
  ] <- 3

  out[
    project != "TCGA-BLCA" &
    g == "G4"
  ] <- 4

  out
}

first_nonmissing <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
    trimws(x) != ""
  ]

  if (length(x) == 0) {
    NA_character_
  } else {
    x[1]
  }
}

# ============================================================
# Primary malignant samples
# ============================================================

d[, primary_malignant :=
  (
    project == "TCGA-LAML" &
    sample_type ==
    "Primary Blood Derived Cancer - Peripheral Blood"
  ) |
  (
    project != "TCGA-LAML" &
    sample_type == "Primary Tumor"
  )
]

expr <- d[
  primary_malignant == TRUE
]

expr[, sample_barcode :=
  substr(sample.submitter_id, 1, 16)
]

expr[, cancer :=
  sub("^TCGA-", "", project)
]

cat(
  "Primary malignant expression rows:",
  nrow(expr),
  "\n"
)

# ============================================================
# Load PanImmune leukocyte fraction
# ============================================================

cat("2. Loading PanImmune leukocyte fractions...\n")

leuk <- fread(
  paste0(
    "TCGA_JAK3_pan_cancer/data/PanImmune/",
    "TCGA_all_leuk_estimate.masked.20170107.tsv"
  ),
  header = FALSE
)

setnames(
  leuk,
  c(
    "cancer",
    "panimmune_barcode",
    "leukocyte_fraction"
  )
)

leuk[, sample_barcode :=
  substr(
    panimmune_barcode,
    1,
    16
  )
]

# Collapse duplicate PanImmune records
leuk_clean <- leuk[
  ,
  .(
    leukocyte_fraction =
      mean(
        leukocyte_fraction,
        na.rm = TRUE
      ),

    panimmune_records = .N
  ),
  by = .(
    cancer,
    sample_barcode
  )
]

# ============================================================
# Match at sample level
# ============================================================

matched <- merge(
  expr,
  leuk_clean,
  by = c(
    "cancer",
    "sample_barcode"
  ),
  all = FALSE,
  sort = FALSE
)

cat(
  "Matched multikinase sample rows:",
  nrow(matched),
  "\n"
)

# ============================================================
# Collapse to patient level
# ============================================================

patient <- matched[
  ,
  .(
    TPM =
      mean(
        TPM,
        na.rm = TRUE
      ),

    log2TPM1 =
      mean(
        log2(TPM + 1),
        na.rm = TRUE
      ),

    leukocyte_fraction =
      mean(
        leukocyte_fraction,
        na.rm = TRUE
      ),

    stage =
      first_nonmissing(
        ajcc_pathologic_stage
      ),

    grade =
      first_nonmissing(
        tumor_grade
      ),

    expression_samples = .N
  ),
  by = .(
    gene,
    project,
    patient_id
  )
]

patient[, stage_num :=
  clean_stage(stage)
]

patient[, grade_num :=
  clean_grade(
    project,
    grade
  )
]

fwrite(
  patient,
  file.path(
    outdir,
    "08_multikinase_PanImmune_patient_data.csv"
  )
)

cat(
  "Matched patient-level rows:",
  nrow(patient),
  "\n"
)

# ============================================================
# 1. EXPRESSION vs LEUKOCYTE FRACTION
# ============================================================

groups <- unique(
  patient[, .(gene, project)]
)

cor_list <- vector(
  "list",
  nrow(groups)
)

for (i in seq_len(nrow(groups))) {

  g <- groups$gene[i]
  cancer <- groups$project[i]

  x <- patient[
    gene == g &
    project == cancer &
    !is.na(log2TPM1) &
    !is.na(leukocyte_fraction)
  ]

  ct <- safe_cor(
    x$log2TPM1,
    x$leukocyte_fraction
  )

  if (nrow(x) >= 20) {

    fit <- lm(
      log2TPM1 ~ leukocyte_fraction,
      data = x
    )

    s <- summary(fit)

    beta <- coef(s)[
      "leukocyte_fraction",
      "Estimate"
    ]

    p_lm <- coef(s)[
      "leukocyte_fraction",
      "Pr(>|t|)"
    ]

    r2 <- s$r.squared

  } else {

    beta <- NA_real_
    p_lm <- NA_real_
    r2 <- NA_real_
  }

  cor_list[[i]] <- data.table(

    gene = g,
    project = cancer,

    N = nrow(x),

    Median_leukocyte_fraction =
      median(
        x$leukocyte_fraction,
        na.rm = TRUE
      ),

    Spearman_rho =
      ct["rho"],

    Spearman_P =
      ct["p"],

    Linear_beta_LF =
      beta,

    Linear_P =
      p_lm,

    R_squared_LF =
      r2,

    Low_N =
      nrow(x) < 30
  )
}

cor_results <- rbindlist(
  cor_list,
  fill = TRUE
)

cor_results[, Spearman_FDR :=
  p.adjust(
    Spearman_P,
    method = "BH"
  ),
  by = gene
]

cor_results[, Linear_FDR :=
  p.adjust(
    Linear_P,
    method = "BH"
  ),
  by = gene
]

cor_results[, Global_Spearman_FDR :=
  p.adjust(
    Spearman_P,
    method = "BH"
  )
]

setorder(
  cor_results,
  gene,
  Spearman_FDR
)

fwrite(
  cor_results,
  file.path(
    outdir,
    "08_leukocyte_correlations.csv"
  )
)

cat("3. Leukocyte correlation analysis complete\n")

# ============================================================
# 2. STAGE AFTER LEUKOCYTE ADJUSTMENT
# ============================================================

stage_list <- list()

counter <- 1L

for (i in seq_len(nrow(groups))) {

  g <- groups$gene[i]
  cancer <- groups$project[i]

  x <- patient[
    gene == g &
    project == cancer &
    !is.na(stage_num) &
    !is.na(log2TPM1) &
    !is.na(leukocyte_fraction)
  ]

  if (nrow(x) == 0) {
    next
  }

  counts <- x[
    ,
    .N,
    by = stage_num
  ]

  eligible <- counts[
    N >= 5,
    stage_num
  ]

  x <- x[
    stage_num %in% eligible
  ]

  if (
    nrow(x) < 20 ||
    uniqueN(x$stage_num) < 2
  ) {
    next
  }

  fit0 <- lm(
    log2TPM1 ~ leukocyte_fraction,
    data = x
  )

  fit_cat <- lm(
    log2TPM1 ~
      leukocyte_fraction +
      factor(stage_num),
    data = x
  )

  a <- anova(
    fit0,
    fit_cat
  )

  stage_cat_p <-
    a$`Pr(>F)`[2]

  fit_trend <- lm(
    log2TPM1 ~
      leukocyte_fraction +
      stage_num,
    data = x
  )

  s <- summary(
    fit_trend
  )

  stage_beta <- coef(s)[
    "stage_num",
    "Estimate"
  ]

  stage_p <- coef(s)[
    "stage_num",
    "Pr(>|t|)"
  ]

  lf_beta <- coef(s)[
    "leukocyte_fraction",
    "Estimate"
  ]

  lf_p <- coef(s)[
    "leukocyte_fraction",
    "Pr(>|t|)"
  ]

  stage_list[[counter]] <-
    data.table(

      gene = g,
      project = cancer,

      N = nrow(x),

      Stage_groups_used =
        uniqueN(
          x$stage_num
        ),

      Stage_levels =
        paste(
          sort(
            unique(
              x$stage_num
            )
          ),
          collapse = ";"
        ),

      Adjusted_stage_categorical_P =
        stage_cat_p,

      Adjusted_stage_trend_beta =
        stage_beta,

      Adjusted_stage_trend_P =
        stage_p,

      Leukocyte_beta =
        lf_beta,

      Leukocyte_P =
        lf_p
    )

  counter <- counter + 1L
}

stage_results <- rbindlist(
  stage_list,
  fill = TRUE
)

stage_results[
  ,
  Adjusted_stage_categorical_FDR :=
    p.adjust(
      Adjusted_stage_categorical_P,
      method = "BH"
    ),
  by = gene
]

stage_results[
  ,
  Adjusted_stage_trend_FDR :=
    p.adjust(
      Adjusted_stage_trend_P,
      method = "BH"
    ),
  by = gene
]

stage_results[
  ,
  Global_adjusted_stage_trend_FDR :=
    p.adjust(
      Adjusted_stage_trend_P,
      method = "BH"
    )
]

fwrite(
  stage_results,
  file.path(
    outdir,
    "08_stage_leukocyte_adjusted.csv"
  )
)

cat("4. Leukocyte-adjusted stage analysis complete\n")

# ============================================================
# 3. GRADE AFTER LEUKOCYTE ADJUSTMENT
# ============================================================

grade_list <- list()

counter <- 1L

for (i in seq_len(nrow(groups))) {

  g <- groups$gene[i]
  cancer <- groups$project[i]

  x <- patient[
    gene == g &
    project == cancer &
    !is.na(grade_num) &
    !is.na(log2TPM1) &
    !is.na(leukocyte_fraction)
  ]

  if (nrow(x) == 0) {
    next
  }

  counts <- x[
    ,
    .N,
    by = grade_num
  ]

  eligible <- counts[
    N >= 5,
    grade_num
  ]

  x <- x[
    grade_num %in% eligible
  ]

  if (
    nrow(x) < 20 ||
    uniqueN(x$grade_num) < 2
  ) {
    next
  }

  fit0 <- lm(
    log2TPM1 ~ leukocyte_fraction,
    data = x
  )

  fit_cat <- lm(
    log2TPM1 ~
      leukocyte_fraction +
      factor(grade_num),
    data = x
  )

  a <- anova(
    fit0,
    fit_cat
  )

  grade_cat_p <-
    a$`Pr(>F)`[2]

  fit_trend <- lm(
    log2TPM1 ~
      leukocyte_fraction +
      grade_num,
    data = x
  )

  s <- summary(
    fit_trend
  )

  grade_beta <- coef(s)[
    "grade_num",
    "Estimate"
  ]

  grade_p <- coef(s)[
    "grade_num",
    "Pr(>|t|)"
  ]

  lf_beta <- coef(s)[
    "leukocyte_fraction",
    "Estimate"
  ]

  lf_p <- coef(s)[
    "leukocyte_fraction",
    "Pr(>|t|)"
  ]

  grade_list[[counter]] <-
    data.table(

      gene = g,
      project = cancer,

      N = nrow(x),

      Grade_groups_used =
        uniqueN(
          x$grade_num
        ),

      Grade_levels =
        paste(
          sort(
            unique(
              x$grade_num
            )
          ),
          collapse = ";"
        ),

      Adjusted_grade_categorical_P =
        grade_cat_p,

      Adjusted_grade_trend_beta =
        grade_beta,

      Adjusted_grade_trend_P =
        grade_p,

      Leukocyte_beta =
        lf_beta,

      Leukocyte_P =
        lf_p
    )

  counter <- counter + 1L
}

grade_results <- rbindlist(
  grade_list,
  fill = TRUE
)

grade_results[
  ,
  Adjusted_grade_categorical_FDR :=
    p.adjust(
      Adjusted_grade_categorical_P,
      method = "BH"
    ),
  by = gene
]

grade_results[
  ,
  Adjusted_grade_trend_FDR :=
    p.adjust(
      Adjusted_grade_trend_P,
      method = "BH"
    ),
  by = gene
]

grade_results[
  ,
  Global_adjusted_grade_trend_FDR :=
    p.adjust(
      Adjusted_grade_trend_P,
      method = "BH"
    )
]

fwrite(
  grade_results,
  file.path(
    outdir,
    "08_grade_leukocyte_adjusted.csv"
  )
)

cat("5. Leukocyte-adjusted grade analysis complete\n")

# ============================================================
# Gene-level comparison
# ============================================================

cor_summary <- cor_results[
  ,
  .(
    Cancers_with_LF = .N,

    Significant_LF_correlations =
      sum(
        Spearman_FDR < 0.05,
        na.rm = TRUE
      ),

    Positive_LF_correlations =
      sum(
        Spearman_FDR < 0.05 &
        Spearman_rho > 0,
        na.rm = TRUE
      ),

    Negative_LF_correlations =
      sum(
        Spearman_FDR < 0.05 &
        Spearman_rho < 0,
        na.rm = TRUE
      )
  ),
  by = gene
]

stage_summary <- stage_results[
  ,
  .(
    Adjusted_stage_cancers = .N,

    Adjusted_stage_categorical_sig =
      sum(
        Adjusted_stage_categorical_FDR < 0.05,
        na.rm = TRUE
      ),

    Adjusted_stage_trend_sig =
      sum(
        Adjusted_stage_trend_FDR < 0.05,
        na.rm = TRUE
      ),

    Positive_adjusted_stage_trends =
      sum(
        Adjusted_stage_trend_FDR < 0.05 &
        Adjusted_stage_trend_beta > 0,
        na.rm = TRUE
      ),

    Negative_adjusted_stage_trends =
      sum(
        Adjusted_stage_trend_FDR < 0.05 &
        Adjusted_stage_trend_beta < 0,
        na.rm = TRUE
      )
  ),
  by = gene
]

grade_summary <- grade_results[
  ,
  .(
    Adjusted_grade_cancers = .N,

    Adjusted_grade_categorical_sig =
      sum(
        Adjusted_grade_categorical_FDR < 0.05,
        na.rm = TRUE
      ),

    Adjusted_grade_trend_sig =
      sum(
        Adjusted_grade_trend_FDR < 0.05,
        na.rm = TRUE
      ),

    Positive_adjusted_grade_trends =
      sum(
        Adjusted_grade_trend_FDR < 0.05 &
        Adjusted_grade_trend_beta > 0,
        na.rm = TRUE
      ),

    Negative_adjusted_grade_trends =
      sum(
        Adjusted_grade_trend_FDR < 0.05 &
        Adjusted_grade_trend_beta < 0,
        na.rm = TRUE
      )
  ),
  by = gene
]

summary_gene <- merge(
  cor_summary,
  stage_summary,
  by = "gene",
  all = TRUE
)

summary_gene <- merge(
  summary_gene,
  grade_summary,
  by = "gene",
  all = TRUE
)

summary_gene[, Reference :=
  gene == "JAK3"
]

fwrite(
  summary_gene,
  file.path(
    outdir,
    "08_PanImmune_gene_comparison_summary.csv"
  )
)

# ============================================================
# Console
# ============================================================

cat("\n========================================\n")
cat("MULTIKINASE PANIMMUNE ADJUSTED ANALYSIS\n")
cat("========================================\n")

cat("\n===== GENE-LEVEL SUMMARY =====\n")
print(summary_gene)

cat("\n===== JAK3 LEUKOCYTE CORRELATIONS =====\n")
print(
  cor_results[
    gene == "JAK3"
  ]
)

cat("\n===== JAK3 ADJUSTED STAGE =====\n")
print(
  stage_results[
    gene == "JAK3"
  ]
)

cat("\n===== JAK3 ADJUSTED GRADE =====\n")
print(
  grade_results[
    gene == "JAK3"
  ]
)

cat("\n08 COMPLETE\n")

