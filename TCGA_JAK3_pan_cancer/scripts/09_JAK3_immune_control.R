# ============================================================
# 13o — TCGA-KIRC JAK3 immune/inflammatory control
#
# Questions:
# 1. Does JAK3 correlate with immune-cell markers?
# 2. Does JAK3 correlate with a multi-marker immune score?
# 3. Does tumor-vs-normal JAK3 elevation remain after
#    adjustment for immune content?
# 4. Does Stage-I-vs-normal JAK3 elevation remain after
#    adjustment for immune content?
#
# Primary immune controls:
#   PTPRC/CD45
#   multi-marker immune score
#
# Immune score:
#   mean of z-scored log2(TPM+1) for:
#   PTPRC, CD3D, CD3E, MS4A1, CD79A,
#   NKG7, LST1, TYROBP, FCER1G, CD68
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Files
# ------------------------------------------------------------

jak3_file <-
  "results/13_TCGA_JAK3/13c_KIRC_JAK3_expression.csv"

immune_file <-
  "results/13_TCGA_priority_genes/13n_KIRC_immune_markers_expression.csv"

clinical_file <-
  "results/13_TCGA_priority_genes/13h_KIRC_clinical_metadata.csv"

outdir <-
  "results/13_TCGA_priority_genes"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(jak3_file),
  file.exists(immune_file),
  file.exists(clinical_file)
)

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------

jak3 <- read.csv(
  jak3_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

immune <- read.csv(
  immune_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

clin <- read.csv(
  clinical_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# Valid records only
# ------------------------------------------------------------

jak3 <- jak3[
  jak3$extraction_status == "OK" &
  !is.na(jak3$TPM),
]

immune <- immune[
  immune$extraction_status == "OK" &
  !is.na(immune$TPM),
]

# ------------------------------------------------------------
# Patient-level JAK3
# Collapse duplicate aliquots using median
# ------------------------------------------------------------

jak3_patient <- aggregate(
  TPM ~ patient_id + sample_type,
  data = jak3,
  FUN = median
)

names(jak3_patient)[
  names(jak3_patient) == "TPM"
] <- "JAK3_TPM"

jak3_patient$JAK3_log2TPM1 <-
  log2(
    jak3_patient$JAK3_TPM + 1
  )

# ------------------------------------------------------------
# Patient-level immune markers
# ------------------------------------------------------------

immune_patient <- aggregate(
  TPM ~ patient_id + sample_type + gene_name,
  data = immune,
  FUN = median
)

markers <- c(
  "PTPRC",
  "CD3D",
  "CD3E",
  "MS4A1",
  "CD79A",
  "NKG7",
  "LST1",
  "TYROBP",
  "FCER1G",
  "CD68"
)

# ------------------------------------------------------------
# Convert immune table to wide format
# ------------------------------------------------------------

base_keys <- unique(
  immune_patient[
    ,
    c(
      "patient_id",
      "sample_type"
    )
  ]
)

immune_wide <- base_keys

for (g in markers) {

  z <- immune_patient[
    immune_patient$gene_name == g,
    c(
      "patient_id",
      "sample_type",
      "TPM"
    )
  ]

  names(z)[3] <- g

  immune_wide <- merge(
    immune_wide,
    z,
    by = c(
      "patient_id",
      "sample_type"
    ),
    all.x = TRUE
  )
}

# ------------------------------------------------------------
# Merge JAK3 + immune markers
# ------------------------------------------------------------

dat <- merge(
  jak3_patient,
  immune_wide,
  by = c(
    "patient_id",
    "sample_type"
  ),
  all = FALSE
)

# ------------------------------------------------------------
# Build log2 immune-marker matrix
# ------------------------------------------------------------

log_marker_names <- paste0(
  markers,
  "_log2TPM1"
)

for (i in seq_along(markers)) {

  dat[[log_marker_names[i]]] <- log2(
    dat[[markers[i]]] + 1
  )
}

# ------------------------------------------------------------
# Standardized multi-marker immune score
#
# Equal weight for each marker after z-standardization.
# JAK3 is NOT included in this score.
# ------------------------------------------------------------

zmat <- matrix(
  NA_real_,
  nrow = nrow(dat),
  ncol = length(markers)
)

colnames(zmat) <- markers

for (i in seq_along(markers)) {

  x <- dat[[log_marker_names[i]]]

  if (
    sum(!is.na(x)) > 1 &&
    sd(x, na.rm = TRUE) > 0
  ) {

    zmat[, i] <-
      as.numeric(
        scale(x)
      )
  }
}

dat$ImmuneScore <-
  rowMeans(
    zmat,
    na.rm = TRUE
  )

dat$PTPRC_log2TPM1 <-
  log2(
    dat$PTPRC + 1
  )

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

dat <- merge(
  dat,
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

# ------------------------------------------------------------
# Save patient-level analysis dataset
# ------------------------------------------------------------

write.csv(
  dat,
  file.path(
    outdir,
    "13o_KIRC_JAK3_immune_patient_data.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 1. JAK3 vs individual immune markers
# ============================================================

cor_results <- list()

contexts <- list(

  Primary_Tumor =
    dat[
      dat$sample_type == "Primary Tumor",
    ],

  Stage_I_Tumor =
    dat[
      dat$sample_type == "Primary Tumor" &
      !is.na(dat$Stage_group) &
      dat$Stage_group == "Stage I",
    ]
)

for (ctx in names(contexts)) {

  d <- contexts[[ctx]]

  for (g in markers) {

    x <- d$JAK3_log2TPM1

    y <- log2(
      d[[g]] + 1
    )

    ok <-
      complete.cases(
        x,
        y
      )

    if (sum(ok) >= 3) {

      ct <- suppressWarnings(
        cor.test(
          x[ok],
          y[ok],
          method = "spearman",
          exact = FALSE
        )
      )

      rho <-
        unname(
          ct$estimate
        )

      p <-
        ct$p.value

    } else {

      rho <- NA_real_
      p <- NA_real_
    }

    cor_results[[
      length(cor_results) + 1
    ]] <- data.frame(

      Context = ctx,

      Marker = g,

      N = sum(ok),

      Spearman_rho = rho,

      P = p,

      stringsAsFactors = FALSE
    )
  }
}

cor_tab <- do.call(
  rbind,
  cor_results
)

# FDR separately within each context
cor_tab$FDR <- NA_real_

for (ctx in unique(cor_tab$Context)) {

  idx <-
    cor_tab$Context == ctx

  cor_tab$FDR[idx] <-
    p.adjust(
      cor_tab$P[idx],
      method = "BH"
    )
}

write.csv(
  cor_tab,
  file.path(
    outdir,
    "13o_KIRC_JAK3_marker_correlations.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 2. JAK3 vs immune score
# ============================================================

score_results <- list()

for (ctx in names(contexts)) {

  d <- contexts[[ctx]]

  ok <-
    complete.cases(
      d$JAK3_log2TPM1,
      d$ImmuneScore
    )

  ct <- suppressWarnings(
    cor.test(
      d$JAK3_log2TPM1[ok],
      d$ImmuneScore[ok],
      method = "spearman",
      exact = FALSE
    )
  )

  score_results[[
    length(score_results) + 1
  ]] <- data.frame(

    Context = ctx,

    N = sum(ok),

    Spearman_rho =
      unname(
        ct$estimate
      ),

    P =
      ct$p.value,

    stringsAsFactors = FALSE
  )
}

score_tab <- do.call(
  rbind,
  score_results
)

score_tab$FDR <-
  p.adjust(
    score_tab$P,
    method = "BH"
  )

write.csv(
  score_tab,
  file.path(
    outdir,
    "13o_KIRC_JAK3_immune_score_correlations.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Helper for linear models
# ============================================================

extract_term <- function(
  fit,
  term,
  comparison,
  adjustment
) {

  sm <- summary(fit)$coefficients

  if (!(term %in% rownames(sm))) {

    return(
      data.frame(
        Comparison = comparison,
        Adjustment = adjustment,
        N = nobs(fit),
        Coefficient = NA_real_,
        SE = NA_real_,
        T = NA_real_,
        P = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  data.frame(

    Comparison =
      comparison,

    Adjustment =
      adjustment,

    N =
      nobs(fit),

    Coefficient =
      sm[
        term,
        "Estimate"
      ],

    SE =
      sm[
        term,
        "Std. Error"
      ],

    T =
      sm[
        term,
        "t value"
      ],

    P =
      sm[
        term,
        "Pr(>|t|)"
      ],

    stringsAsFactors = FALSE
  )
}

model_results <- list()

# ============================================================
# 3. ALL TUMORS VS NORMAL
# ============================================================

tn <- dat[
  dat$sample_type %in%
    c(
      "Primary Tumor",
      "Solid Tissue Normal"
    ),
]

tn$Tumor_vs_Normal <-
  ifelse(
    tn$sample_type == "Primary Tumor",
    1,
    0
  )

# Unadjusted linear model
fit <- lm(
  JAK3_log2TPM1 ~
    Tumor_vs_Normal,
  data = tn
)

model_results[[
  length(model_results) + 1
]] <- extract_term(
  fit,
  "Tumor_vs_Normal",
  "All primary tumors vs normal",
  "Unadjusted"
)

# PTPRC-adjusted
fit <- lm(
  JAK3_log2TPM1 ~
    Tumor_vs_Normal +
    PTPRC_log2TPM1,
  data = tn
)

model_results[[
  length(model_results) + 1
]] <- extract_term(
  fit,
  "Tumor_vs_Normal",
  "All primary tumors vs normal",
  "Adjusted for PTPRC/CD45"
)

# Immune-score adjusted
fit <- lm(
  JAK3_log2TPM1 ~
    Tumor_vs_Normal +
    ImmuneScore,
  data = tn
)

model_results[[
  length(model_results) + 1
]] <- extract_term(
  fit,
  "Tumor_vs_Normal",
  "All primary tumors vs normal",
  "Adjusted for multi-marker immune score"
)

# ============================================================
# 4. STAGE I TUMORS VS NORMAL
# ============================================================

early <- dat[
  (
    dat$sample_type == "Solid Tissue Normal"
  ) |
  (
    dat$sample_type == "Primary Tumor" &
    !is.na(dat$Stage_group) &
    dat$Stage_group == "Stage I"
  ),
]

early$StageI_vs_Normal <-
  ifelse(
    early$sample_type == "Primary Tumor",
    1,
    0
  )

# Unadjusted
fit <- lm(
  JAK3_log2TPM1 ~
    StageI_vs_Normal,
  data = early
)

model_results[[
  length(model_results) + 1
]] <- extract_term(
  fit,
  "StageI_vs_Normal",
  "Stage I tumors vs normal",
  "Unadjusted"
)

# PTPRC-adjusted
fit <- lm(
  JAK3_log2TPM1 ~
    StageI_vs_Normal +
    PTPRC_log2TPM1,
  data = early
)

model_results[[
  length(model_results) + 1
]] <- extract_term(
  fit,
  "StageI_vs_Normal",
  "Stage I tumors vs normal",
  "Adjusted for PTPRC/CD45"
)

# Immune-score adjusted
fit <- lm(
  JAK3_log2TPM1 ~
    StageI_vs_Normal +
    ImmuneScore,
  data = early
)

model_results[[
  length(model_results) + 1
]] <- extract_term(
  fit,
  "StageI_vs_Normal",
  "Stage I tumors vs normal",
  "Adjusted for multi-marker immune score"
)

model_tab <- do.call(
  rbind,
  model_results
)

# FDR across the six prespecified models
model_tab$FDR <-
  p.adjust(
    model_tab$P,
    method = "BH"
  )

write.csv(
  model_tab,
  file.path(
    outdir,
    "13o_KIRC_JAK3_immune_adjusted_models.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Console report
# ============================================================

cat("\n========================================\n")
cat("KIRC JAK3 IMMUNE CONTROL AUDIT\n")
cat("========================================\n\n")

cat(
  "Patient-level Primary Tumor:",
  sum(
    dat$sample_type == "Primary Tumor"
  ),
  "\n"
)

cat(
  "Patient-level Normal:",
  sum(
    dat$sample_type == "Solid Tissue Normal"
  ),
  "\n"
)

cat(
  "Stage I tumors:",
  sum(
    dat$sample_type == "Primary Tumor" &
    dat$Stage_group == "Stage I",
    na.rm = TRUE
  ),
  "\n"
)

cat("\n========================================\n")
cat("JAK3 vs MULTI-MARKER IMMUNE SCORE\n")
cat("========================================\n\n")

print(
  score_tab,
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("JAK3 vs INDIVIDUAL IMMUNE MARKERS\n")
cat("========================================\n\n")

print(
  cor_tab,
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("JAK3 IMMUNE-ADJUSTED MODELS\n")
cat("========================================\n\n")

print(
  model_tab,
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("13o COMPLETE\n")
cat("========================================\n")
