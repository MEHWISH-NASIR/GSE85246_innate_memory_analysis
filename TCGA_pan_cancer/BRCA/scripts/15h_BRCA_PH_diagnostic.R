# ============================================================
# 15h — TCGA-BRCA proportional hazards diagnostic
# Determine which adjusted-model covariate violates PH
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(survival)

expr_file <-
  "results/15_TCGA_BRCA/15a_BRCA_JAK3_EPHB2_expression.csv"

surv_file <-
  "results/15_TCGA_BRCA/15e_BRCA_survival_metadata.csv"

expr <- read.csv(
  expr_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

survdat <- read.csv(
  surv_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# Patient-level tumor expression
# ------------------------------------------------------------

expr <- expr[
  expr$extraction_status == "OK" &
  expr$sample_type == "Primary Tumor" &
  !is.na(expr$TPM),
]

patient <- aggregate(
  TPM ~ gene_name + patient_id,
  data = expr,
  FUN = median
)

patient$log2TPM1 <- log2(patient$TPM + 1)

# ------------------------------------------------------------
# Stage
# ------------------------------------------------------------

normalize_stage <- function(x) {

  x <- toupper(trimws(as.character(x)))

  out <- rep(NA_character_, length(x))

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

  out[
    grepl("^STAGE 0", x)
  ] <- "Stage 0"

  out
}

survdat$Stage_group <-
  normalize_stage(
    survdat$pathologic_stage
  )

survdat$Advanced_stage <- ifelse(
  survdat$Stage_group %in%
    c("Stage III", "Stage IV"),
  "Advanced",
  ifelse(
    survdat$Stage_group %in%
      c("Stage 0", "Stage I", "Stage II"),
    "Early",
    NA
  )
)

survdat$Advanced_stage <- factor(
  survdat$Advanced_stage,
  levels = c("Early", "Advanced")
)

# ------------------------------------------------------------
# Merge
# ------------------------------------------------------------

dat <- merge(
  patient,
  survdat,
  by = "patient_id",
  all.x = TRUE
)

dat <- dat[
  !is.na(dat$OS_days) &
  !is.na(dat$survival_event),
]

# ------------------------------------------------------------
# Diagnostic
# ------------------------------------------------------------

for (gene in c("JAK3", "EPHB2")) {

  cat("\n========================================\n")
  cat(gene, " ADJUSTED MODEL PH TEST\n")
  cat("========================================\n\n")

  z <- dat[
    dat$gene_name == gene &
    complete.cases(
      dat$OS_days,
      dat$survival_event,
      dat$log2TPM1,
      dat$age_at_diagnosis_years,
      dat$Advanced_stage
    ),
  ]

  fit <- coxph(
    Surv(
      OS_days,
      survival_event
    ) ~
      log2TPM1 +
      age_at_diagnosis_years +
      Advanced_stage,
    data = z,
    x = TRUE
  )

  ph <- cox.zph(
    fit,
    transform = "km"
  )

  print(ph)

  cat("\nModel coefficients:\n")
  print(
    summary(fit)$coefficients
  )
}

cat("\n========================================\n")
cat("15h COMPLETE\n")
cat("========================================\n")
