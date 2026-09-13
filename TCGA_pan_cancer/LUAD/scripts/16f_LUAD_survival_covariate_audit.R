# ============================================================
# 16f — TCGA-LUAD survival covariate audit
# JAK3 + EPHB2
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

expr_file <-
  "results/16_TCGA_LUAD/16a_LUAD_JAK3_EPHB2_expression.csv"

surv_file <-
  "results/16_TCGA_LUAD/16e_LUAD_survival_metadata.csv"

expr <- read.csv(
  expr_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

surv <- read.csv(
  surv_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# Patient-level primary-tumor expression
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

patient$log2TPM1 <-
  log2(patient$TPM + 1)

# ------------------------------------------------------------
# Stage grouping
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

surv$Stage_group <-
  normalize_stage(
    surv$pathologic_stage
  )

surv$Advanced_stage <- ifelse(
  surv$Stage_group %in%
    c("Stage III", "Stage IV"),
  "Advanced",
  ifelse(
    surv$Stage_group %in%
      c("Stage I", "Stage II"),
    "Early",
    NA
  )
)

# ------------------------------------------------------------
# Merge
# ------------------------------------------------------------

dat <- merge(
  patient,
  surv,
  by = "patient_id",
  all.x = TRUE
)

dat <- dat[
  !is.na(dat$OS_days) &
  !is.na(dat$survival_event),
]

cat("\n========================================\n")
cat("LUAD SURVIVAL COVARIATE AUDIT\n")
cat("========================================\n")

for (g in c("JAK3", "EPHB2")) {

  z <- dat[
    dat$gene_name == g,
  ]

  cat("\n----------------------------------------\n")
  cat(g, "\n")
  cat("----------------------------------------\n")

  cat(
    "Survival patients:",
    nrow(z),
    "\n"
  )

  cat(
    "Deaths:",
    sum(z$survival_event == 1),
    "\n"
  )

  cat("\nSex x event:\n")

  print(
    table(
      z$gender,
      z$survival_event,
      useNA = "ifany"
    )
  )

  cat("\nStage group x event:\n")

  print(
    table(
      z$Stage_group,
      z$survival_event,
      useNA = "ifany"
    )
  )

  cat("\nEarly/advanced stage x event:\n")

  print(
    table(
      z$Advanced_stage,
      z$survival_event,
      useNA = "ifany"
    )
  )

  complete <- complete.cases(
    z$OS_days,
    z$survival_event,
    z$log2TPM1,
    z$age_at_diagnosis_years,
    z$gender,
    z$Advanced_stage
  )

  zc <- z[
    complete,
  ]

  cat(
    "\nComplete cases for age + sex + stage model:",
    nrow(zc),
    "\n"
  )

  cat(
    "Deaths in complete model:",
    sum(zc$survival_event == 1),
    "\n"
  )

  cat("\nSex x event in complete model:\n")

  print(
    table(
      zc$gender,
      zc$survival_event
    )
  )

  cat("\nStage x event in complete model:\n")

  print(
    table(
      zc$Advanced_stage,
      zc$survival_event
    )
  )
}

cat("\n========================================\n")
cat("16f AUDIT COMPLETE\n")
cat("========================================\n")
