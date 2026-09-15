# ============================================================
# 15i — TCGA-BRCA final stage-stratified survival model
#
# Reason:
# Advanced_stage violated the proportional hazards assumption
# in the ordinary adjusted Cox model.
#
# Final primary model:
#   Surv(OS) ~ log2(TPM+1) + age + strata(Advanced_stage)
#
# Grade unavailable.
# Sex excluded because only one male death occurred in the
# complete-case survival cohort.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(survival)

expr_file <-
  "results/15_TCGA_BRCA/15a_BRCA_JAK3_EPHB2_expression.csv"

surv_file <-
  "results/15_TCGA_BRCA/15e_BRCA_survival_metadata.csv"

prev_file <-
  "results/15_TCGA_BRCA/15b_BRCA_prevalence_summary.csv"

outdir <- "results/15_TCGA_BRCA"

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
  levels = c(
    "Early",
    "Advanced"
  )
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
  !is.na(dat$survival_event) &
  dat$OS_days >= 0,
]

# ------------------------------------------------------------
# Assign normal-derived group for secondary sensitivity
# ------------------------------------------------------------

dat$Aberrant <- NA

for (gene in c("JAK3", "EPHB2")) {

  p <- prev[
    prev$Gene == gene,
  ]

  idx <- dat$gene_name == gene

  if (p$Direction == "UP") {

    dat$Aberrant[idx] <-
      dat$TPM[idx] >
      p$Prevalence_threshold_TPM

  } else {

    dat$Aberrant[idx] <-
      dat$TPM[idx] <
      p$Prevalence_threshold_TPM
  }
}

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------

extract_term <- function(fit, term) {

  s <- summary(fit)

  i <- which(
    rownames(s$coefficients) == term
  )

  if (length(i) != 1) {

    return(
      c(
        HR = NA,
        Lower95 = NA,
        Upper95 = NA,
        P = NA
      )
    )
  }

  c(
    HR =
      s$conf.int[i, "exp(coef)"],

    Lower95 =
      s$conf.int[i, "lower .95"],

    Upper95 =
      s$conf.int[i, "upper .95"],

    P =
      s$coefficients[i, "Pr(>|z|)"]
  )
}

# ------------------------------------------------------------
# Continuous primary model
# ------------------------------------------------------------

continuous_results <- list()
group_results <- list()

for (gene in c("JAK3", "EPHB2")) {

  cat("\n========================================\n")
  cat(gene, "\n")
  cat("========================================\n")

  z <- dat[
    dat$gene_name == gene,
  ]

  za <- z[
    complete.cases(
      z$OS_days,
      z$survival_event,
      z$log2TPM1,
      z$age_at_diagnosis_years,
      z$Advanced_stage
    ),
  ]

  # ==========================================================
  # FINAL PRIMARY MODEL
  # ==========================================================

  fit <- coxph(
    Surv(
      OS_days,
      survival_event
    ) ~
      log2TPM1 +
      age_at_diagnosis_years +
      strata(Advanced_stage),
    data = za,
    x = TRUE
  )

  gene_result <- extract_term(
    fit,
    "log2TPM1"
  )

  ph <- cox.zph(
    fit,
    transform = "km"
  )

  gene_ph <- if (
    "log2TPM1" %in%
    rownames(ph$table)
  ) {
    ph$table[
      "log2TPM1",
      "p"
    ]
  } else {
    NA_real_
  }

  global_ph <- if (
    "GLOBAL" %in%
    rownames(ph$table)
  ) {
    ph$table[
      "GLOBAL",
      "p"
    ]
  } else {
    NA_real_
  }

  continuous_results[[gene]] <- data.frame(
    Gene = gene,

    Model_n =
      nrow(za),

    Deaths =
      sum(
        za$survival_event == 1
      ),

    Stratified_adjusted_HR =
      gene_result["HR"],

    Lower95 =
      gene_result["Lower95"],

    Upper95 =
      gene_result["Upper95"],

    P =
      gene_result["P"],

    Gene_PH_P =
      gene_ph,

    Global_PH_P =
      global_ph,

    Adjustment =
      paste0(
        "Age + strata(Early/Advanced stage); ",
        "stage stratified because PH assumption was violated"
      ),

    stringsAsFactors = FALSE
  )

  cat("\nFinal continuous model:\n")
  print(summary(fit)$coefficients)

  cat("\nPH test after stage stratification:\n")
  print(ph)

  # ==========================================================
  # SECONDARY GROUP SENSITIVITY
  #
  # Age-adjusted + stage-stratified version of the predefined
  # normal-derived aberrant-expression grouping.
  # ==========================================================

  zg <- z[
    complete.cases(
      z$OS_days,
      z$survival_event,
      z$Aberrant,
      z$age_at_diagnosis_years,
      z$Advanced_stage
    ),
  ]

  zg$Aberrant_group <- factor(
    ifelse(
      zg$Aberrant,
      "Aberrant",
      "Reference"
    ),
    levels = c(
      "Reference",
      "Aberrant"
    )
  )

  fit_group <- coxph(
    Surv(
      OS_days,
      survival_event
    ) ~
      Aberrant_group +
      age_at_diagnosis_years +
      strata(Advanced_stage),
    data = zg,
    x = TRUE
  )

  gr <- extract_term(
    fit_group,
    "Aberrant_groupAberrant"
  )

  phg <- cox.zph(
    fit_group,
    transform = "km"
  )

  group_ph <- if (
    "Aberrant_group" %in%
    rownames(phg$table)
  ) {
    phg$table[
      "Aberrant_group",
      "p"
    ]
  } else {
    NA_real_
  }

  group_results[[gene]] <- data.frame(
    Gene = gene,

    Model_n =
      nrow(zg),

    Deaths =
      sum(
        zg$survival_event == 1
      ),

    Aberrant_group_HR =
      gr["HR"],

    Lower95 =
      gr["Lower95"],

    Upper95 =
      gr["Upper95"],

    P =
      gr["P"],

    Group_PH_P =
      group_ph,

    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Combine + BH correction
# ------------------------------------------------------------

continuous_results <- do.call(
  rbind,
  continuous_results
)

group_results <- do.call(
  rbind,
  group_results
)

continuous_results$FDR <- p.adjust(
  continuous_results$P,
  method = "BH"
)

group_results$FDR <- p.adjust(
  group_results$P,
  method = "BH"
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  continuous_results,
  file.path(
    outdir,
    "15i_BRCA_final_stage_stratified_survival.csv"
  ),
  row.names = FALSE
)

write.csv(
  group_results,
  file.path(
    outdir,
    "15i_BRCA_group_stage_stratified_survival.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("FINAL BRCA STAGE-STRATIFIED SURVIVAL\n")
cat("========================================\n\n")

print(
  continuous_results,
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("SECONDARY GROUP SENSITIVITY\n")
cat("AGE-ADJUSTED + STAGE-STRATIFIED\n")
cat("========================================\n\n")

print(
  group_results,
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("15i COMPLETE\n")
cat("========================================\n")
