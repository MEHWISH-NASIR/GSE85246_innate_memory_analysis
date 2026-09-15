# ============================================================
# 13k — TCGA-KIRC survival analysis
# JAK3 + EPHB2
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' is required.")
}

library(survival)

expr_file <-
  "results/13_TCGA_priority_genes/13i_KIRC_patient_expression_clinical.csv"

surv_file <-
  "results/13_TCGA_priority_genes/13j_KIRC_survival_metadata.csv"

prev_file <-
  "results/13_TCGA_priority_genes/13i_KIRC_directional_prevalence_summary.csv"

outdir <- "results/13_TCGA_priority_genes"

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

thr <- read.csv(
  prev_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dat <- merge(
  expr,
  surv[
    ,
    c(
      "patient_id",
      "survival_event",
      "OS_days",
      "OS_years",
      "survival_usable"
    )
  ],
  by = "patient_id",
  all.x = TRUE
)

dat$survival_usable <-
  tolower(as.character(dat$survival_usable)) == "true"

# ------------------------------------------------------------
# Clinical covariates
# ------------------------------------------------------------

dat$Advanced_stage <- ifelse(
  dat$Stage_group %in% c("Stage III", "Stage IV"),
  1,
  ifelse(
    dat$Stage_group %in% c("Stage I", "Stage II"),
    0,
    NA
  )
)

dat$High_grade <- ifelse(
  dat$Grade_group %in% c("G3", "G4"),
  1,
  ifelse(
    dat$Grade_group %in% c("G1", "G2"),
    0,
    NA
  )
)

dat$gender <- factor(dat$gender)

# ------------------------------------------------------------
# Assign normal-derived aberrant-expression status
# ------------------------------------------------------------

dat$Aberrant_group <- NA_character_

for (g in c("JAK3", "EPHB2")) {

  t <- thr[thr$Gene == g, ]

  idx <- dat$Gene == g

  if (t$Direction == "UP") {

    dat$Aberrant_group[idx] <-
      ifelse(
        dat$TPM[idx] >
          t$Prevalence_threshold_TPM,
        "Aberrant",
        "Reference"
      )

  } else {

    dat$Aberrant_group[idx] <-
      ifelse(
        dat$TPM[idx] <
          t$Prevalence_threshold_TPM,
        "Aberrant",
        "Reference"
      )
  }
}

dat$Aberrant_group <- factor(
  dat$Aberrant_group,
  levels = c("Reference", "Aberrant")
)

# ------------------------------------------------------------
# Analyse each gene
# ------------------------------------------------------------

genes <- c("JAK3", "EPHB2")

results <- list()
groups <- list()

for (g in genes) {

  z <- dat[
    dat$Gene == g &
    dat$survival_usable &
    !is.na(dat$OS_days) &
    !is.na(dat$survival_event),
  ]

  # ----------------------------------------------------------
  # Continuous-expression univariable Cox
  # ----------------------------------------------------------

  fit_uni <- coxph(
    Surv(OS_days, survival_event) ~ log2TPM1,
    data = z
  )

  su <- summary(fit_uni)

  uni_hr <-
    su$coefficients["log2TPM1", "exp(coef)"]

  uni_p <-
    su$coefficients["log2TPM1", "Pr(>|z|)"]

  uni_ci_low <-
    su$conf.int["log2TPM1", "lower .95"]

  uni_ci_high <-
    su$conf.int["log2TPM1", "upper .95"]

  # ----------------------------------------------------------
  # Multivariable Cox
  # age + sex + stage + grade
  # ----------------------------------------------------------

  za <- z[
    complete.cases(
      z[
        ,
        c(
          "OS_days",
          "survival_event",
          "log2TPM1",
          "age_at_diagnosis_years",
          "gender",
          "Advanced_stage",
          "High_grade"
        )
      ]
    ),
  ]

  za$gender <- droplevels(
    factor(za$gender)
  )

  fit_adj <- coxph(
    Surv(OS_days, survival_event) ~
      log2TPM1 +
      age_at_diagnosis_years +
      gender +
      Advanced_stage +
      High_grade,
    data = za
  )

  sa <- summary(fit_adj)

  adj_hr <-
    sa$coefficients["log2TPM1", "exp(coef)"]

  adj_p <-
    sa$coefficients["log2TPM1", "Pr(>|z|)"]

  adj_ci_low <-
    sa$conf.int["log2TPM1", "lower .95"]

  adj_ci_high <-
    sa$conf.int["log2TPM1", "upper .95"]

  zph <- cox.zph(fit_adj)

  ph_p <-
    zph$table[
      "log2TPM1",
      "p"
    ]

  # ----------------------------------------------------------
  # Aberrant-expression survival comparison
  # ----------------------------------------------------------

  zg <- z[
    !is.na(z$Aberrant_group),
  ]

  zg$Aberrant_group <-
    droplevels(
      factor(
        zg$Aberrant_group,
        levels = c(
          "Reference",
          "Aberrant"
        )
      )
    )

  fit_group <- coxph(
    Surv(OS_days, survival_event) ~
      Aberrant_group,
    data = zg
  )

  sg <- summary(fit_group)

  group_hr <-
    sg$coefficients[
      "Aberrant_groupAberrant",
      "exp(coef)"
    ]

  group_p <-
    sg$coefficients[
      "Aberrant_groupAberrant",
      "Pr(>|z|)"
    ]

  group_ci_low <-
    sg$conf.int[
      "Aberrant_groupAberrant",
      "lower .95"
    ]

  group_ci_high <-
    sg$conf.int[
      "Aberrant_groupAberrant",
      "upper .95"
    ]

  lr <- survdiff(
    Surv(OS_days, survival_event) ~
      Aberrant_group,
    data = zg
  )

  logrank_p <-
    pchisq(
      lr$chisq,
      df = length(lr$n) - 1,
      lower.tail = FALSE
    )

  # ----------------------------------------------------------
  # Main result row
  # ----------------------------------------------------------

  results[[g]] <- data.frame(
    Gene = g,

    Survival_patients =
      nrow(z),

    Deaths =
      sum(z$survival_event == 1),

    Adjusted_model_n =
      nrow(za),

    Adjusted_model_deaths =
      sum(za$survival_event == 1),

    Continuous_HR =
      uni_hr,

    Continuous_CI_low =
      uni_ci_low,

    Continuous_CI_high =
      uni_ci_high,

    Continuous_P =
      uni_p,

    Adjusted_HR =
      adj_hr,

    Adjusted_CI_low =
      adj_ci_low,

    Adjusted_CI_high =
      adj_ci_high,

    Adjusted_P =
      adj_p,

    PH_assumption_P =
      ph_p,

    Aberrant_status_HR =
      group_hr,

    Aberrant_status_CI_low =
      group_ci_low,

    Aberrant_status_CI_high =
      group_ci_high,

    Aberrant_status_P =
      group_p,

    Logrank_P =
      logrank_p,

    stringsAsFactors = FALSE
  )

  # ----------------------------------------------------------
  # Group counts
  # ----------------------------------------------------------

  for (
    gr in c(
      "Reference",
      "Aberrant"
    )
  ) {

    zz <- zg[
      zg$Aberrant_group == gr,
    ]

    groups[[
      paste(g, gr, sep = "_")
    ]] <- data.frame(
      Gene = g,
      Group = gr,
      N = nrow(zz),

      Deaths =
        sum(
          zz$survival_event == 1
        ),

      Median_OS_years =
        if (nrow(zz) > 0) {

          sf <- survfit(
            Surv(OS_days, survival_event) ~ 1,
            data = zz
          )

          as.numeric(
            summary(sf)$table["median"]
          ) / 365.25

        } else {

          NA_real_
        },

      stringsAsFactors = FALSE
    )
  }
}

res <- do.call(
  rbind,
  results
)

grp <- do.call(
  rbind,
  groups
)

# FDR across the two genes
res$Continuous_FDR <-
  p.adjust(
    res$Continuous_P,
    method = "BH"
  )

res$Adjusted_FDR <-
  p.adjust(
    res$Adjusted_P,
    method = "BH"
  )

res$Aberrant_status_FDR <-
  p.adjust(
    res$Aberrant_status_P,
    method = "BH"
  )

res$Logrank_FDR <-
  p.adjust(
    res$Logrank_P,
    method = "BH"
  )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  res,
  file.path(
    outdir,
    "13k_KIRC_survival_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  grp,
  file.path(
    outdir,
    "13k_KIRC_survival_groups.csv"
  ),
  row.names = FALSE
)

write.csv(
  dat,
  file.path(
    outdir,
    "13k_KIRC_survival_patient_data.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("KIRC SURVIVAL ASSOCIATION\n")
cat("========================================\n\n")

print(
  res,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("NORMAL-DERIVED EXPRESSION GROUPS\n")
cat("========================================\n\n")

print(
  grp,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("13k COMPLETE\n")
cat("========================================\n")
