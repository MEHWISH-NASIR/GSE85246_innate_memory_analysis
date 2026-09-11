# ============================================================
# 14f — TCGA-BLCA survival analysis
# JAK3 + EPHB2
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' is required.")
}

library(survival)

expr_file <-
  "results/14_TCGA_BLCA/14d_BLCA_patient_expression_clinical.csv"

surv_file <-
  "results/14_TCGA_BLCA/14e_BLCA_survival_metadata.csv"

prev_file <-
  "results/14_TCGA_BLCA/14b_BLCA_prevalence_summary.csv"

outdir <- "results/14_TCGA_BLCA"

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

prev <- read.csv(
  prev_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# Merge survival
# ------------------------------------------------------------

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
  tolower(
    as.character(dat$survival_usable)
  ) == "true"

# ------------------------------------------------------------
# Clinical covariates
# ------------------------------------------------------------

# Stage 0/I/II versus III/IV
dat$Advanced_stage <- ifelse(
  dat$Stage_group %in%
    c("Stage III", "Stage IV"),
  1,
  ifelse(
    dat$Stage_group %in%
      c("Stage 0", "Stage I", "Stage II"),
    0,
    NA
  )
)

# BLCA grade: Low vs High
dat$High_grade <- ifelse(
  dat$Grade_group == "High Grade",
  1,
  ifelse(
    dat$Grade_group == "Low Grade",
    0,
    NA
  )
)

dat$gender <- factor(dat$gender)

# ------------------------------------------------------------
# Normal-derived aberrant-expression status
# Secondary analysis only
# ------------------------------------------------------------

dat$Aberrant_group <- NA_character_

for (g in c("JAK3", "EPHB2")) {

  p <- prev[
    prev$Gene == g,
  ]

  idx <- dat$gene_name == g

  if (p$Direction == "UP") {

    dat$Aberrant_group[idx] <- ifelse(
      dat$TPM[idx] >
        p$Prevalence_threshold_TPM,
      "Aberrant",
      "Reference"
    )

  } else {

    dat$Aberrant_group[idx] <- ifelse(
      dat$TPM[idx] <
        p$Prevalence_threshold_TPM,
      "Aberrant",
      "Reference"
    )
  }
}

dat$Aberrant_group <- factor(
  dat$Aberrant_group,
  levels = c(
    "Reference",
    "Aberrant"
  )
)

# ------------------------------------------------------------
# Analyse each gene
# ------------------------------------------------------------

genes <- c("JAK3", "EPHB2")

results <- list()
groups <- list()

for (g in genes) {

  z <- dat[
    dat$gene_name == g &
    dat$survival_usable &
    !is.na(dat$OS_days) &
    !is.na(dat$survival_event),
  ]

  # ==========================================================
  # Primary: continuous expression, univariable
  # ==========================================================

  fit_uni <- coxph(
    Surv(OS_days, survival_event) ~ log2TPM1,
    data = z
  )

  su <- summary(fit_uni)

  # ==========================================================
  # Primary: continuous expression, adjusted model
  # age + sex + stage + grade
  # ==========================================================

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
          "Advanced_stage"
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
      Advanced_stage,
    data = za
  )

  sa <- summary(fit_adj)

  # proportional hazards test
  ph <- cox.zph(fit_adj)

  # ==========================================================
  # Secondary: normal-derived aberrant group
  # ==========================================================

  zg <- z[
    !is.na(z$Aberrant_group),
  ]

  zg$Aberrant_group <- droplevels(
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

  lr <- survdiff(
    Surv(OS_days, survival_event) ~
      Aberrant_group,
    data = zg
  )

  logrank_p <- pchisq(
    lr$chisq,
    df = length(lr$n) - 1,
    lower.tail = FALSE
  )

  # ==========================================================
  # Result row
  # ==========================================================

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
      su$coefficients[
        "log2TPM1",
        "exp(coef)"
      ],

    Continuous_CI_low =
      su$conf.int[
        "log2TPM1",
        "lower .95"
      ],

    Continuous_CI_high =
      su$conf.int[
        "log2TPM1",
        "upper .95"
      ],

    Continuous_P =
      su$coefficients[
        "log2TPM1",
        "Pr(>|z|)"
      ],

    Adjusted_HR =
      sa$coefficients[
        "log2TPM1",
        "exp(coef)"
      ],

    Adjusted_CI_low =
      sa$conf.int[
        "log2TPM1",
        "lower .95"
      ],

    Adjusted_CI_high =
      sa$conf.int[
        "log2TPM1",
        "upper .95"
      ],

    Adjusted_P =
      sa$coefficients[
        "log2TPM1",
        "Pr(>|z|)"
      ],

    PH_assumption_P =
      ph$table[
        "log2TPM1",
        "p"
      ],

    Aberrant_status_HR =
      sg$coefficients[
        "Aberrant_groupAberrant",
        "exp(coef)"
      ],

    Aberrant_status_CI_low =
      sg$conf.int[
        "Aberrant_groupAberrant",
        "lower .95"
      ],

    Aberrant_status_CI_high =
      sg$conf.int[
        "Aberrant_groupAberrant",
        "upper .95"
      ],

    Aberrant_status_P =
      sg$coefficients[
        "Aberrant_groupAberrant",
        "Pr(>|z|)"
      ],

    Logrank_P =
      logrank_p,

    stringsAsFactors = FALSE
  )

  # ==========================================================
  # Survival groups
  # ==========================================================

  for (gr in c("Reference", "Aberrant")) {

    zz <- zg[
      zg$Aberrant_group == gr,
    ]

    med <- NA_real_

    if (nrow(zz) > 0) {

      sf <- survfit(
        Surv(OS_days, survival_event) ~ 1,
        data = zz
      )

      med_days <- as.numeric(
        summary(sf)$table["median"]
      )

      if (!is.na(med_days)) {
        med <- med_days / 365.25
      }
    }

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
        med,

      stringsAsFactors = FALSE
    )
  }
}

# ------------------------------------------------------------
# Combine + FDR across the two genes
# ------------------------------------------------------------

res <- do.call(
  rbind,
  results
)

grp <- do.call(
  rbind,
  groups
)

res$Continuous_FDR <- p.adjust(
  res$Continuous_P,
  method = "BH"
)

res$Adjusted_FDR <- p.adjust(
  res$Adjusted_P,
  method = "BH"
)

res$Aberrant_status_FDR <- p.adjust(
  res$Aberrant_status_P,
  method = "BH"
)

res$Logrank_FDR <- p.adjust(
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
    "14f_BLCA_survival_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  grp,
  file.path(
    outdir,
    "14f_BLCA_survival_groups.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("BLCA SURVIVAL ASSOCIATION\n")
cat("========================================\n\n")

print(
  res,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("BLCA NORMAL-DERIVED EXPRESSION GROUPS\n")
cat("========================================\n\n")

print(
  grp,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("14f COMPLETE\n")
cat("========================================\n")
