# ============================================================
# 17g — TCGA-LUSC survival analysis
# JAK3 + EPHB2
#
# Primary predictor:
#   continuous log2(TPM + 1)
#
# Primary adjusted model:
#   gene expression + age + sex + early/advanced stage
#
# Grade unavailable in TCGA-LUSC.
#
# Secondary analysis:
#   normal-derived aberrant-expression group
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

library(survival)

expr_file <-
  "results/17_TCGA_LUSC/17a_LUSC_JAK3_EPHB2_expression.csv"

surv_file <-
  "results/17_TCGA_LUSC/17e_LUSC_survival_metadata.csv"

prev_file <-
  "results/17_TCGA_LUSC/17b_LUSC_prevalence_summary.csv"

outdir <- "results/17_TCGA_LUSC"

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

patient$log2TPM1 <- log2(patient$TPM + 1)

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
      c("Stage I", "Stage II"),
    "Early",
    NA
  )
)

survdat$Advanced_stage <- factor(
  survdat$Advanced_stage,
  levels = c("Early", "Advanced")
)

survdat$gender <- factor(
  tolower(survdat$gender),
  levels = c("female", "male")
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
# Normal-derived aberrant-expression status
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

dat$Aberrant_label <- ifelse(
  dat$Aberrant,
  "Aberrant",
  "Reference"
)

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

get_ph <- function(fit, term) {

  pht <- tryCatch(
    cox.zph(
      fit,
      transform = "km"
    ),
    error = function(e) NULL
  )

  if (is.null(pht)) {
    return(
      c(
        Term_PH_P = NA,
        Global_PH_P = NA
      )
    )
  }

  term_p <- if (
    term %in% rownames(pht$table)
  ) {
    pht$table[term, "p"]
  } else {
    NA_real_
  }

  global_p <- if (
    "GLOBAL" %in% rownames(pht$table)
  ) {
    pht$table["GLOBAL", "p"]
  } else {
    NA_real_
  }

  c(
    Term_PH_P = term_p,
    Global_PH_P = global_p
  )
}

# ------------------------------------------------------------
# Analysis
# ------------------------------------------------------------

genes <- c(
  "JAK3",
  "EPHB2"
)

continuous_results <- list()
group_results <- list()
group_summary <- list()

for (gene in genes) {

  cat("\n========================================\n")
  cat(gene, "\n")
  cat("========================================\n")

  z <- dat[
    dat$gene_name == gene,
  ]

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

  # ==========================================================
  # Unadjusted continuous model
  # ==========================================================

  zu <- z[
    complete.cases(
      z$OS_days,
      z$survival_event,
      z$log2TPM1
    ),
  ]

  fit_unadj <- coxph(
    Surv(
      OS_days,
      survival_event
    ) ~ log2TPM1,
    data = zu
  )

  unadj <- extract_term(
    fit_unadj,
    "log2TPM1"
  )

  # ==========================================================
  # Primary adjusted model
  # ==========================================================

  za <- z[
    complete.cases(
      z$OS_days,
      z$survival_event,
      z$log2TPM1,
      z$age_at_diagnosis_years,
      z$gender,
      z$Advanced_stage
    ),
  ]

  fit_adj <- coxph(
    Surv(
      OS_days,
      survival_event
    ) ~
      log2TPM1 +
      age_at_diagnosis_years +
      gender +
      Advanced_stage,
    data = za,
    x = TRUE
  )

  adj <- extract_term(
    fit_adj,
    "log2TPM1"
  )

  ph <- get_ph(
    fit_adj,
    "log2TPM1"
  )

  continuous_results[[gene]] <- data.frame(
    Gene = gene,

    Survival_patient_n =
      nrow(z),

    Survival_deaths =
      sum(z$survival_event == 1),

    Continuous_model_n =
      nrow(zu),

    Continuous_HR =
      unadj["HR"],

    Continuous_lower95 =
      unadj["Lower95"],

    Continuous_upper95 =
      unadj["Upper95"],

    Continuous_P =
      unadj["P"],

    Adjusted_model_n =
      nrow(za),

    Adjusted_deaths =
      sum(za$survival_event == 1),

    Adjusted_HR =
      adj["HR"],

    Adjusted_lower95 =
      adj["Lower95"],

    Adjusted_upper95 =
      adj["Upper95"],

    Adjusted_P =
      adj["P"],

    Gene_PH_P =
      ph["Term_PH_P"],

    Global_PH_P =
      ph["Global_PH_P"],

    Adjustment =
      "Age + sex + early/advanced pathologic stage",

    Grade_adjustment =
      "Unavailable: TCGA-LUSC GDC tumor_grade field unpopulated",

    stringsAsFactors = FALSE
  )

  cat(
    "Adjusted model:",
    nrow(za),
    "patients /",
    sum(za$survival_event == 1),
    "deaths\n"
  )

  cat("\nAdjusted-model PH test:\n")

  print(
    cox.zph(
      fit_adj,
      transform = "km"
    )
  )

  # ==========================================================
  # Secondary normal-derived expression group
  # ==========================================================

  zg <- z[
    complete.cases(
      z$OS_days,
      z$survival_event,
      z$Aberrant
    ),
  ]

  zg$Aberrant_label <- factor(
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
    ) ~ Aberrant_label,
    data = zg
  )

  gr <- extract_term(
    fit_group,
    "Aberrant_labelAberrant"
  )

  sf <- survdiff(
    Surv(
      OS_days,
      survival_event
    ) ~ Aberrant_label,
    data = zg
  )

  logrank_p <-
    1 - pchisq(
      sf$chisq,
      df = length(sf$n) - 1
    )

  group_results[[gene]] <- data.frame(
    Gene = gene,

    Group_model_n =
      nrow(zg),

    Group_deaths =
      sum(zg$survival_event == 1),

    Aberrant_HR =
      gr["HR"],

    Aberrant_lower95 =
      gr["Lower95"],

    Aberrant_upper95 =
      gr["Upper95"],

    Cox_P =
      gr["P"],

    Logrank_P =
      logrank_p,

    stringsAsFactors = FALSE
  )

  # ==========================================================
  # Descriptive KM group summary
  # ==========================================================

  for (
    group in c(
      "Reference",
      "Aberrant"
    )
  ) {

    zz <- zg[
      zg$Aberrant_label == group,
    ]

    km <- survfit(
      Surv(
        OS_days,
        survival_event
      ) ~ 1,
      data = zz
    )

    med_days <- as.numeric(
      summary(km)$table["median"]
    )

    group_summary[[
      paste(gene, group, sep = "_")
    ]] <- data.frame(
      Gene = gene,
      Group = group,

      N =
        nrow(zz),

      Deaths =
        sum(zz$survival_event == 1),

      Median_OS_years =
        ifelse(
          is.na(med_days),
          NA_real_,
          med_days / 365.25
        ),

      stringsAsFactors = FALSE
    )
  }
}

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------

continuous_results <- do.call(
  rbind,
  continuous_results
)

group_results <- do.call(
  rbind,
  group_results
)

group_summary <- do.call(
  rbind,
  group_summary
)

# ------------------------------------------------------------
# BH correction across two genes within each test family
# ------------------------------------------------------------

continuous_results$Continuous_FDR <- p.adjust(
  continuous_results$Continuous_P,
  method = "BH"
)

continuous_results$Adjusted_FDR <- p.adjust(
  continuous_results$Adjusted_P,
  method = "BH"
)

group_results$Cox_FDR <- p.adjust(
  group_results$Cox_P,
  method = "BH"
)

group_results$Logrank_FDR <- p.adjust(
  group_results$Logrank_P,
  method = "BH"
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  continuous_results,
  file.path(
    outdir,
    "17g_LUSC_survival_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  group_results,
  file.path(
    outdir,
    "17g_LUSC_survival_group_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  group_summary,
  file.path(
    outdir,
    "17g_LUSC_survival_groups.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("LUSC CONTINUOUS SURVIVAL ASSOCIATION\n")
cat("========================================\n\n")

print(
  continuous_results,
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("LUSC NORMAL-DERIVED GROUP SURVIVAL\n")
cat("========================================\n\n")

print(
  group_results,
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("LUSC SURVIVAL GROUP SUMMARY\n")
cat("========================================\n\n")

print(
  group_summary,
  row.names = FALSE,
  digits = 6
)

cat("\n========================================\n")
cat("17g COMPLETE\n")
cat("========================================\n")
