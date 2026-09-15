# ============================================================
# 15g — TCGA-BRCA survival analysis
# JAK3 + EPHB2
#
# Primary predictor:
#   continuous log2(TPM + 1)
#
# Adjusted model:
#   age + early/advanced stage
#
# Grade:
#   unavailable in TCGA-BRCA GDC clinical metadata
#
# Sex:
#   excluded from primary model because only 1 male death
#   occurred in the complete-case survival cohort, making the
#   sex coefficient poorly estimable.
#
# Secondary analysis:
#   normal-derived aberrant-expression group
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages(
    "survival",
    repos = "https://cloud.r-project.org"
  )
}

library(survival)

expr_file <-
  "results/15_TCGA_BRCA/15a_BRCA_JAK3_EPHB2_expression.csv"

surv_file <-
  "results/15_TCGA_BRCA/15e_BRCA_survival_metadata.csv"

prev_file <-
  "results/15_TCGA_BRCA/15b_BRCA_prevalence_summary.csv"

outdir <- "results/15_TCGA_BRCA"

# ------------------------------------------------------------
# Read
# ------------------------------------------------------------

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
# Merge expression and survival
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
# Assign normal-derived aberrant-expression group
# ------------------------------------------------------------

dat$Aberrant <- NA

for (gene in c("JAK3", "EPHB2")) {

  pp <- prev[
    prev$Gene == gene,
  ]

  idx <- dat$gene_name == gene

  if (pp$Direction == "UP") {

    dat$Aberrant[idx] <-
      dat$TPM[idx] >
      pp$Prevalence_threshold_TPM

  } else {

    dat$Aberrant[idx] <-
      dat$TPM[idx] <
      pp$Prevalence_threshold_TPM
  }
}

dat$Aberrant_label <- ifelse(
  dat$Aberrant,
  "Aberrant",
  "Reference"
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

extract_gene_cox <- function(fit, term) {

  s <- summary(fit)

  rn <- rownames(s$coefficients)

  i <- which(rn == term)

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

safe_ph_p <- function(fit, term) {

  z <- tryCatch(
    cox.zph(fit),
    error = function(e) NULL
  )

  if (is.null(z)) {
    return(NA_real_)
  }

  tab <- z$table

  if (!(term %in% rownames(tab))) {
    return(NA_real_)
  }

  as.numeric(
    tab[term, "p"]
  )
}

safe_global_ph_p <- function(fit) {

  z <- tryCatch(
    cox.zph(fit),
    error = function(e) NULL
  )

  if (is.null(z)) {
    return(NA_real_)
  }

  tab <- z$table

  if (!("GLOBAL" %in% rownames(tab))) {
    return(NA_real_)
  }

  as.numeric(
    tab["GLOBAL", "p"]
  )
}

# ------------------------------------------------------------
# Analysis
# ------------------------------------------------------------

genes <- c(
  "JAK3",
  "EPHB2"
)

results <- list()
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
  # 1. Continuous expression — unadjusted
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

  unadj <- extract_gene_cox(
    fit_unadj,
    "log2TPM1"
  )

  # ==========================================================
  # 2. Continuous expression — adjusted
  #
  # Primary adjusted BRCA model:
  # age + stage
  #
  # Sex excluded because the complete cohort contains only
  # one male death.
  # Grade unavailable.
  # ==========================================================

  za <- z[
    complete.cases(
      z$OS_days,
      z$survival_event,
      z$log2TPM1,
      z$age_at_diagnosis_years,
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
      Advanced_stage,
    data = za
  )

  adj <- extract_gene_cox(
    fit_adj,
    "log2TPM1"
  )

  ph_gene <- safe_ph_p(
    fit_adj,
    "log2TPM1"
  )

  ph_global <- safe_global_ph_p(
    fit_adj
  )

  results[[gene]] <- data.frame(
    Gene = gene,

    Survival_patient_n =
      nrow(z),

    Survival_deaths =
      sum(
        z$survival_event == 1
      ),

    Continuous_model_n =
      nrow(zu),

    Continuous_deaths =
      sum(
        zu$survival_event == 1
      ),

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
      sum(
        za$survival_event == 1
      ),

    Adjusted_HR =
      adj["HR"],

    Adjusted_lower95 =
      adj["Lower95"],

    Adjusted_upper95 =
      adj["Upper95"],

    Adjusted_P =
      adj["P"],

    Gene_PH_P =
      ph_gene,

    Global_PH_P =
      ph_global,

    Adjustment =
      "Age + early/advanced pathologic stage",

    Sex_adjustment =
      "Excluded: only one male death in complete-case cohort",

    Grade_adjustment =
      "Unavailable: TCGA-BRCA GDC tumor_grade field unpopulated",

    stringsAsFactors = FALSE
  )

  # ==========================================================
  # 3. Secondary normal-derived group analysis
  # ==========================================================

  zg <- z[
    complete.cases(
      z$OS_days,
      z$survival_event,
      z$Aberrant
    ),
  ]

  zg$Aberrant_label <- factor(
    zg$Aberrant_label,
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

  gr <- extract_gene_cox(
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
      sum(
        zg$survival_event == 1
      ),

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
  # Group descriptive survival
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
        sum(
          zz$survival_event == 1
        ),

      Median_OS_years =
        ifelse(
          is.na(med_days),
          NA_real_,
          med_days / 365.25
        ),

      stringsAsFactors = FALSE
    )
  }

  cat(
    "Adjusted model:",
    nrow(za),
    "patients /",
    sum(za$survival_event == 1),
    "deaths\n"
  )
}

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------

results <- do.call(
  rbind,
  results
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
# Multiple testing correction across the two genes
# within each survival test family
# ------------------------------------------------------------

results$Continuous_FDR <- p.adjust(
  results$Continuous_P,
  method = "BH"
)

results$Adjusted_FDR <- p.adjust(
  results$Adjusted_P,
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
  results,
  file.path(
    outdir,
    "15g_BRCA_survival_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  group_results,
  file.path(
    outdir,
    "15g_BRCA_survival_group_association.csv"
  ),
  row.names = FALSE
)

write.csv(
  group_summary,
  file.path(
    outdir,
    "15g_BRCA_survival_groups.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("BRCA CONTINUOUS SURVIVAL ASSOCIATION\n")
cat("========================================\n\n")

print(
  results,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("BRCA NORMAL-DERIVED GROUP SURVIVAL\n")
cat("========================================\n\n")

print(
  group_results,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("BRCA SURVIVAL GROUP SUMMARY\n")
cat("========================================\n\n")

print(
  group_summary,
  row.names = FALSE,
  digits = 5
)

cat("\n========================================\n")
cat("15g COMPLETE\n")
cat("========================================\n")
