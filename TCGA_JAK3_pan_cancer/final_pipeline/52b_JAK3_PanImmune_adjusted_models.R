library(dplyr)

d <- read.csv(
  "TCGA_JAK3_pan_cancer/results/JAK3_PanImmune_patient_data_complete.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

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

  case_when(
    grepl("^STAGE IV", x)  ~ 4,
    grepl("^STAGE III", x) ~ 3,
    grepl("^STAGE II", x)  ~ 2,
    grepl("^STAGE I", x)   ~ 1,
    TRUE ~ NA_real_
  )
}

clean_grade <- function(project, grade) {

  g <- toupper(trimws(as.character(grade)))

  case_when(

    project == "TCGA-BLCA" &
      g == "LOW GRADE" ~ 1,

    project == "TCGA-BLCA" &
      g == "HIGH GRADE" ~ 2,

    project != "TCGA-BLCA" &
      g == "G1" ~ 1,

    project != "TCGA-BLCA" &
      g == "G2" ~ 2,

    project != "TCGA-BLCA" &
      g == "G3" ~ 3,

    project != "TCGA-BLCA" &
      g == "G4" ~ 4,

    TRUE ~ NA_real_
  )
}

d <- d %>%
  mutate(
    stage_num = clean_stage(stage),

    grade_num =
      clean_grade(project, grade)
  )

# ============================================================
# 1. JAK3 - LEUKOCYTE CORRELATION
# ============================================================

projects <- sort(unique(d$project))

cor_results <- lapply(
  projects,
  function(p) {

    x <- d %>%
      filter(
        project == p,
        !is.na(JAK3_log2TPM1),
        !is.na(leukocyte_fraction)
      )

    ctest <- safe_cor(
      x$JAK3_log2TPM1,
      x$leukocyte_fraction
    )

    if (nrow(x) >= 20) {

      fit <- lm(
        JAK3_log2TPM1 ~ leukocyte_fraction,
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

    data.frame(
      project = p,

      N = nrow(x),

      Median_leukocyte_fraction =
        median(
          x$leukocyte_fraction,
          na.rm = TRUE
        ),

      Spearman_rho =
        ctest["rho"],

      Spearman_P =
        ctest["p"],

      Linear_beta_LF =
        beta,

      Linear_P =
        p_lm,

      R_squared_LF =
        r2,

      Low_N =
        nrow(x) < 30,

      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

cor_results$Spearman_FDR <-
  p.adjust(
    cor_results$Spearman_P,
    method = "BH"
  )

cor_results$Linear_FDR <-
  p.adjust(
    cor_results$Linear_P,
    method = "BH"
  )

cor_results <- cor_results %>%
  arrange(Spearman_FDR)

# ============================================================
# 2. STAGE TREND AFTER LEUKOCYTE ADJUSTMENT
# ============================================================

stage_results <- list()

for (p in projects) {

  x <- d %>%
    filter(
      project == p,
      !is.na(stage_num),
      !is.na(JAK3_log2TPM1),
      !is.na(leukocyte_fraction)
    )

  # Remove stage groups with <5 patients
  eligible_groups <- x %>%
    count(stage_num) %>%
    filter(n >= 5) %>%
    pull(stage_num)

  x <- x %>%
    filter(stage_num %in% eligible_groups)

  if (
    nrow(x) < 20 ||
    n_distinct(x$stage_num) < 2
  ) {
    next
  }

  # Categorical stage effect adjusted for LF
  fit0 <- lm(
    JAK3_log2TPM1 ~ leukocyte_fraction,
    data = x
  )

  fit_cat <- lm(
    JAK3_log2TPM1 ~
      leukocyte_fraction +
      factor(stage_num),
    data = x
  )

  a <- anova(fit0, fit_cat)

  stage_cat_p <- a$`Pr(>F)`[2]

  # Ordered stage trend adjusted for LF
  fit_trend <- lm(
    JAK3_log2TPM1 ~
      leukocyte_fraction +
      stage_num,
    data = x
  )

  s <- summary(fit_trend)

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

  stage_results[[p]] <- data.frame(
    project = p,

    N = nrow(x),

    Stage_groups_used =
      n_distinct(x$stage_num),

    Stage_levels =
      paste(
        sort(unique(x$stage_num)),
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
      lf_p,

    stringsAsFactors = FALSE
  )
}

stage_results <- bind_rows(stage_results)

if (nrow(stage_results) > 0) {

  stage_results$Adjusted_stage_categorical_FDR <-
    p.adjust(
      stage_results$Adjusted_stage_categorical_P,
      method = "BH"
    )

  stage_results$Adjusted_stage_trend_FDR <-
    p.adjust(
      stage_results$Adjusted_stage_trend_P,
      method = "BH"
    )

  stage_results <- stage_results %>%
    arrange(Adjusted_stage_trend_FDR)
}

# ============================================================
# 3. GRADE TREND AFTER LEUKOCYTE ADJUSTMENT
# ============================================================

grade_results <- list()

for (p in projects) {

  x <- d %>%
    filter(
      project == p,
      !is.na(grade_num),
      !is.na(JAK3_log2TPM1),
      !is.na(leukocyte_fraction)
    )

  eligible_groups <- x %>%
    count(grade_num) %>%
    filter(n >= 5) %>%
    pull(grade_num)

  x <- x %>%
    filter(grade_num %in% eligible_groups)

  if (
    nrow(x) < 20 ||
    n_distinct(x$grade_num) < 2
  ) {
    next
  }

  fit0 <- lm(
    JAK3_log2TPM1 ~ leukocyte_fraction,
    data = x
  )

  fit_cat <- lm(
    JAK3_log2TPM1 ~
      leukocyte_fraction +
      factor(grade_num),
    data = x
  )

  a <- anova(fit0, fit_cat)

  grade_cat_p <- a$`Pr(>F)`[2]

  fit_trend <- lm(
    JAK3_log2TPM1 ~
      leukocyte_fraction +
      grade_num,
    data = x
  )

  s <- summary(fit_trend)

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

  grade_results[[p]] <- data.frame(
    project = p,

    N = nrow(x),

    Grade_groups_used =
      n_distinct(x$grade_num),

    Grade_levels =
      paste(
        sort(unique(x$grade_num)),
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
      lf_p,

    stringsAsFactors = FALSE
  )
}

grade_results <- bind_rows(grade_results)

if (nrow(grade_results) > 0) {

  grade_results$Adjusted_grade_categorical_FDR <-
    p.adjust(
      grade_results$Adjusted_grade_categorical_P,
      method = "BH"
    )

  grade_results$Adjusted_grade_trend_FDR <-
    p.adjust(
      grade_results$Adjusted_grade_trend_P,
      method = "BH"
    )

  grade_results <- grade_results %>%
    arrange(Adjusted_grade_trend_FDR)
}

# ============================================================
# SAVE
# ============================================================

write.csv(
  cor_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_PanImmune_leukocyte_correlations.csv",
  row.names = FALSE
)

write.csv(
  stage_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_PanImmune_stage_adjusted.csv",
  row.names = FALSE
)

write.csv(
  grade_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_PanImmune_grade_adjusted.csv",
  row.names = FALSE
)

# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("JAK3 PANIMMUNE LEUKOCYTE-ADJUSTED ANALYSIS\n")
cat("========================================\n")

cat(
  "\nCancers with leukocyte data:",
  nrow(cor_results),
  "\n"
)

cat(
  "JAK3-LF Spearman FDR < 0.05:",
  sum(
    cor_results$Spearman_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "JAK3-LF positive correlations:",
  sum(
    cor_results$Spearman_rho > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nStage cancers adjusted for LF:",
  nrow(stage_results),
  "\n"
)

cat(
  "Adjusted stage categorical FDR < 0.05:",
  sum(
    stage_results$Adjusted_stage_categorical_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Adjusted stage trend FDR < 0.05:",
  sum(
    stage_results$Adjusted_stage_trend_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nGrade cancers adjusted for LF:",
  nrow(grade_results),
  "\n"
)

cat(
  "Adjusted grade categorical FDR < 0.05:",
  sum(
    grade_results$Adjusted_grade_categorical_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Adjusted grade trend FDR < 0.05:",
  sum(
    grade_results$Adjusted_grade_trend_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat("\n===== TOP JAK3-LEUKOCYTE CORRELATIONS =====\n")

print(
  head(
    as.data.frame(cor_results),
    15
  ),
  row.names = FALSE
)

cat("\n===== LEUKOCYTE-ADJUSTED STAGE =====\n")

print(
  as.data.frame(stage_results),
  row.names = FALSE
)

cat("\n===== LEUKOCYTE-ADJUSTED GRADE =====\n")

print(
  as.data.frame(grade_results),
  row.names = FALSE
)

cat("\n52b PANIMMUNE ADJUSTED ANALYSIS COMPLETE\n")
