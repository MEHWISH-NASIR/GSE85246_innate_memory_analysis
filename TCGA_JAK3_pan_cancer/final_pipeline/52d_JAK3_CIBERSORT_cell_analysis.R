library(dplyr)
library(tidyr)

# ============================================================
# INPUTS
# ============================================================

cib <- read.csv(
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_quality_matched_samples.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

lf <- read.csv(
  "TCGA_JAK3_pan_cancer/results/JAK3_PanImmune_patient_data_complete.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

all_expr <- read.csv(
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_expression_clinical_complete_33projects.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

ciber_file <-
  "TCGA_JAK3_pan_cancer/data/PanImmune/TCGA.Kallisto.fullIDs.cibersort.relative.tsv"

ciber_raw <- read.delim(
  ciber_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cell_cols <- names(ciber_raw)[3:24]

# ============================================================
# HELPERS
# ============================================================

safe_spearman <- function(x, y) {

  z <- complete.cases(x, y)

  if (sum(z) < 10) {
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

# ============================================================
# QUALITY-FILTERED CIBERSORT:
# ONE ROW PER PATIENT
# ============================================================

patient_cib <- cib %>%
  group_by(
    project,
    cases.submitter_id
  ) %>%
  summarise(
    JAK3_TPM =
      mean(TPM, na.rm = TRUE),

    JAK3_log2TPM1 =
      mean(log2(TPM + 1), na.rm = TRUE),

    across(
      all_of(cell_cols),
      ~ mean(.x, na.rm = TRUE)
    ),

    CIBERSORT_records = n(),

    .groups = "drop"
  )

# ============================================================
# ADD LEUKOCYTE FRACTION
# ============================================================

lf_small <- lf %>%
  select(
    project,
    cases.submitter_id,
    leukocyte_fraction
  ) %>%
  distinct(
    project,
    cases.submitter_id,
    .keep_all = TRUE
  )

patient_cib <- patient_cib %>%
  left_join(
    lf_small,
    by = c(
      "project",
      "cases.submitter_id"
    )
  )

write.csv(
  patient_cib,
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_patient_level_complete.csv",
  row.names = FALSE
)

# ============================================================
# LONG FORMAT
# ============================================================

long <- patient_cib %>%
  pivot_longer(
    cols = all_of(cell_cols),
    names_to = "Cell_type",
    values_to = "Relative_fraction"
  ) %>%
  mutate(
    Tissue_fraction_proxy =
      leukocyte_fraction *
      Relative_fraction
  )

# ============================================================
# 1. QUALITY-FILTERED RELATIVE-FRACTION CORRELATIONS
# ============================================================

relative_results <- long %>%
  group_by(
    project,
    Cell_type
  ) %>%
  group_modify(
    ~ {
      test <- safe_spearman(
        .x$JAK3_log2TPM1,
        .x$Relative_fraction
      )

      data.frame(
        N = sum(
          complete.cases(
            .x$JAK3_log2TPM1,
            .x$Relative_fraction
          )
        ),

        rho = test["rho"],
        P = test["p"]
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    Low_N = N < 30
  )

# Global pan-cancer multiple testing
relative_results$Global_FDR <-
  p.adjust(
    relative_results$P,
    method = "BH"
  )

# Also provide cancer-specific 22-cell FDR
relative_results <- relative_results %>%
  group_by(project) %>%
  mutate(
    Within_cancer_FDR =
      p.adjust(
        P,
        method = "BH"
      )
  ) %>%
  ungroup()

# ============================================================
# 2. TISSUE-FRACTION PROXY
#    LF × CIBERSORT relative fraction
# ============================================================

tissue_results <- long %>%
  filter(
    !is.na(leukocyte_fraction)
  ) %>%
  group_by(
    project,
    Cell_type
  ) %>%
  group_modify(
    ~ {
      test <- safe_spearman(
        .x$JAK3_log2TPM1,
        .x$Tissue_fraction_proxy
      )

      data.frame(
        N = sum(
          complete.cases(
            .x$JAK3_log2TPM1,
            .x$Tissue_fraction_proxy
          )
        ),

        rho = test["rho"],
        P = test["p"]
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    Low_N = N < 30
  )

tissue_results$Global_FDR <-
  p.adjust(
    tissue_results$P,
    method = "BH"
  )

tissue_results <- tissue_results %>%
  group_by(project) %>%
  mutate(
    Within_cancer_FDR =
      p.adjust(
        P,
        method = "BH"
      )
  ) %>%
  ungroup()

# ============================================================
# 3. CELL COMPOSITION AFTER TOTAL LF ADJUSTMENT
#
# JAK3 ~ leukocyte fraction + relative cell fraction
# One immune cell at a time to avoid compositional
# multicollinearity from including all 22 together.
# ============================================================

adjusted_results <- long %>%
  filter(
    !is.na(leukocyte_fraction)
  ) %>%
  group_by(
    project,
    Cell_type
  ) %>%
  group_modify(
    ~ {

      x <- .x %>%
        filter(
          complete.cases(
            JAK3_log2TPM1,
            leukocyte_fraction,
            Relative_fraction
          )
        )

      if (nrow(x) < 30 ||
          sd(x$Relative_fraction) == 0) {

        return(
          data.frame(
            N = nrow(x),
            Cell_beta = NA_real_,
            Cell_P = NA_real_,
            LF_beta = NA_real_,
            LF_P = NA_real_
          )
        )
      }

      fit <- tryCatch(
        lm(
          JAK3_log2TPM1 ~
            leukocyte_fraction +
            Relative_fraction,
          data = x
        ),
        error = function(e) NULL
      )

      if (is.null(fit)) {

        return(
          data.frame(
            N = nrow(x),
            Cell_beta = NA_real_,
            Cell_P = NA_real_,
            LF_beta = NA_real_,
            LF_P = NA_real_
          )
        )
      }

      s <- summary(fit)$coefficients

      data.frame(
        N = nrow(x),

        Cell_beta =
          s[
            "Relative_fraction",
            "Estimate"
          ],

        Cell_P =
          s[
            "Relative_fraction",
            "Pr(>|t|)"
          ],

        LF_beta =
          s[
            "leukocyte_fraction",
            "Estimate"
          ],

        LF_P =
          s[
            "leukocyte_fraction",
            "Pr(>|t|)"
          ]
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    Low_N = N < 30
  )

adjusted_results$Global_Cell_FDR <-
  p.adjust(
    adjusted_results$Cell_P,
    method = "BH"
  )

adjusted_results <- adjusted_results %>%
  group_by(project) %>%
  mutate(
    Within_cancer_Cell_FDR =
      p.adjust(
        Cell_P,
        method = "BH"
      )
  ) %>%
  ungroup()

# ============================================================
# 4. ALL-CIBERSORT SAMPLES SENSITIVITY ANALYSIS
# ============================================================

all_ciber <- ciber_raw %>%
  mutate(
    SampleID_hyphen =
      gsub("\\.", "-", SampleID),

    sample_barcode =
      substr(
        SampleID_hyphen,
        1,
        16
      )
  ) %>%
  arrange(
    CancerType,
    sample_barcode,
    P.value,
    desc(Correlation),
    RMSE
  ) %>%
  group_by(
    CancerType,
    sample_barcode
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

expr <- all_expr %>%
  mutate(
    primary_malignant = case_when(

      project == "TCGA-LAML" &
        sample_type ==
        "Primary Blood Derived Cancer - Peripheral Blood" ~ TRUE,

      project != "TCGA-LAML" &
        sample_type == "Primary Tumor" ~ TRUE,

      TRUE ~ FALSE
    ),

    cancer =
      sub("^TCGA-", "", project),

    sample_barcode =
      substr(
        sample.submitter_id,
        1,
        16
      )
  ) %>%
  filter(primary_malignant)

all_match <- expr %>%
  inner_join(
    all_ciber,
    by = c(
      "cancer" = "CancerType",
      "sample_barcode"
    )
  ) %>%
  group_by(
    project,
    cases.submitter_id
  ) %>%
  summarise(
    JAK3_log2TPM1 =
      mean(
        log2(TPM + 1),
        na.rm = TRUE
      ),

    across(
      all_of(cell_cols),
      ~ mean(.x, na.rm = TRUE)
    ),

    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(cell_cols),
    names_to = "Cell_type",
    values_to = "Relative_fraction"
  )

sensitivity_results <- all_match %>%
  group_by(
    project,
    Cell_type
  ) %>%
  group_modify(
    ~ {
      test <- safe_spearman(
        .x$JAK3_log2TPM1,
        .x$Relative_fraction
      )

      data.frame(
        N = sum(
          complete.cases(
            .x$JAK3_log2TPM1,
            .x$Relative_fraction
          )
        ),

        rho = test["rho"],
        P = test["p"]
      )
    }
  ) %>%
  ungroup()

sensitivity_results$Global_FDR <-
  p.adjust(
    sensitivity_results$P,
    method = "BH"
  )

# ============================================================
# SAVE
# ============================================================

write.csv(
  relative_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_relative_correlations.csv",
  row.names = FALSE
)

write.csv(
  tissue_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_tissue_fraction_correlations.csv",
  row.names = FALSE
)

write.csv(
  adjusted_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_LF_adjusted_cell_models.csv",
  row.names = FALSE
)

write.csv(
  sensitivity_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_all_samples_sensitivity.csv",
  row.names = FALSE
)

# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("JAK3 PAN-CANCER CIBERSORT ANALYSIS\n")
cat("========================================\n")

cat(
  "\nQuality-filtered cancers:",
  n_distinct(relative_results$project),
  "\n"
)

cat(
  "Cancer-cell tests:",
  nrow(relative_results),
  "\n"
)

cat(
  "Relative-fraction global FDR < 0.05:",
  sum(
    relative_results$Global_FDR < 0.05 &
    !relative_results$Low_N,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Tissue-fraction global FDR < 0.05:",
  sum(
    tissue_results$Global_FDR < 0.05 &
    !tissue_results$Low_N,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "LF-adjusted cell-composition global FDR < 0.05:",
  sum(
    adjusted_results$Global_Cell_FDR < 0.05 &
    !adjusted_results$Low_N,
    na.rm = TRUE
  ),
  "\n"
)

cat("\n===== TOP TISSUE-FRACTION ASSOCIATIONS =====\n")

print(
  head(
    tissue_results %>%
      filter(!Low_N) %>%
      arrange(Global_FDR),
    25
  ),
  row.names = FALSE
)

cat("\n===== TOP CELL EFFECTS AFTER LF ADJUSTMENT =====\n")

print(
  head(
    adjusted_results %>%
      filter(!Low_N) %>%
      arrange(Global_Cell_FDR),
    25
  ),
  row.names = FALSE
)

cat("\n52d CIBERSORT CELL ANALYSIS COMPLETE\n")
