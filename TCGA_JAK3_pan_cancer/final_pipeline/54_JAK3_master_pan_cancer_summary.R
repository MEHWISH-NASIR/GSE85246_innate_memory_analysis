library(dplyr)

# ============================================================
# 54 — JAK3 MASTER PAN-CANCER SUMMARY
#
# Integrates the validated canonical JAK3 analyses:
#
# - analysis availability / sample sizes
# - tumor vs normal
# - matched tumor-normal
# - Stage-I vs normal
# - aberrant-expression prevalence
# - stage statistics
# - grade statistics
# - TCGA-CDR survival
# - PanImmune leukocyte fraction
# - leukocyte-adjusted stage / grade
# - CIBERSORT robustness
# - PTPRC / immune-score correlations
# - immune-adjusted tumor-normal models
# - immune-adjusted Stage-I-normal models
#
# One row per TCGA project.
# No overall ranking or arbitrary score is generated.
# ============================================================

root <- "TCGA_JAK3_pan_cancer/results"

outdir <- file.path(
  root,
  "master"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# HELPERS
# ============================================================

read_required <- function(path) {

  if (!file.exists(path)) {
    stop(
      paste(
        "Required file missing:",
        path
      )
    )
  }

  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

check_project_unique <- function(dat, label) {

  if (!("project" %in% names(dat))) {
    stop(
      paste(
        label,
        "has no project column."
      )
    )
  }

  dup <- dat$project[
    duplicated(dat$project)
  ]

  if (length(dup) > 0) {
    stop(
      paste(
        "Duplicate projects in",
        label,
        ":",
        paste(
          unique(dup),
          collapse = ", "
        )
      )
    )
  }

  invisible(TRUE)
}

# Helper for immune-adjusted model extraction
extract_model <- function(
  dat,
  comparison,
  adjustment,
  prefix
) {

  z <- dat[
    dat$Comparison == comparison &
    dat$Adjustment == adjustment,
    c(
      "project",
      "N",
      "Tumor_n",
      "Normal_n",
      "Coefficient",
      "P",
      "Within_cancer_FDR",
      "Model_family_FDR",
      "Global_FDR",
      "Low_N"
    )
  ]

  check_project_unique(
    z,
    paste(
      comparison,
      adjustment
    )
  )

  names(z)[-1] <-
    paste0(
      prefix,
      names(z)[-1]
    )

  z
}

# Helper for selected immune correlations
extract_correlation <- function(
  dat,
  context,
  feature,
  prefix
) {

  z <- dat[
    dat$Context == context &
    dat$Feature == feature,
    c(
      "project",
      "N",
      "Spearman_rho",
      "P",
      "Within_cancer_context_FDR",
      "Feature_family_FDR",
      "Global_FDR"
    )
  ]

  check_project_unique(
    z,
    paste(
      context,
      feature
    )
  )

  names(z)[-1] <-
    paste0(
      prefix,
      names(z)[-1]
    )

  z
}

# ============================================================
# FILES
# ============================================================

availability_file <- file.path(
  root,
  "JAK3_complete_33cancer_analysis_availability.csv"
)

tn_file <- file.path(
  root,
  "JAK3_complete_tumor_normal.csv"
)

matched_file <- file.path(
  root,
  "JAK3_complete_matched_tumor_normal.csv"
)

stageI_file <- file.path(
  root,
  "JAK3_complete_stageI_vs_normal.csv"
)

prevalence_file <- file.path(
  root,
  "JAK3_complete_aberrant_prevalence_P95.csv"
)

stage_file <- file.path(
  root,
  "JAK3_stage_statistics_complete.csv"
)

grade_file <- file.path(
  root,
  "JAK3_grade_statistics_complete.csv"
)

survival_file <- file.path(
  root,
  "JAK3_TCGA_CDR_survival_complete.csv"
)

lf_file <- file.path(
  root,
  "JAK3_PanImmune_leukocyte_correlations.csv"
)

lf_stage_file <- file.path(
  root,
  "JAK3_PanImmune_stage_adjusted.csv"
)

lf_grade_file <- file.path(
  root,
  "JAK3_PanImmune_grade_adjusted.csv"
)

cibersort_file <- file.path(
  root,
  "JAK3_CIBERSORT_robustness_by_cancer.csv"
)

immune_models_file <- file.path(
  root,
  "immune_markers",
  "JAK3_pan_cancer_immune_adjusted_models.csv"
)

immune_cor_file <- file.path(
  root,
  "immune_markers",
  "JAK3_pan_cancer_immune_marker_correlations.csv"
)

# ============================================================
# READ
# ============================================================

availability <- read_required(
  availability_file
)

tn <- read_required(
  tn_file
)

matched <- read_required(
  matched_file
)

stageI <- read_required(
  stageI_file
)

prevalence <- read_required(
  prevalence_file
)

stage <- read_required(
  stage_file
)

grade <- read_required(
  grade_file
)

survival <- read_required(
  survival_file
)

lf <- read_required(
  lf_file
)

lf_stage <- read_required(
  lf_stage_file
)

lf_grade <- read_required(
  lf_grade_file
)

cibersort <- read_required(
  cibersort_file
)

immune_models <- read_required(
  immune_models_file
)

immune_cor <- read_required(
  immune_cor_file
)

# ============================================================
# BASE AUDIT
# ============================================================

check_project_unique(
  availability,
  "availability"
)

if (nrow(availability) != 33) {
  stop(
    paste(
      "Expected 33 projects in availability table; found",
      nrow(availability)
    )
  )
}

cat("\n========================================\n")
cat("JAK3 MASTER PAN-CANCER INTEGRATION\n")
cat("========================================\n")

cat(
  "\nBase projects:",
  nrow(availability),
  "\n"
)

# ============================================================
# TUMOR VS NORMAL
# ============================================================

tn2 <- tn %>%
  transmute(

    project,

    TN_tumor_n =
      tumor_n,

    TN_normal_n =
      normal_n,

    TN_tumor_median_TPM =
      tumor_median_TPM,

    TN_normal_median_TPM =
      normal_median_TPM,

    TN_median_log2TPM1_difference =
      Median_log2TPM1_difference,

    TN_direction =
      Direction,

    TN_P =
      Wilcoxon_P,

    TN_FDR =
      Wilcoxon_FDR,

    TN_small_normal_n =
      Small_normal_n
  )

check_project_unique(
  tn2,
  "tumor-normal"
)

# ============================================================
# MATCHED TUMOR-NORMAL
# ============================================================

matched2 <- matched %>%
  transmute(

    project,

    Matched_pairs_n =
      matched_pairs,

    Matched_median_log2_difference =
      median_log2_difference,

    Matched_tumor_higher_n =
      tumor_higher_n,

    Matched_tumor_lower_or_equal_n =
      tumor_lower_or_equal_n,

    Matched_P =
      Paired_Wilcoxon_P,

    Matched_FDR =
      Paired_Wilcoxon_FDR,

    Matched_small_n =
      Small_matched_n
  )

check_project_unique(
  matched2,
  "matched tumor-normal"
)

# ============================================================
# STAGE I VS NORMAL
# ============================================================

stageI2 <- stageI %>%
  transmute(

    project,

    StageI_definition =
      Early_stage_definition,

    StageI_tumor_n =
      Stage_I_n,

    StageI_normal_n =
      normal_n,

    StageI_median_TPM =
      Stage_I_median_TPM,

    StageI_normal_median_TPM =
      normal_median_TPM,

    StageI_median_log2TPM1_difference =
      Median_log2TPM1_difference,

    StageI_direction =
      Direction,

    StageI_P =
      Wilcoxon_P,

    StageI_FDR =
      Wilcoxon_FDR,

    StageI_small_n =
      Small_stageI_n,

    StageI_small_normal_n =
      Small_normal_n
  )

check_project_unique(
  stageI2,
  "Stage-I vs normal"
)

# ============================================================
# ABERRANT PREVALENCE
# ============================================================

prev2 <- prevalence %>%
  transmute(

    project,

    Prevalence_tumor_n =
      tumor_n,

    Prevalence_normal_n =
      normal_n,

    Normal_P95_TPM =
      Normal_P95_TPM,

    Aberrant_tumor_n =
      Aberrant_tumor_n,

    Aberrant_tumor_percent =
      Aberrant_tumor_percent,

    Prevalence_small_normal_n =
      Small_normal_n
  )

check_project_unique(
  prev2,
  "prevalence"
)

# ============================================================
# STAGE STATISTICS
# ============================================================

stage2 <- stage %>%
  transmute(

    project,

    Stage_patients_used =
      Patients_used,

    Stage_groups_used =
      Stage_groups_used,

    Stage_groups =
      Groups_used,

    Stage_all_groups =
      All_stage_groups,

    Stage_sparse_groups_excluded =
      Sparse_groups_excluded,

    Stage_KW_P =
      Kruskal_Wallis_P,

    Stage_KW_FDR =
      Kruskal_Wallis_FDR,

    Stage_Spearman_rho =
      Spearman_rho,

    Stage_Spearman_P =
      Spearman_P,

    Stage_trend_FDR =
      Stage_trend_FDR
  )

check_project_unique(
  stage2,
  "stage statistics"
)

# ============================================================
# GRADE STATISTICS
# ============================================================

grade2 <- grade %>%
  transmute(

    project,

    Grade_patients_used =
      Patients_used,

    Grade_groups_used =
      Grade_groups_used,

    Grade_groups =
      Groups_used,

    Grade_all_valid_groups =
      All_valid_grade_groups,

    Grade_sparse_groups_excluded =
      Sparse_valid_groups_excluded,

    Grade_unmapped_labels =
      Unmapped_grade_labels,

    Grade_KW_P =
      Kruskal_Wallis_P,

    Grade_KW_FDR =
      Kruskal_Wallis_FDR,

    Grade_Spearman_rho =
      Spearman_rho,

    Grade_Spearman_P =
      Spearman_P,

    Grade_trend_FDR =
      Grade_trend_FDR
  )

check_project_unique(
  grade2,
  "grade statistics"
)

# ============================================================
# SURVIVAL
# ============================================================

survival2 <- survival %>%
  transmute(

    project,

    Survival_endpoint =
      Endpoint,

    Survival_recommendation =
      Recommendation,

    Survival_N =
      N,

    Survival_events =
      Events,

    Survival_HR =
      HR,

    Survival_CI95_low =
      CI95_low,

    Survival_CI95_high =
      CI95_high,

    Survival_Cox_P =
      Cox_P,

    Survival_Cox_FDR =
      Cox_FDR,

    Survival_PH_P =
      PH_P,

    Survival_PH_violation =
      PH_violation,

    Survival_low_event_n =
      Low_event_n,

    Survival_status =
      Analysis_status,

    Survival_direction =
      Direction
  )

check_project_unique(
  survival2,
  "survival"
)

# ============================================================
# PANIMMUNE LEUKOCYTE FRACTION
# ============================================================

lf2 <- lf %>%
  transmute(

    project,

    LF_N =
      N,

    LF_median =
      Median_leukocyte_fraction,

    LF_Spearman_rho =
      Spearman_rho,

    LF_Spearman_P =
      Spearman_P,

    LF_Spearman_FDR =
      Spearman_FDR,

    LF_linear_beta =
      Linear_beta_LF,

    LF_linear_P =
      Linear_P,

    LF_linear_FDR =
      Linear_FDR,

    LF_R_squared =
      R_squared_LF,

    LF_low_N =
      Low_N
  )

check_project_unique(
  lf2,
  "PanImmune leukocyte correlation"
)

# ============================================================
# PANIMMUNE STAGE-ADJUSTED
# ============================================================

lf_stage2 <- lf_stage %>%
  transmute(

    project,

    LF_stage_N =
      N,

    LF_stage_groups =
      Stage_groups_used,

    LF_stage_levels =
      Stage_levels,

    LF_adjusted_stage_categorical_P =
      Adjusted_stage_categorical_P,

    LF_adjusted_stage_categorical_FDR =
      Adjusted_stage_categorical_FDR,

    LF_adjusted_stage_trend_beta =
      Adjusted_stage_trend_beta,

    LF_adjusted_stage_trend_P =
      Adjusted_stage_trend_P,

    LF_adjusted_stage_trend_FDR =
      Adjusted_stage_trend_FDR,

    LF_stage_leukocyte_beta =
      Leukocyte_beta,

    LF_stage_leukocyte_P =
      Leukocyte_P
  )

check_project_unique(
  lf_stage2,
  "PanImmune stage-adjusted"
)

# ============================================================
# PANIMMUNE GRADE-ADJUSTED
# ============================================================

lf_grade2 <- lf_grade %>%
  transmute(

    project,

    LF_grade_N =
      N,

    LF_grade_groups =
      Grade_groups_used,

    LF_grade_levels =
      Grade_levels,

    LF_adjusted_grade_categorical_P =
      Adjusted_grade_categorical_P,

    LF_adjusted_grade_categorical_FDR =
      Adjusted_grade_categorical_FDR,

    LF_adjusted_grade_trend_beta =
      Adjusted_grade_trend_beta,

    LF_adjusted_grade_trend_P =
      Adjusted_grade_trend_P,

    LF_adjusted_grade_trend_FDR =
      Adjusted_grade_trend_FDR,

    LF_grade_leukocyte_beta =
      Leukocyte_beta,

    LF_grade_leukocyte_P =
      Leukocyte_P
  )

check_project_unique(
  lf_grade2,
  "PanImmune grade-adjusted"
)

# ============================================================
# CIBERSORT ROBUSTNESS
# ============================================================

cib2 <- cibersort %>%
  transmute(

    project,

    CIBERSORT_tested_cells =
      Tested_cells,

    CIBERSORT_primary_significant =
      Primary_significant,

    CIBERSORT_robust_relative_signals =
      Robust_relative_signals,

    CIBERSORT_LF_adjusted_signals =
      LF_adjusted_signals,

    CIBERSORT_strongest_signals =
      Strongest_signals
  )

check_project_unique(
  cib2,
  "CIBERSORT robustness"
)

# ============================================================
# IMMUNE-MARKER CORRELATIONS
#
# Keep PTPRC and aggregate ImmuneScore in master table.
# Full 10-marker results remain in the detailed source CSV.
# ============================================================

cor_primary_ptprc <- extract_correlation(
  immune_cor,
  "Primary tumors",
  "PTPRC",
  "Cor_Primary_PTPRC_"
)

cor_primary_score <- extract_correlation(
  immune_cor,
  "Primary tumors",
  "ImmuneScore",
  "Cor_Primary_ImmuneScore_"
)

cor_stageI_ptprc <- extract_correlation(
  immune_cor,
  "Stage I tumors",
  "PTPRC",
  "Cor_StageI_PTPRC_"
)

cor_stageI_score <- extract_correlation(
  immune_cor,
  "Stage I tumors",
  "ImmuneScore",
  "Cor_StageI_ImmuneScore_"
)

# ============================================================
# IMMUNE-ADJUSTED TUMOR-NORMAL MODELS
# ============================================================

ia_tn_unadj <- extract_model(
  immune_models,
  "All primary tumors vs normal",
  "Unadjusted",
  "IA_TN_Unadjusted_"
)

ia_tn_ptprc <- extract_model(
  immune_models,
  "All primary tumors vs normal",
  "Adjusted for PTPRC/CD45",
  "IA_TN_PTPRC_"
)

ia_tn_score <- extract_model(
  immune_models,
  "All primary tumors vs normal",
  "Adjusted for multi-marker immune score",
  "IA_TN_ImmuneScore_"
)

# ============================================================
# IMMUNE-ADJUSTED STAGE-I-NORMAL MODELS
# ============================================================

ia_stageI_unadj <- extract_model(
  immune_models,
  "Stage I tumors vs normal",
  "Unadjusted",
  "IA_StageI_Unadjusted_"
)

ia_stageI_ptprc <- extract_model(
  immune_models,
  "Stage I tumors vs normal",
  "Adjusted for PTPRC/CD45",
  "IA_StageI_PTPRC_"
)

ia_stageI_score <- extract_model(
  immune_models,
  "Stage I tumors vs normal",
  "Adjusted for multi-marker immune score",
  "IA_StageI_ImmuneScore_"
)

# ============================================================
# BUILD MASTER TABLE
# ============================================================

master <- availability %>%

  left_join(
    tn2,
    by = "project"
  ) %>%

  left_join(
    matched2,
    by = "project"
  ) %>%

  left_join(
    stageI2,
    by = "project"
  ) %>%

  left_join(
    prev2,
    by = "project"
  ) %>%

  left_join(
    stage2,
    by = "project"
  ) %>%

  left_join(
    grade2,
    by = "project"
  ) %>%

  left_join(
    survival2,
    by = "project"
  ) %>%

  left_join(
    lf2,
    by = "project"
  ) %>%

  left_join(
    lf_stage2,
    by = "project"
  ) %>%

  left_join(
    lf_grade2,
    by = "project"
  ) %>%

  left_join(
    cib2,
    by = "project"
  ) %>%

  left_join(
    cor_primary_ptprc,
    by = "project"
  ) %>%

  left_join(
    cor_primary_score,
    by = "project"
  ) %>%

  left_join(
    cor_stageI_ptprc,
    by = "project"
  ) %>%

  left_join(
    cor_stageI_score,
    by = "project"
  ) %>%

  left_join(
    ia_tn_unadj,
    by = "project"
  ) %>%

  left_join(
    ia_tn_ptprc,
    by = "project"
  ) %>%

  left_join(
    ia_tn_score,
    by = "project"
  ) %>%

  left_join(
    ia_stageI_unadj,
    by = "project"
  ) %>%

  left_join(
    ia_stageI_ptprc,
    by = "project"
  ) %>%

  left_join(
    ia_stageI_score,
    by = "project"
  )

# ============================================================
# FINAL MASTER AUDIT
# ============================================================

if (nrow(master) != 33) {
  stop(
    paste(
      "Master table should contain 33 projects; found",
      nrow(master)
    )
  )
}

if (anyDuplicated(master$project)) {
  stop("Duplicated projects detected in final master table.")
}

master <- master %>%
  arrange(project)

# ============================================================
# FACTUAL SUMMARY FLAGS
#
# These are NOT rankings.
# They only make the key statistical results easier to filter.
# ============================================================

flags <- master %>%
  transmute(

    project,

    Tumor_vs_normal_available,

    TN_significant_UP =
      coalesce(
        TN_FDR < 0.05 &
        TN_direction == "UP",
        FALSE
      ),

    TN_significant_DOWN =
      coalesce(
        TN_FDR < 0.05 &
        TN_direction == "DOWN",
        FALSE
      ),

    Matched_available,

    Matched_significant_positive =
      coalesce(
        Matched_FDR < 0.05 &
        Matched_median_log2_difference > 0,
        FALSE
      ),

    Early_stage_vs_normal_available,

    StageI_significant_UP =
      coalesce(
        StageI_FDR < 0.05 &
        StageI_direction == "UP",
        FALSE
      ),

    StageI_significant_DOWN =
      coalesce(
        StageI_FDR < 0.05 &
        StageI_direction == "DOWN",
        FALSE
      ),

    Stage_categorical_significant =
      coalesce(
        Stage_KW_FDR < 0.05,
        FALSE
      ),

    Stage_monotonic_trend_significant =
      coalesce(
        Stage_trend_FDR < 0.05,
        FALSE
      ),

    Grade_categorical_significant =
      coalesce(
        Grade_KW_FDR < 0.05,
        FALSE
      ),

    Grade_monotonic_trend_significant =
      coalesce(
        Grade_trend_FDR < 0.05,
        FALSE
      ),

    Survival_significant =
      coalesce(
        Survival_status == "Analysed" &
        Survival_Cox_FDR < 0.05,
        FALSE
      ),

    JAK3_leukocyte_correlation_significant =
      coalesce(
        LF_Spearman_FDR < 0.05,
        FALSE
      ),

    LF_adjusted_stage_categorical_significant =
      coalesce(
        LF_adjusted_stage_categorical_FDR < 0.05,
        FALSE
      ),

    LF_adjusted_stage_trend_significant =
      coalesce(
        LF_adjusted_stage_trend_FDR < 0.05,
        FALSE
      ),

    LF_adjusted_grade_categorical_significant =
      coalesce(
        LF_adjusted_grade_categorical_FDR < 0.05,
        FALSE
      ),

    LF_adjusted_grade_trend_significant =
      coalesce(
        LF_adjusted_grade_trend_FDR < 0.05,
        FALSE
      ),

    CIBERSORT_has_strongest_signal =
      coalesce(
        CIBERSORT_strongest_signals > 0,
        FALSE
      ),

    TN_PTPRC_adjusted_positive_persists =
      coalesce(
        IA_TN_PTPRC_Model_family_FDR < 0.05 &
        IA_TN_PTPRC_Coefficient > 0 &
        !IA_TN_PTPRC_Low_N,
        FALSE
      ),

    TN_ImmuneScore_adjusted_positive_persists =
      coalesce(
        IA_TN_ImmuneScore_Model_family_FDR < 0.05 &
        IA_TN_ImmuneScore_Coefficient > 0 &
        !IA_TN_ImmuneScore_Low_N,
        FALSE
      ),

    StageI_PTPRC_adjusted_positive_persists =
      coalesce(
        IA_StageI_PTPRC_Model_family_FDR < 0.05 &
        IA_StageI_PTPRC_Coefficient > 0 &
        !IA_StageI_PTPRC_Low_N,
        FALSE
      ),

    StageI_ImmuneScore_adjusted_positive_persists =
      coalesce(
        IA_StageI_ImmuneScore_Model_family_FDR < 0.05 &
        IA_StageI_ImmuneScore_Coefficient > 0 &
        !IA_StageI_ImmuneScore_Low_N,
        FALSE
      ),

    Small_normal_n,
    Small_matched_n,
    Small_stageI_n
  )

# ============================================================
# SAVE
# ============================================================

master_file <- file.path(
  outdir,
  "JAK3_master_pan_cancer_summary.csv"
)

flags_file <- file.path(
  outdir,
  "JAK3_master_key_flags.csv"
)

write.csv(
  master,
  master_file,
  row.names = FALSE
)

write.csv(
  flags,
  flags_file,
  row.names = FALSE
)

# ============================================================
# SOURCE MANIFEST
# ============================================================

source_manifest <- data.frame(

  Module = c(
    "Availability",
    "Tumor vs normal",
    "Matched tumor-normal",
    "Stage-I vs normal",
    "Aberrant prevalence",
    "Stage statistics",
    "Grade statistics",
    "Survival",
    "PanImmune leukocyte fraction",
    "PanImmune stage-adjusted",
    "PanImmune grade-adjusted",
    "CIBERSORT robustness",
    "Immune-marker correlations",
    "Immune-adjusted tumor/Stage-I models"
  ),

  File = c(
    availability_file,
    tn_file,
    matched_file,
    stageI_file,
    prevalence_file,
    stage_file,
    grade_file,
    survival_file,
    lf_file,
    lf_stage_file,
    lf_grade_file,
    cibersort_file,
    immune_cor_file,
    immune_models_file
  ),

  stringsAsFactors = FALSE
)

write.csv(
  source_manifest,
  file.path(
    outdir,
    "JAK3_master_source_manifest.csv"
  ),
  row.names = FALSE
)

# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("54 JAK3 MASTER SUMMARY\n")
cat("========================================\n")

cat(
  "\nProjects in master table:",
  nrow(master),
  "\n"
)

cat(
  "\nSignificant tumor-normal UP:",
  sum(
    flags$TN_significant_UP,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significant matched positive:",
  sum(
    flags$Matched_significant_positive,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significant Stage-I UP:",
  sum(
    flags$StageI_significant_UP,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significant stage monotonic trends:",
  sum(
    flags$Stage_monotonic_trend_significant,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significant grade monotonic trends:",
  sum(
    flags$Grade_monotonic_trend_significant,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significant survival associations:",
  sum(
    flags$Survival_significant,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significant JAK3-leukocyte correlations:",
  sum(
    flags$JAK3_leukocyte_correlation_significant,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Cancers with >=1 strongest CIBERSORT signal:",
  sum(
    flags$CIBERSORT_has_strongest_signal,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nTumor-normal positive association persisting after PTPRC adjustment:",
  sum(
    flags$TN_PTPRC_adjusted_positive_persists,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Tumor-normal positive association persisting after ImmuneScore adjustment:",
  sum(
    flags$TN_ImmuneScore_adjusted_positive_persists,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Stage-I positive association persisting after PTPRC adjustment:",
  sum(
    flags$StageI_PTPRC_adjusted_positive_persists,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Stage-I positive association persisting after ImmuneScore adjustment:",
  sum(
    flags$StageI_ImmuneScore_adjusted_positive_persists,
    na.rm = TRUE
  ),
  "\n"
)

# ============================================================
# KIRC AUDIT
# ============================================================

cat("\n===== KIRC MASTER CHECK =====\n")

kirc <- master[
  master$project == "TCGA-KIRC",
  c(
    "project",
    "TN_median_log2TPM1_difference",
    "TN_FDR",
    "Matched_median_log2_difference",
    "Matched_FDR",
    "StageI_median_log2TPM1_difference",
    "StageI_FDR",
    "Stage_Spearman_rho",
    "Stage_trend_FDR",
    "Grade_Spearman_rho",
    "Grade_trend_FDR",
    "Survival_HR",
    "Survival_Cox_FDR",
    "LF_Spearman_rho",
    "LF_Spearman_FDR",
    "CIBERSORT_strongest_signals",
    "IA_TN_PTPRC_Coefficient",
    "IA_TN_PTPRC_Model_family_FDR",
    "IA_TN_ImmuneScore_Coefficient",
    "IA_TN_ImmuneScore_Model_family_FDR",
    "IA_StageI_PTPRC_Coefficient",
    "IA_StageI_PTPRC_Model_family_FDR",
    "IA_StageI_ImmuneScore_Coefficient",
    "IA_StageI_ImmuneScore_Model_family_FDR"
  )
]

print(
  as.data.frame(kirc),
  row.names = FALSE
)

cat("\nSaved:\n")
cat(" ", master_file, "\n")
cat(" ", flags_file, "\n")
cat(
  " ",
  file.path(
    outdir,
    "JAK3_master_source_manifest.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("54 COMPLETE\n")
cat("========================================\n")
