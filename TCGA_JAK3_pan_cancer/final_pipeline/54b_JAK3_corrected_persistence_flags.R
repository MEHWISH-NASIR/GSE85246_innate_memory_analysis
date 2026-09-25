library(dplyr)

# ============================================================
# 54b — CORRECT JAK3 IMMUNE-ADJUSTMENT PERSISTENCE FLAGS
#
# "Persists" requires:
#   1. unadjusted OLS association significant + positive
#   2. corresponding adjusted OLS association significant + positive
#   3. no Low_N flag
#
# This avoids comparing canonical Wilcoxon significance directly
# with adjusted linear-model significance.
# ============================================================

master_file <-
  "TCGA_JAK3_pan_cancer/results/master/JAK3_master_pan_cancer_summary.csv"

outfile <-
  "TCGA_JAK3_pan_cancer/results/master/JAK3_master_key_flags_corrected.csv"

stopifnot(file.exists(master_file))

x <- read.csv(
  master_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required <- c(
  "project",

  "TN_FDR",
  "TN_direction",
  "StageI_FDR",
  "StageI_direction",

  "IA_TN_Unadjusted_Coefficient",
  "IA_TN_Unadjusted_Model_family_FDR",
  "IA_TN_Unadjusted_Low_N",

  "IA_TN_PTPRC_Coefficient",
  "IA_TN_PTPRC_Model_family_FDR",
  "IA_TN_PTPRC_Low_N",

  "IA_TN_ImmuneScore_Coefficient",
  "IA_TN_ImmuneScore_Model_family_FDR",
  "IA_TN_ImmuneScore_Low_N",

  "IA_StageI_Unadjusted_Coefficient",
  "IA_StageI_Unadjusted_Model_family_FDR",
  "IA_StageI_Unadjusted_Low_N",

  "IA_StageI_PTPRC_Coefficient",
  "IA_StageI_PTPRC_Model_family_FDR",
  "IA_StageI_PTPRC_Low_N",

  "IA_StageI_ImmuneScore_Coefficient",
  "IA_StageI_ImmuneScore_Model_family_FDR",
  "IA_StageI_ImmuneScore_Low_N"
)

missing_cols <- setdiff(required, names(x))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

sig_positive <- function(coef, fdr, low_n) {

  !is.na(coef) &
  !is.na(fdr) &
  coef > 0 &
  fdr < 0.05 &
  (is.na(low_n) | !low_n)
}

flags <- x %>%
  transmute(

    project,

    # Canonical non-parametric results
    Canonical_TN_significant_UP =
      coalesce(
        TN_FDR < 0.05 &
        TN_direction == "UP",
        FALSE
      ),

    Canonical_StageI_significant_UP =
      coalesce(
        StageI_FDR < 0.05 &
        StageI_direction == "UP",
        FALSE
      ),

    # Same-model-family OLS results
    TN_OLS_unadjusted_positive =
      sig_positive(
        IA_TN_Unadjusted_Coefficient,
        IA_TN_Unadjusted_Model_family_FDR,
        IA_TN_Unadjusted_Low_N
      ),

    TN_PTPRC_adjusted_positive =
      sig_positive(
        IA_TN_PTPRC_Coefficient,
        IA_TN_PTPRC_Model_family_FDR,
        IA_TN_PTPRC_Low_N
      ),

    TN_ImmuneScore_adjusted_positive =
      sig_positive(
        IA_TN_ImmuneScore_Coefficient,
        IA_TN_ImmuneScore_Model_family_FDR,
        IA_TN_ImmuneScore_Low_N
      ),

    StageI_OLS_unadjusted_positive =
      sig_positive(
        IA_StageI_Unadjusted_Coefficient,
        IA_StageI_Unadjusted_Model_family_FDR,
        IA_StageI_Unadjusted_Low_N
      ),

    StageI_PTPRC_adjusted_positive =
      sig_positive(
        IA_StageI_PTPRC_Coefficient,
        IA_StageI_PTPRC_Model_family_FDR,
        IA_StageI_PTPRC_Low_N
      ),

    StageI_ImmuneScore_adjusted_positive =
      sig_positive(
        IA_StageI_ImmuneScore_Coefficient,
        IA_StageI_ImmuneScore_Model_family_FDR,
        IA_StageI_ImmuneScore_Low_N
      )
  )

flags <- flags %>%
  mutate(

    TN_positive_persists_after_PTPRC =
      TN_OLS_unadjusted_positive &
      TN_PTPRC_adjusted_positive,

    TN_positive_persists_after_ImmuneScore =
      TN_OLS_unadjusted_positive &
      TN_ImmuneScore_adjusted_positive,

    StageI_positive_persists_after_PTPRC =
      StageI_OLS_unadjusted_positive &
      StageI_PTPRC_adjusted_positive,

    StageI_positive_persists_after_ImmuneScore =
      StageI_OLS_unadjusted_positive &
      StageI_ImmuneScore_adjusted_positive,

    # Cross-method support kept separate and explicitly named
    Canonical_TN_UP_and_ImmuneScore_adjusted_positive =
      Canonical_TN_significant_UP &
      TN_ImmuneScore_adjusted_positive,

    Canonical_StageI_UP_and_ImmuneScore_adjusted_positive =
      Canonical_StageI_significant_UP &
      StageI_ImmuneScore_adjusted_positive
  )

write.csv(
  flags,
  outfile,
  row.names = FALSE
)

cat("\n========================================\n")
cat("54b CORRECTED PERSISTENCE AUDIT\n")
cat("========================================\n")

cat(
  "\nCanonical Wilcoxon tumor-normal UP:",
  sum(flags$Canonical_TN_significant_UP),
  "\n"
)

cat(
  "Unadjusted OLS tumor-normal positive:",
  sum(flags$TN_OLS_unadjusted_positive),
  "\n"
)

cat(
  "Persists after PTPRC:",
  sum(flags$TN_positive_persists_after_PTPRC),
  "\n"
)

cat(
  "Persists after ImmuneScore:",
  sum(flags$TN_positive_persists_after_ImmuneScore),
  "\n"
)

cat(
  "\nCanonical Wilcoxon Stage-I UP:",
  sum(flags$Canonical_StageI_significant_UP),
  "\n"
)

cat(
  "Unadjusted OLS Stage-I positive:",
  sum(flags$StageI_OLS_unadjusted_positive),
  "\n"
)

cat(
  "Persists after PTPRC:",
  sum(flags$StageI_positive_persists_after_PTPRC),
  "\n"
)

cat(
  "Persists after ImmuneScore:",
  sum(flags$StageI_positive_persists_after_ImmuneScore),
  "\n"
)

cat("\n===== KIRC CHECK =====\n")

print(
  flags[
    flags$project == "TCGA-KIRC",
  ],
  row.names = FALSE
)

cat("\nSaved:\n ", outfile, "\n")

cat("\n========================================\n")
cat("54b COMPLETE\n")
cat("========================================\n")
