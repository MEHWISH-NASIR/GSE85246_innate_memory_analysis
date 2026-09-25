library(dplyr)

primary <- read.csv(
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_relative_correlations.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sens <- read.csv(
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_all_samples_sensitivity.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

adj <- read.csv(
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_LF_adjusted_cell_models.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Primary vs all-sample sensitivity concordance
# ------------------------------------------------------------

robust <- primary %>%
  rename(
    Primary_N = N,
    Primary_rho = rho,
    Primary_P = P,
    Primary_Global_FDR = Global_FDR,
    Primary_Within_cancer_FDR = Within_cancer_FDR,
    Primary_Low_N = Low_N
  ) %>%

  left_join(
    sens %>%
      rename(
        Sensitivity_N = N,
        Sensitivity_rho = rho,
        Sensitivity_P = P,
        Sensitivity_Global_FDR = Global_FDR
      ),
    by = c(
      "project",
      "Cell_type"
    )
  ) %>%

  mutate(
    Same_direction =
      case_when(
        is.na(Primary_rho) |
          is.na(Sensitivity_rho) ~ NA,

        sign(Primary_rho) ==
          sign(Sensitivity_rho) ~ TRUE,

        TRUE ~ FALSE
      ),

    Primary_significant =
      !Primary_Low_N &
      !is.na(Primary_Global_FDR) &
      Primary_Global_FDR < 0.05,

    Sensitivity_significant =
      !is.na(Sensitivity_Global_FDR) &
      Sensitivity_Global_FDR < 0.05,

    Robust_relative_signal =
      Primary_significant &
      Sensitivity_significant &
      Same_direction %in% TRUE
  )

# ------------------------------------------------------------
# Add leukocyte-fraction-adjusted cell result
# ------------------------------------------------------------

robust <- robust %>%
  left_join(
    adj %>%
      select(
        project,
        Cell_type,
        Adjusted_N = N,
        Cell_beta,
        Cell_P,
        Global_Cell_FDR,
        Within_cancer_Cell_FDR,
        Adjusted_Low_N = Low_N
      ),
    by = c(
      "project",
      "Cell_type"
    )
  ) %>%

  mutate(
    LF_adjusted_significant =
      !Adjusted_Low_N &
      !is.na(Global_Cell_FDR) &
      Global_Cell_FDR < 0.05,

    Relative_and_adjusted_consistent =
      case_when(
        is.na(Primary_rho) |
          is.na(Cell_beta) ~ NA,

        sign(Primary_rho) ==
          sign(Cell_beta) ~ TRUE,

        TRUE ~ FALSE
      ),

    Strongest_CIBERSORT_evidence =
      Robust_relative_signal &
      LF_adjusted_significant &
      Relative_and_adjusted_consistent %in% TRUE
  )

# ------------------------------------------------------------
# Summary by cancer
# ------------------------------------------------------------

cancer_summary <- robust %>%
  group_by(project) %>%
  summarise(
    Tested_cells = n(),

    Primary_significant =
      sum(
        Primary_significant,
        na.rm = TRUE
      ),

    Robust_relative_signals =
      sum(
        Robust_relative_signal,
        na.rm = TRUE
      ),

    LF_adjusted_signals =
      sum(
        LF_adjusted_significant,
        na.rm = TRUE
      ),

    Strongest_signals =
      sum(
        Strongest_CIBERSORT_evidence,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  arrange(
    desc(Strongest_signals),
    desc(Robust_relative_signals)
  )

# ------------------------------------------------------------
# Strongest individual associations
# ------------------------------------------------------------

strong <- robust %>%
  filter(
    Strongest_CIBERSORT_evidence
  ) %>%
  arrange(
    Global_Cell_FDR,
    Primary_Global_FDR
  )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  robust,
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_robustness_complete.csv",
  row.names = FALSE
)

write.csv(
  cancer_summary,
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_robustness_by_cancer.csv",
  row.names = FALSE
)

write.csv(
  strong,
  "TCGA_JAK3_pan_cancer/results/JAK3_CIBERSORT_strongest_signals.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

cat("\n========================================\n")
cat("JAK3 CIBERSORT ROBUSTNESS SUMMARY\n")
cat("========================================\n")

cat(
  "\nPrimary globally significant:",
  sum(
    robust$Primary_significant,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Primary + sensitivity significant, same direction:",
  sum(
    robust$Robust_relative_signal,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "LF-adjusted globally significant:",
  sum(
    robust$LF_adjusted_significant,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Strongest evidence across all three checks:",
  sum(
    robust$Strongest_CIBERSORT_evidence,
    na.rm = TRUE
  ),
  "\n"
)

cat("\n===== BY CANCER =====\n")
print(
  as.data.frame(cancer_summary),
  row.names = FALSE
)

cat("\n===== STRONGEST CELL-TYPE SIGNALS =====\n")

print(
  head(
    as.data.frame(strong),
    40
  ),
  row.names = FALSE
)

cat("\n52e CIBERSORT ROBUSTNESS SUMMARY COMPLETE\n")
