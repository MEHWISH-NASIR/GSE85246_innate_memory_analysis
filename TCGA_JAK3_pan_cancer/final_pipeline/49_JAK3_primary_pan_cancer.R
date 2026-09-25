library(dplyr)

infile <- paste0(
  "TCGA_JAK3_pan_cancer/results/expression/",
  "JAK3_expression_clinical_complete_33projects.csv"
)

d <- read.csv(
  infile,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & trimws(as.character(x)) != ""]
  if (length(x) == 0) NA_character_ else as.character(x[1])
}

safe_wilcox_unpaired <- function(x, y) {
  tryCatch(
    wilcox.test(
      x, y,
      paired = FALSE,
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

safe_wilcox_paired <- function(x, y) {
  tryCatch(
    wilcox.test(
      x, y,
      paired = TRUE,
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

# ------------------------------------------------------------
# Canonical sample classes
# ------------------------------------------------------------

d <- d %>%
  mutate(
    analysis_class = case_when(
      project == "TCGA-LAML" &
        sample_type ==
        "Primary Blood Derived Cancer - Peripheral Blood" ~
        "Primary malignant",

      project != "TCGA-LAML" &
        sample_type == "Primary Tumor" ~
        "Primary malignant",

      sample_type == "Solid Tissue Normal" ~
        "Normal",

      TRUE ~ "Other"
    )
  )

# ------------------------------------------------------------
# Collapse to one row per patient
# ------------------------------------------------------------

tumor <- d %>%
  filter(analysis_class == "Primary malignant") %>%
  group_by(project, cases.submitter_id) %>%
  summarise(
    TPM = mean(TPM, na.rm = TRUE),

    stage = first_nonmissing(
      ajcc_pathologic_stage
    ),

    expression_samples = n(),
    .groups = "drop"
  ) %>%
  mutate(
    log2TPM1 = log2(TPM + 1),

    stage_clean = toupper(
      trimws(stage)
    ),

    is_stage_I =
      !is.na(stage_clean) &
      grepl(
        "^STAGE I([ABC])?$",
        stage_clean
      )
  )

normal <- d %>%
  filter(analysis_class == "Normal") %>%
  group_by(project, cases.submitter_id) %>%
  summarise(
    TPM = mean(TPM, na.rm = TRUE),
    expression_samples = n(),
    .groups = "drop"
  ) %>%
  mutate(
    log2TPM1 = log2(TPM + 1)
  )

# save patient-level datasets
write.csv(
  tumor,
  "TCGA_JAK3_pan_cancer/results/JAK3_primary_tumor_patient_level.csv",
  row.names = FALSE
)

write.csv(
  normal,
  "TCGA_JAK3_pan_cancer/results/JAK3_normal_patient_level.csv",
  row.names = FALSE
)

# ============================================================
# 1. TUMOR VS NORMAL
# ============================================================

projects_tn <- intersect(
  unique(tumor$project),
  unique(normal$project)
)

tn_results <- lapply(
  projects_tn,
  function(p) {

    t <- tumor %>%
      filter(project == p)

    n <- normal %>%
      filter(project == p)

    data.frame(
      project = p,

      tumor_n = nrow(t),
      normal_n = nrow(n),

      tumor_median_TPM =
        median(t$TPM, na.rm = TRUE),

      normal_median_TPM =
        median(n$TPM, na.rm = TRUE),

      Median_log2TPM1_difference =
        median(t$log2TPM1, na.rm = TRUE) -
        median(n$log2TPM1, na.rm = TRUE),

      Direction =
        ifelse(
          median(t$log2TPM1, na.rm = TRUE) >
            median(n$log2TPM1, na.rm = TRUE),
          "UP",
          "DOWN"
        ),

      Wilcoxon_P =
        safe_wilcox_unpaired(
          t$log2TPM1,
          n$log2TPM1
        ),

      Small_normal_n =
        nrow(n) < 10,

      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

tn_results$Wilcoxon_FDR <- p.adjust(
  tn_results$Wilcoxon_P,
  method = "BH"
)

tn_results <- tn_results %>%
  arrange(Wilcoxon_FDR)

# ============================================================
# 2. MATCHED TUMOR-NORMAL
# ============================================================

matched <- inner_join(
  tumor %>%
    select(
      project,
      cases.submitter_id,
      tumor_TPM = TPM,
      tumor_log2TPM1 = log2TPM1
    ),

  normal %>%
    select(
      project,
      cases.submitter_id,
      normal_TPM = TPM,
      normal_log2TPM1 = log2TPM1
    ),

  by = c(
    "project",
    "cases.submitter_id"
  )
) %>%
  mutate(
    log2_difference =
      tumor_log2TPM1 -
      normal_log2TPM1,

    tumor_higher =
      tumor_log2TPM1 >
      normal_log2TPM1
  )

write.csv(
  matched,
  "TCGA_JAK3_pan_cancer/results/JAK3_matched_pairs_complete.csv",
  row.names = FALSE
)

matched_results <- matched %>%
  group_by(project) %>%
  summarise(
    matched_pairs = n(),

    median_log2_difference =
      median(
        log2_difference,
        na.rm = TRUE
      ),

    tumor_higher_n =
      sum(tumor_higher, na.rm = TRUE),

    tumor_lower_or_equal_n =
      sum(!tumor_higher, na.rm = TRUE),

    Paired_Wilcoxon_P =
      safe_wilcox_paired(
        tumor_log2TPM1,
        normal_log2TPM1
      ),

    Small_matched_n =
      n() < 10,

    .groups = "drop"
  )

matched_results$Paired_Wilcoxon_FDR <-
  p.adjust(
    matched_results$Paired_Wilcoxon_P,
    method = "BH"
  )

matched_results <- matched_results %>%
  arrange(Paired_Wilcoxon_FDR)

# ============================================================
# 3. EXACT STAGE-I VS NORMAL
# ============================================================

early_results <- lapply(
  projects_tn,
  function(p) {

    e <- tumor %>%
      filter(
        project == p,
        is_stage_I
      )

    n <- normal %>%
      filter(project == p)

    if (nrow(e) == 0) {
      return(NULL)
    }

    data.frame(
      project = p,

      Early_stage_definition =
        "AJCC Stage I/IA/IB/IC only",

      Stage_I_n = nrow(e),
      normal_n = nrow(n),

      Stage_I_median_TPM =
        median(e$TPM, na.rm = TRUE),

      normal_median_TPM =
        median(n$TPM, na.rm = TRUE),

      Median_log2TPM1_difference =
        median(
          e$log2TPM1,
          na.rm = TRUE
        ) -
        median(
          n$log2TPM1,
          na.rm = TRUE
        ),

      Direction =
        ifelse(
          median(e$log2TPM1, na.rm = TRUE) >
            median(n$log2TPM1, na.rm = TRUE),
          "UP",
          "DOWN"
        ),

      Wilcoxon_P =
        safe_wilcox_unpaired(
          e$log2TPM1,
          n$log2TPM1
        ),

      Small_stageI_n =
        nrow(e) < 10,

      Small_normal_n =
        nrow(n) < 10,

      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

early_results$Wilcoxon_FDR <-
  p.adjust(
    early_results$Wilcoxon_P,
    method = "BH"
  )

early_results <- early_results %>%
  arrange(Wilcoxon_FDR)

# ============================================================
# 4. NORMAL-DERIVED P95 PREVALENCE
# ============================================================

prevalence_results <- lapply(
  projects_tn,
  function(p) {

    t <- tumor %>%
      filter(project == p)

    n <- normal %>%
      filter(project == p)

    threshold <- as.numeric(
      quantile(
        n$TPM,
        probs = 0.95,
        na.rm = TRUE,
        names = FALSE
      )
    )

    aberrant_n <-
      sum(
        t$TPM > threshold,
        na.rm = TRUE
      )

    data.frame(
      project = p,

      tumor_n = nrow(t),
      normal_n = nrow(n),

      Normal_P95_TPM = threshold,

      Aberrant_tumor_n =
        aberrant_n,

      Aberrant_tumor_percent =
        100 *
        aberrant_n /
        nrow(t),

      Small_normal_n =
        nrow(n) < 10,

      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

# ============================================================
# SAVE
# ============================================================

write.csv(
  tn_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_complete_tumor_normal.csv",
  row.names = FALSE
)

write.csv(
  matched_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_complete_matched_tumor_normal.csv",
  row.names = FALSE
)

write.csv(
  early_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_complete_stageI_vs_normal.csv",
  row.names = FALSE
)

write.csv(
  prevalence_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_complete_aberrant_prevalence_P95.csv",
  row.names = FALSE
)

# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("JAK3 PRIMARY PAN-CANCER ANALYSIS\n")
cat("========================================\n")

cat("\nTumor-normal cancers:",
    nrow(tn_results), "\n")

cat("Matched cancers:",
    nrow(matched_results), "\n")

cat("Stage-I vs normal cancers:",
    nrow(early_results), "\n")

cat("P95 prevalence cancers:",
    nrow(prevalence_results), "\n")

cat(
  "\nTumor-normal FDR < 0.05:",
  sum(
    tn_results$Wilcoxon_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Matched FDR < 0.05:",
  sum(
    matched_results$Paired_Wilcoxon_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Stage-I FDR < 0.05:",
  sum(
    early_results$Wilcoxon_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat("\n===== TOP TUMOR-NORMAL =====\n")
print(
  head(tn_results, 10),
  row.names = FALSE
)

cat("\n===== TOP MATCHED =====\n")
print(
  head(matched_results, 10),
  row.names = FALSE
)

cat("\n===== TOP STAGE-I =====\n")
print(
  head(early_results, 10),
  row.names = FALSE
)

cat("\n49 PRIMARY PAN-CANCER ANALYSIS COMPLETE\n")
