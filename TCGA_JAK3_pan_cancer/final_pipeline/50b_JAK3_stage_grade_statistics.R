library(dplyr)

d <- read.csv(
  "TCGA_JAK3_pan_cancer/results/JAK3_stage_grade_patient_level.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

d <- d %>%
  mutate(
    log2TPM1 = log2(TPM + 1),

    stage_num = case_when(
      stage_simple == "Stage I"   ~ 1,
      stage_simple == "Stage II"  ~ 2,
      stage_simple == "Stage III" ~ 3,
      stage_simple == "Stage IV"  ~ 4,
      TRUE ~ NA_real_
    ),

    grade_num = case_when(

      # BLCA uses its own binary grade convention
      project == "TCGA-BLCA" &
        grade_raw == "Low Grade"  ~ 1,

      project == "TCGA-BLCA" &
        grade_raw == "High Grade" ~ 2,

      # Standard G1-G4 scale for other projects
      project != "TCGA-BLCA" &
        grade_raw == "G1" ~ 1,

      project != "TCGA-BLCA" &
        grade_raw == "G2" ~ 2,

      project != "TCGA-BLCA" &
        grade_raw == "G3" ~ 3,

      project != "TCGA-BLCA" &
        grade_raw == "G4" ~ 4,

      TRUE ~ NA_real_
    ),

    grade_harmonized = case_when(
      project == "TCGA-BLCA" & grade_num == 1 ~ "Low Grade",
      project == "TCGA-BLCA" & grade_num == 2 ~ "High Grade",

      project != "TCGA-BLCA" & grade_num == 1 ~ "G1",
      project != "TCGA-BLCA" & grade_num == 2 ~ "G2",
      project != "TCGA-BLCA" & grade_num == 3 ~ "G3",
      project != "TCGA-BLCA" & grade_num == 4 ~ "G4",

      TRUE ~ NA_character_
    )
  )

safe_kw <- function(x, g) {
  tryCatch(
    kruskal.test(x ~ factor(g))$p.value,
    error = function(e) NA_real_
  )
}

safe_spearman <- function(x, y) {

  z <- complete.cases(x, y)

  if (sum(z) < 5 ||
      length(unique(x[z])) < 2) {
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
# STAGE DESCRIPTIVE
# ============================================================

stage_desc <- d %>%
  filter(!is.na(stage_num)) %>%
  group_by(
    project,
    stage_simple,
    stage_num
  ) %>%
  summarise(
    n = n(),

    median_TPM =
      median(TPM, na.rm = TRUE),

    median_log2TPM1 =
      median(log2TPM1, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  arrange(project, stage_num)

write.csv(
  stage_desc,
  "TCGA_JAK3_pan_cancer/results/JAK3_stage_descriptive_complete.csv",
  row.names = FALSE
)

# ============================================================
# STAGE INFERENTIAL
# Require >=5 patients within a stage group
# ============================================================

stage_eligible <- d %>%
  filter(!is.na(stage_num)) %>%
  group_by(
    project,
    stage_simple,
    stage_num
  ) %>%
  mutate(
    stage_group_n = n()
  ) %>%
  ungroup() %>%
  filter(stage_group_n >= 5)

stage_projects <- stage_eligible %>%
  group_by(project) %>%
  summarise(
    groups_used = n_distinct(stage_num),
    patients_used = n(),
    .groups = "drop"
  ) %>%
  filter(groups_used >= 2)

stage_results <- lapply(
  stage_projects$project,
  function(p) {

    x <- stage_eligible %>%
      filter(project == p)

    sp <- safe_spearman(
      x$stage_num,
      x$log2TPM1
    )

    original <- stage_desc %>%
      filter(project == p)

    data.frame(
      project = p,

      Patients_used = nrow(x),

      Stage_groups_used =
        n_distinct(x$stage_num),

      Groups_used =
        paste(
          sort(unique(x$stage_simple)),
          collapse = "; "
        ),

      All_stage_groups =
        paste(
          original$stage_simple,
          original$n,
          sep = ":n=",
          collapse = "; "
        ),

      Sparse_groups_excluded =
        any(original$n < 5),

      Kruskal_Wallis_P =
        safe_kw(
          x$log2TPM1,
          x$stage_num
        ),

      Spearman_rho =
        sp["rho"],

      Spearman_P =
        sp["p"],

      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

stage_results$Kruskal_Wallis_FDR <-
  p.adjust(
    stage_results$Kruskal_Wallis_P,
    method = "BH"
  )

stage_results$Stage_trend_FDR <-
  p.adjust(
    stage_results$Spearman_P,
    method = "BH"
  )

stage_results <- stage_results %>%
  arrange(Stage_trend_FDR)

write.csv(
  stage_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_stage_statistics_complete.csv",
  row.names = FALSE
)

# ============================================================
# GRADE LABEL AUDIT / EXCLUDED LABELS
# ============================================================

grade_excluded <- d %>%
  filter(
    !is.na(grade_raw),
    is.na(grade_num)
  ) %>%
  count(
    project,
    grade_raw,
    name = "patients"
  ) %>%
  arrange(project, grade_raw)

write.csv(
  grade_excluded,
  "TCGA_JAK3_pan_cancer/results/JAK3_grade_labels_excluded.csv",
  row.names = FALSE
)

# ============================================================
# GRADE DESCRIPTIVE
# ============================================================

grade_desc <- d %>%
  filter(!is.na(grade_num)) %>%
  group_by(
    project,
    grade_harmonized,
    grade_num
  ) %>%
  summarise(
    n = n(),

    median_TPM =
      median(TPM, na.rm = TRUE),

    median_log2TPM1 =
      median(log2TPM1, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  arrange(project, grade_num)

write.csv(
  grade_desc,
  "TCGA_JAK3_pan_cancer/results/JAK3_grade_descriptive_complete.csv",
  row.names = FALSE
)

# ============================================================
# GRADE INFERENTIAL
# Require >=5 patients within a grade group
# ============================================================

grade_eligible <- d %>%
  filter(!is.na(grade_num)) %>%
  group_by(
    project,
    grade_harmonized,
    grade_num
  ) %>%
  mutate(
    grade_group_n = n()
  ) %>%
  ungroup() %>%
  filter(grade_group_n >= 5)

grade_projects <- grade_eligible %>%
  group_by(project) %>%
  summarise(
    groups_used = n_distinct(grade_num),
    patients_used = n(),
    .groups = "drop"
  ) %>%
  filter(groups_used >= 2)

grade_results <- lapply(
  grade_projects$project,
  function(p) {

    x <- grade_eligible %>%
      filter(project == p)

    sp <- safe_spearman(
      x$grade_num,
      x$log2TPM1
    )

    original <- grade_desc %>%
      filter(project == p)

    excluded_labels <- grade_excluded %>%
      filter(project == p)

    data.frame(
      project = p,

      Patients_used = nrow(x),

      Grade_groups_used =
        n_distinct(x$grade_num),

      Groups_used =
        paste(
          sort(unique(x$grade_harmonized)),
          collapse = "; "
        ),

      All_valid_grade_groups =
        paste(
          original$grade_harmonized,
          original$n,
          sep = ":n=",
          collapse = "; "
        ),

      Sparse_valid_groups_excluded =
        any(original$n < 5),

      Unmapped_grade_labels =
        ifelse(
          nrow(excluded_labels) == 0,
          "",
          paste(
            excluded_labels$grade_raw,
            excluded_labels$patients,
            sep = ":n=",
            collapse = "; "
          )
        ),

      Kruskal_Wallis_P =
        safe_kw(
          x$log2TPM1,
          x$grade_num
        ),

      Spearman_rho =
        sp["rho"],

      Spearman_P =
        sp["p"],

      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

grade_results$Kruskal_Wallis_FDR <-
  p.adjust(
    grade_results$Kruskal_Wallis_P,
    method = "BH"
  )

grade_results$Grade_trend_FDR <-
  p.adjust(
    grade_results$Spearman_P,
    method = "BH"
  )

grade_results <- grade_results %>%
  arrange(Grade_trend_FDR)

write.csv(
  grade_results,
  "TCGA_JAK3_pan_cancer/results/JAK3_grade_statistics_complete.csv",
  row.names = FALSE
)

# ============================================================
# CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("JAK3 STAGE / GRADE STATISTICS\n")
cat("========================================\n")

cat(
  "\nStage cancers with robust inferential data:",
  nrow(stage_results),
  "\n"
)

cat(
  "Stage KW FDR < 0.05:",
  sum(
    stage_results$Kruskal_Wallis_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Stage trend FDR < 0.05:",
  sum(
    stage_results$Stage_trend_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nGrade cancers with robust inferential data:",
  nrow(grade_results),
  "\n"
)

cat(
  "Grade KW FDR < 0.05:",
  sum(
    grade_results$Kruskal_Wallis_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Grade trend FDR < 0.05:",
  sum(
    grade_results$Grade_trend_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat("\n===== STAGE RESULTS =====\n")
print(
  as.data.frame(stage_results),
  row.names = FALSE
)

cat("\n===== GRADE RESULTS =====\n")
print(
  as.data.frame(grade_results),
  row.names = FALSE
)

cat("\n===== EXCLUDED GRADE LABELS =====\n")
print(
  as.data.frame(grade_excluded),
  row.names = FALSE
)

cat("\n50b STAGE/GRADE STATISTICS COMPLETE\n")
