
library(data.table)
setDTthreads(1)

root <- "TCGA_multikinase_analysis"

outdir <- file.path(root, "results/stage_grade")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("1. Loading multikinase clinical dataset...\n")

d <- fread(
  file.path(
    root,
    "results/expression",
    "03_multikinase_expression_clinical_complete.csv"
  ),
  nThread = 1
)

# ============================================================
# Helpers
# ============================================================

first_nonmissing <- function(x) {

  x <- as.character(x)
  x <- x[
    !is.na(x) &
    trimws(x) != ""
  ]

  if (length(x) == 0) {
    NA_character_
  } else {
    x[1]
  }
}

clean_missing <- function(x) {

  x <- trimws(as.character(x))

  bad <- tolower(x) %in% c(
    "",
    "na",
    "n/a",
    "unknown",
    "not reported",
    "not available",
    "[not available]",
    "--"
  )

  x[bad] <- NA_character_

  x
}

safe_kw <- function(x, g) {

  tryCatch(
    kruskal.test(
      x ~ factor(g)
    )$p.value,
    error = function(e) NA_real_
  )
}

safe_spearman <- function(x, y) {

  z <- complete.cases(x, y)

  if (
    sum(z) < 5 ||
    length(unique(x[z])) < 2
  ) {
    return(c(
      rho = NA_real_,
      p = NA_real_
    ))
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
    return(c(
      rho = NA_real_,
      p = NA_real_
    ))
  }

  c(
    rho = unname(out$estimate),
    p = out$p.value
  )
}

# ============================================================
# Canonical primary malignant samples
# ============================================================

d[, primary_malignant :=
  (
    project == "TCGA-LAML" &
    sample_type ==
    "Primary Blood Derived Cancer - Peripheral Blood"
  ) |
  (
    project != "TCGA-LAML" &
    sample_type == "Primary Tumor"
  )
]

# ============================================================
# Patient-level collapse
# ============================================================

patient <- d[
  primary_malignant == TRUE,
  .(
    TPM = mean(TPM, na.rm = TRUE),

    stage_raw =
      first_nonmissing(
        ajcc_pathologic_stage
      ),

    grade_raw =
      first_nonmissing(
        tumor_grade
      ),

    expression_samples = .N
  ),
  by = .(
    gene,
    project,
    patient_id
  )
]

patient[, stage_raw :=
  clean_missing(stage_raw)
]

patient[, grade_raw :=
  clean_missing(grade_raw)
]

patient[, stage_upper :=
  toupper(stage_raw)
]

# Order matters: IV before III before II before I
patient[, stage_simple :=
  fifelse(
    grepl("^STAGE IV", stage_upper),
    "Stage IV",

    fifelse(
      grepl("^STAGE III", stage_upper),
      "Stage III",

      fifelse(
        grepl("^STAGE II", stage_upper),
        "Stage II",

        fifelse(
          grepl("^STAGE I", stage_upper),
          "Stage I",
          NA_character_
        )
      )
    )
  )
]

patient[, log2TPM1 :=
  log2(TPM + 1)
]

patient[, stage_num :=
  fifelse(
    stage_simple == "Stage I", 1,
    fifelse(
      stage_simple == "Stage II", 2,
      fifelse(
        stage_simple == "Stage III", 3,
        fifelse(
          stage_simple == "Stage IV", 4,
          NA_real_
        )
      )
    )
  )
]

# ============================================================
# Grade harmonization
# ============================================================

patient[, grade_num := NA_real_]

patient[
  project == "TCGA-BLCA" &
  grade_raw == "Low Grade",
  grade_num := 1
]

patient[
  project == "TCGA-BLCA" &
  grade_raw == "High Grade",
  grade_num := 2
]

patient[
  project != "TCGA-BLCA" &
  grade_raw == "G1",
  grade_num := 1
]

patient[
  project != "TCGA-BLCA" &
  grade_raw == "G2",
  grade_num := 2
]

patient[
  project != "TCGA-BLCA" &
  grade_raw == "G3",
  grade_num := 3
]

patient[
  project != "TCGA-BLCA" &
  grade_raw == "G4",
  grade_num := 4
]

patient[, grade_harmonized :=
  fifelse(
    project == "TCGA-BLCA" &
    grade_num == 1,
    "Low Grade",

    fifelse(
      project == "TCGA-BLCA" &
      grade_num == 2,
      "High Grade",

      fifelse(
        project != "TCGA-BLCA" &
        grade_num == 1,
        "G1",

        fifelse(
          project != "TCGA-BLCA" &
          grade_num == 2,
          "G2",

          fifelse(
            project != "TCGA-BLCA" &
            grade_num == 3,
            "G3",

            fifelse(
              project != "TCGA-BLCA" &
              grade_num == 4,
              "G4",
              NA_character_
            )
          )
        )
      )
    )
  )
]

fwrite(
  patient,
  file.path(
    outdir,
    "06_multikinase_stage_grade_patient_level.csv"
  )
)

cat("2. Patient-level stage/grade dataset prepared\n")

# ============================================================
# STAGE DESCRIPTIVE
# ============================================================

stage_desc <- patient[
  !is.na(stage_num),
  .(
    n = .N,
    median_TPM =
      median(TPM, na.rm = TRUE),
    median_log2TPM1 =
      median(log2TPM1, na.rm = TRUE)
  ),
  by = .(
    gene,
    project,
    stage_simple,
    stage_num
  )
]

setorder(
  stage_desc,
  gene,
  project,
  stage_num
)

fwrite(
  stage_desc,
  file.path(
    outdir,
    "06_stage_descriptive.csv"
  )
)

# ============================================================
# STAGE INFERENTIAL
# ============================================================

stage_eligible <- copy(
  patient[!is.na(stage_num)]
)

stage_eligible[
  ,
  stage_group_n := .N,
  by = .(
    gene,
    project,
    stage_simple,
    stage_num
  )
]

stage_eligible <-
  stage_eligible[
    stage_group_n >= 5
  ]

stage_projects <- stage_eligible[
  ,
  .(
    groups_used =
      uniqueN(stage_num),
    patients_used = .N
  ),
  by = .(
    gene,
    project
  )
][
  groups_used >= 2
]

stage_list <- vector(
  "list",
  nrow(stage_projects)
)

for (i in seq_len(nrow(stage_projects))) {

  g <- stage_projects$gene[i]
  cancer <- stage_projects$project[i]

  x <- stage_eligible[
    gene == g &
    project == cancer
  ]

  sp <- safe_spearman(
    x$stage_num,
    x$log2TPM1
  )

  original <- stage_desc[
    gene == g &
    project == cancer
  ]

  stage_list[[i]] <- data.table(

    gene = g,
    project = cancer,

    Patients_used =
      nrow(x),

    Stage_groups_used =
      uniqueN(x$stage_num),

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
      sp["p"]
  )
}

stage_results <- rbindlist(
  stage_list,
  fill = TRUE
)

stage_results[, Kruskal_Wallis_FDR :=
  p.adjust(
    Kruskal_Wallis_P,
    method = "BH"
  ),
  by = gene
]

stage_results[, Stage_trend_FDR :=
  p.adjust(
    Spearman_P,
    method = "BH"
  ),
  by = gene
]

stage_results[, Global_KW_FDR :=
  p.adjust(
    Kruskal_Wallis_P,
    method = "BH"
  )
]

stage_results[, Global_stage_trend_FDR :=
  p.adjust(
    Spearman_P,
    method = "BH"
  )
]

setorder(
  stage_results,
  gene,
  Stage_trend_FDR
)

fwrite(
  stage_results,
  file.path(
    outdir,
    "06_stage_statistics.csv"
  )
)

cat(
  "3. Stage analysis complete:",
  nrow(stage_results),
  "gene-cancer tests\n"
)

# ============================================================
# GRADE LABEL AUDIT
# ============================================================

grade_excluded <- patient[
  !is.na(grade_raw) &
  is.na(grade_num),
  .(
    patients = .N
  ),
  by = .(
    gene,
    project,
    grade_raw
  )
]

fwrite(
  grade_excluded,
  file.path(
    outdir,
    "06_grade_labels_excluded.csv"
  )
)

# ============================================================
# GRADE DESCRIPTIVE
# ============================================================

grade_desc <- patient[
  !is.na(grade_num),
  .(
    n = .N,

    median_TPM =
      median(TPM, na.rm = TRUE),

    median_log2TPM1 =
      median(log2TPM1, na.rm = TRUE)
  ),
  by = .(
    gene,
    project,
    grade_harmonized,
    grade_num
  )
]

setorder(
  grade_desc,
  gene,
  project,
  grade_num
)

fwrite(
  grade_desc,
  file.path(
    outdir,
    "06_grade_descriptive.csv"
  )
)

# ============================================================
# GRADE INFERENTIAL
# ============================================================

grade_eligible <-
  copy(
    patient[!is.na(grade_num)]
  )

grade_eligible[
  ,
  grade_group_n := .N,
  by = .(
    gene,
    project,
    grade_harmonized,
    grade_num
  )
]

grade_eligible <-
  grade_eligible[
    grade_group_n >= 5
  ]

grade_projects <- grade_eligible[
  ,
  .(
    groups_used =
      uniqueN(grade_num),
    patients_used = .N
  ),
  by = .(
    gene,
    project
  )
][
  groups_used >= 2
]

grade_list <- vector(
  "list",
  nrow(grade_projects)
)

for (i in seq_len(nrow(grade_projects))) {

  g <- grade_projects$gene[i]
  cancer <- grade_projects$project[i]

  x <- grade_eligible[
    gene == g &
    project == cancer
  ]

  sp <- safe_spearman(
    x$grade_num,
    x$log2TPM1
  )

  original <- grade_desc[
    gene == g &
    project == cancer
  ]

  excluded_labels <- grade_excluded[
    gene == g &
    project == cancer
  ]

  grade_list[[i]] <- data.table(

    gene = g,
    project = cancer,

    Patients_used =
      nrow(x),

    Grade_groups_used =
      uniqueN(x$grade_num),

    Groups_used =
      paste(
        sort(
          unique(
            x$grade_harmonized
          )
        ),
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
      if (
        nrow(excluded_labels) == 0
      ) {
        ""
      } else {
        paste(
          excluded_labels$grade_raw,
          excluded_labels$patients,
          sep = ":n=",
          collapse = "; "
        )
      },

    Kruskal_Wallis_P =
      safe_kw(
        x$log2TPM1,
        x$grade_num
      ),

    Spearman_rho =
      sp["rho"],

    Spearman_P =
      sp["p"]
  )
}

grade_results <- rbindlist(
  grade_list,
  fill = TRUE
)

grade_results[, Kruskal_Wallis_FDR :=
  p.adjust(
    Kruskal_Wallis_P,
    method = "BH"
  ),
  by = gene
]

grade_results[, Grade_trend_FDR :=
  p.adjust(
    Spearman_P,
    method = "BH"
  ),
  by = gene
]

grade_results[, Global_KW_FDR :=
  p.adjust(
    Kruskal_Wallis_P,
    method = "BH"
  )
]

grade_results[, Global_grade_trend_FDR :=
  p.adjust(
    Spearman_P,
    method = "BH"
  )
]

setorder(
  grade_results,
  gene,
  Grade_trend_FDR
)

fwrite(
  grade_results,
  file.path(
    outdir,
    "06_grade_statistics.csv"
  )
)

cat(
  "4. Grade analysis complete:",
  nrow(grade_results),
  "gene-cancer tests\n"
)

# ============================================================
# GENE-LEVEL COMPARISON
# ============================================================

stage_summary <- stage_results[
  ,
  .(
    Stage_cancers_tested = .N,

    Stage_KW_significant =
      sum(
        Kruskal_Wallis_FDR < 0.05,
        na.rm = TRUE
      ),

    Stage_trend_significant =
      sum(
        Stage_trend_FDR < 0.05,
        na.rm = TRUE
      ),

    Positive_stage_trends =
      sum(
        Stage_trend_FDR < 0.05 &
        Spearman_rho > 0,
        na.rm = TRUE
      ),

    Negative_stage_trends =
      sum(
        Stage_trend_FDR < 0.05 &
        Spearman_rho < 0,
        na.rm = TRUE
      )
  ),
  by = gene
]

grade_summary <- grade_results[
  ,
  .(
    Grade_cancers_tested = .N,

    Grade_KW_significant =
      sum(
        Kruskal_Wallis_FDR < 0.05,
        na.rm = TRUE
      ),

    Grade_trend_significant =
      sum(
        Grade_trend_FDR < 0.05,
        na.rm = TRUE
      ),

    Positive_grade_trends =
      sum(
        Grade_trend_FDR < 0.05 &
        Spearman_rho > 0,
        na.rm = TRUE
      ),

    Negative_grade_trends =
      sum(
        Grade_trend_FDR < 0.05 &
        Spearman_rho < 0,
        na.rm = TRUE
      )
  ),
  by = gene
]

summary_gene <- merge(
  stage_summary,
  grade_summary,
  by = "gene",
  all = TRUE
)

summary_gene[, Reference :=
  gene == "JAK3"
]

fwrite(
  summary_gene,
  file.path(
    outdir,
    "06_stage_grade_gene_comparison_summary.csv"
  )
)

# ============================================================
# Console
# ============================================================

cat("\n========================================\n")
cat("MULTIKINASE STAGE / GRADE ANALYSIS\n")
cat("========================================\n")

cat("\n===== GENE-LEVEL COMPARISON =====\n")
print(summary_gene)

cat("\n===== JAK3 STAGE RESULTS =====\n")

print(
  stage_results[
    gene == "JAK3"
  ]
)

cat("\n===== JAK3 GRADE RESULTS =====\n")

print(
  grade_results[
    gene == "JAK3"
  ]
)

cat("\n06 COMPLETE\n")

