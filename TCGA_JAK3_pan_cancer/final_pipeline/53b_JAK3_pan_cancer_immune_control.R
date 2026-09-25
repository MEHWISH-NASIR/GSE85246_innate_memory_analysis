library(dplyr)

# ============================================================
# 53b — JAK3 PAN-CANCER IMMUNE-MARKER CONTROL
#
# Questions:
# 1. Does JAK3 correlate with immune markers in tumors?
# 2. Does JAK3 correlate with the 10-marker immune score?
# 3. Does tumor-vs-normal elevation persist after adjustment?
# 4. Does Stage-I-vs-normal elevation persist after adjustment?
#
# Immune score:
# Mean of z-scored log2(TPM+1):
# PTPRC, CD3D, CD3E, MS4A1, CD79A,
# NKG7, LST1, TYROBP, FCER1G, CD68
#
# Patient-level duplicate samples are collapsed by MEAN,
# consistent with the current pan-cancer workflow.
# ============================================================

markers <- c(
  "PTPRC",
  "CD3D",
  "CD3E",
  "MS4A1",
  "CD79A",
  "NKG7",
  "LST1",
  "TYROBP",
  "FCER1G",
  "CD68"
)

immune_file <-
  "TCGA_JAK3_pan_cancer/results/immune_markers/TCGA_pan_cancer_10_immune_markers_TPM.csv"

jak3_file <-
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_expression_clinical_complete_33projects.csv"

outdir <-
  "TCGA_JAK3_pan_cancer/results/immune_markers"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(immune_file),
  file.exists(jak3_file)
)

# ============================================================
# HELPERS
# ============================================================

find_col <- function(dat, candidates, label) {

  hit <- candidates[candidates %in% names(dat)]

  if (length(hit) == 0) {
    stop(
      paste(
        "Could not identify column:",
        label,
        "\nAvailable columns:\n",
        paste(names(dat), collapse = ", ")
      )
    )
  }

  hit[1]
}

first_nonmissing <- function(x) {

  x <- as.character(x)

  x <- x[
    !is.na(x) &
    trimws(x) != ""
  ]

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[1]
}

normalize_stage <- function(x) {

  x <- toupper(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA_character_,
    length(x)
  )

  out[
    grepl("^STAGE IV([ABC])?$", x)
  ] <- "Stage IV"

  out[
    grepl("^STAGE III([ABC])?$", x)
  ] <- "Stage III"

  out[
    grepl("^STAGE II([ABC])?$", x)
  ] <- "Stage II"

  out[
    grepl("^STAGE I([ABC])?$", x)
  ] <- "Stage I"

  out
}

safe_spearman <- function(x, y) {

  z <- complete.cases(x, y)

  n <- sum(z)

  if (
    n < 10 ||
    sd(x[z]) == 0 ||
    sd(y[z]) == 0
  ) {

    return(
      c(
        N = n,
        rho = NA_real_,
        P = NA_real_
      )
    )
  }

  ct <- suppressWarnings(
    cor.test(
      x[z],
      y[z],
      method = "spearman",
      exact = FALSE
    )
  )

  c(
    N = n,
    rho = unname(ct$estimate),
    P = ct$p.value
  )
}

# ============================================================
# READ DATA
# ============================================================

immune <- read.csv(
  immune_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

jak3 <- read.csv(
  jak3_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("JAK3 PAN-CANCER IMMUNE CONTROL\n")
cat("========================================\n")

cat("\nImmune-marker rows:", nrow(immune), "\n")
cat("JAK3 rows:", nrow(jak3), "\n")

# ============================================================
# IDENTIFY COLUMNS
# ============================================================

i_project <- find_col(
  immune,
  c(
    "project",
    "TCGA_project",
    "project_id",
    "cases.project.project_id"
  ),
  "immune project"
)

i_patient <- find_col(
  immune,
  c(
    "cases.submitter_id",
    "patient_id",
    "cases"
  ),
  "immune patient"
)

i_sample_type <- find_col(
  immune,
  c(
    "sample_type",
    "samples.sample_type"
  ),
  "immune sample type"
)

i_tpm <- find_col(
  immune,
  c("TPM"),
  "immune TPM"
)

j_project <- find_col(
  jak3,
  c(
    "project",
    "TCGA_project",
    "project_id",
    "cases.project.project_id"
  ),
  "JAK3 project"
)

j_patient <- find_col(
  jak3,
  c(
    "cases.submitter_id",
    "patient_id",
    "cases"
  ),
  "JAK3 patient"
)

j_sample_type <- find_col(
  jak3,
  c(
    "sample_type",
    "samples.sample_type"
  ),
  "JAK3 sample type"
)

j_tpm <- find_col(
  jak3,
  c(
    "TPM",
    "JAK3_TPM"
  ),
  "JAK3 TPM"
)

j_stage <- find_col(
  jak3,
  c(
    "ajcc_pathologic_stage",
    "pathologic_stage",
    "tumor_stage",
    "stage"
  ),
  "pathologic stage"
)

# ============================================================
# STRICT IMMUNE TABLE AUDIT
# ============================================================

if (nrow(immune) != 115050) {
  stop("Immune table does not contain exactly 115,050 rows.")
}

if (
  "extraction_status" %in% names(immune) &&
  any(immune$extraction_status != "OK")
) {
  stop("Non-OK immune-marker extraction records detected.")
}

immune_key <- paste(
  immune$file_name,
  immune$gene_name,
  sep = "__"
)

if (anyDuplicated(immune_key)) {
  stop("Duplicated file_name + gene_name records detected.")
}

# ============================================================
# STANDARDIZE IMMUNE TABLE
# ============================================================

immune_std <- data.frame(

  project =
    as.character(
      immune[[i_project]]
    ),

  patient_id =
    as.character(
      immune[[i_patient]]
    ),

  sample_type =
    as.character(
      immune[[i_sample_type]]
    ),

  gene_name =
    as.character(
      immune$gene_name
    ),

  TPM =
    as.numeric(
      immune[[i_tpm]]
    ),

  stringsAsFactors = FALSE
)

immune_std <- immune_std[
  immune_std$gene_name %in% markers,
]

# ============================================================
# PATIENT-LEVEL IMMUNE MARKERS
#
# Current pan-cancer convention:
# duplicate aliquots -> MEAN
# ============================================================

immune_patient <- immune_std %>%
  group_by(
    project,
    patient_id,
    sample_type,
    gene_name
  ) %>%
  summarise(
    TPM = mean(TPM, na.rm = TRUE),
    .groups = "drop"
  )

base_keys <- unique(
  immune_patient[
    ,
    c(
      "project",
      "patient_id",
      "sample_type"
    )
  ]
)

immune_wide <- base_keys

for (g in markers) {

  z <- immune_patient[
    immune_patient$gene_name == g,
    c(
      "project",
      "patient_id",
      "sample_type",
      "TPM"
    )
  ]

  names(z)[4] <- g

  immune_wide <- merge(
    immune_wide,
    z,
    by = c(
      "project",
      "patient_id",
      "sample_type"
    ),
    all.x = TRUE
  )
}

# ============================================================
# LOG2(TPM+1)
# ============================================================

log_marker_names <- paste0(
  markers,
  "_log2TPM1"
)

for (i in seq_along(markers)) {

  immune_wide[[log_marker_names[i]]] <-
    log2(
      immune_wide[[markers[i]]] + 1
    )
}

# ============================================================
# IMMUNE SCORE
#
# IMPORTANT:
# z-standardization is performed WITHIN EACH CANCER,
# using its full Primary Tumor + Solid Tissue Normal cohort.
#
# This preserves the logic of the old KIRC analysis while
# avoiding cross-cancer tissue-expression differences.
# ============================================================

immune_wide$ImmuneScore <- NA_real_

for (p in unique(immune_wide$project)) {

  idx <-
    immune_wide$project == p &
    immune_wide$sample_type %in%
      c(
        "Primary Tumor",
        "Solid Tissue Normal"
      )

  zmat <- matrix(
    NA_real_,
    nrow = sum(idx),
    ncol = length(markers)
  )

  for (i in seq_along(markers)) {

    x <-
      immune_wide[
        idx,
        log_marker_names[i]
      ]

    if (
      sum(!is.na(x)) > 1 &&
      sd(x, na.rm = TRUE) > 0
    ) {

      zmat[, i] <-
        as.numeric(
          scale(x)
        )
    }
  }

  immune_wide$ImmuneScore[idx] <-
    rowMeans(
      zmat,
      na.rm = TRUE
    )
}

immune_wide$PTPRC_log2TPM1 <-
  log2(
    immune_wide$PTPRC + 1
  )

write.csv(
  immune_wide,
  file.path(
    outdir,
    "TCGA_pan_cancer_immune_marker_patient_data.csv"
  ),
  row.names = FALSE
)

# ============================================================
# STANDARDIZE JAK3 TABLE
# ============================================================

jak3_std <- data.frame(

  project =
    as.character(
      jak3[[j_project]]
    ),

  patient_id =
    as.character(
      jak3[[j_patient]]
    ),

  sample_type =
    as.character(
      jak3[[j_sample_type]]
    ),

  JAK3_TPM =
    as.numeric(
      jak3[[j_tpm]]
    ),

  Stage_raw =
    as.character(
      jak3[[j_stage]]
    ),

  stringsAsFactors = FALSE
)

# ============================================================
# PATIENT-LEVEL JAK3
# ============================================================

jak3_patient <- jak3_std %>%
  group_by(
    project,
    patient_id,
    sample_type
  ) %>%
  summarise(

    JAK3_TPM =
      mean(
        JAK3_TPM,
        na.rm = TRUE
      ),

    Stage_raw =
      first_nonmissing(
        Stage_raw
      ),

    .groups = "drop"
  )

jak3_patient$JAK3_log2TPM1 <-
  log2(
    jak3_patient$JAK3_TPM + 1
  )

jak3_patient$Stage_group <-
  normalize_stage(
    jak3_patient$Stage_raw
  )

# ============================================================
# MERGE JAK3 + IMMUNE MARKERS
# ============================================================

dat <- merge(
  jak3_patient,
  immune_wide,
  by = c(
    "project",
    "patient_id",
    "sample_type"
  ),
  all = FALSE
)

cat(
  "Merged patient/sample-type records:",
  nrow(dat),
  "\n"
)

cat(
  "Projects represented:",
  length(unique(dat$project)),
  "\n"
)

# ============================================================
# IMMUNE MARKER CORRELATIONS
#
# Primary tumor and Stage-I tumor contexts
# ============================================================

cor_results <- list()

features <- c(
  log_marker_names,
  "ImmuneScore"
)

feature_labels <- c(
  markers,
  "ImmuneScore"
)

for (p in sort(unique(dat$project))) {

  dp <- dat[
    dat$project == p,
  ]

  contexts <- list(

    "Primary tumors" =
      dp[
        dp$sample_type ==
          "Primary Tumor",
      ],

    "Stage I tumors" =
      dp[
        dp$sample_type ==
          "Primary Tumor" &
        dp$Stage_group ==
          "Stage I",
      ]
  )

  for (ctx in names(contexts)) {

    d <- contexts[[ctx]]

    for (i in seq_along(features)) {

      rr <- safe_spearman(
        d$JAK3_log2TPM1,
        d[[features[i]]]
      )

      cor_results[[
        length(cor_results) + 1
      ]] <- data.frame(

        project = p,
        Context = ctx,
        Feature = feature_labels[i],
        N = rr["N"],
        Spearman_rho = rr["rho"],
        P = rr["P"],

        stringsAsFactors = FALSE
      )
    }
  }
}

cor_tab <- bind_rows(
  cor_results
)

cor_tab <- cor_tab %>%
  group_by(
    project,
    Context
  ) %>%
  mutate(
    Within_cancer_context_FDR =
      p.adjust(
        P,
        method = "BH"
      )
  ) %>%
  ungroup() %>%
  group_by(
    Context,
    Feature
  ) %>%
  mutate(
    Feature_family_FDR =
      p.adjust(
        P,
        method = "BH"
      )
  ) %>%
  ungroup()

cor_tab$Global_FDR <-
  p.adjust(
    cor_tab$P,
    method = "BH"
  )

write.csv(
  cor_tab,
  file.path(
    outdir,
    "JAK3_pan_cancer_immune_marker_correlations.csv"
  ),
  row.names = FALSE
)

# ============================================================
# LINEAR MODEL HELPER
# ============================================================

run_model <- function(
  d,
  comparison,
  adjustment
) {

  tumor_n <-
    sum(
      d$Group == 1,
      na.rm = TRUE
    )

  normal_n <-
    sum(
      d$Group == 0,
      na.rm = TRUE
    )

  if (
    tumor_n < 1 ||
    normal_n < 1 ||
    nrow(d) < 3
  ) {
    return(NULL)
  }

  form <-
    switch(

      adjustment,

      "Unadjusted" =
        JAK3_log2TPM1 ~ Group,

      "Adjusted for PTPRC/CD45" =
        JAK3_log2TPM1 ~
          Group +
          PTPRC_log2TPM1,

      "Adjusted for multi-marker immune score" =
        JAK3_log2TPM1 ~
          Group +
          ImmuneScore
    )

  fit <- tryCatch(
    lm(
      form,
      data = d
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(NULL)
  }

  sm <- summary(fit)$coefficients

  if (!("Group" %in% rownames(sm))) {
    return(NULL)
  }

  data.frame(

    Comparison =
      comparison,

    Adjustment =
      adjustment,

    N =
      nobs(fit),

    Tumor_n =
      tumor_n,

    Normal_n =
      normal_n,

    Coefficient =
      sm[
        "Group",
        "Estimate"
      ],

    SE =
      sm[
        "Group",
        "Std. Error"
      ],

    T =
      sm[
        "Group",
        "t value"
      ],

    P =
      sm[
        "Group",
        "Pr(>|t|)"
      ],

    Low_N =
      tumor_n < 5 |
      normal_n < 5,

    stringsAsFactors = FALSE
  )
}

# ============================================================
# PAN-CANCER ADJUSTED MODELS
# ============================================================

model_results <- list()

adjustments <- c(
  "Unadjusted",
  "Adjusted for PTPRC/CD45",
  "Adjusted for multi-marker immune score"
)

for (p in sort(unique(dat$project))) {

  d <- dat[
    dat$project == p &
    dat$sample_type %in%
      c(
        "Primary Tumor",
        "Solid Tissue Normal"
      ),
  ]

  # ----------------------------------------------------------
  # ALL PRIMARY TUMORS VS NORMAL
  # ----------------------------------------------------------

  if (
    any(
      d$sample_type ==
        "Primary Tumor"
    ) &&
    any(
      d$sample_type ==
        "Solid Tissue Normal"
    )
  ) {

    tn <- d

    tn$Group <-
      ifelse(
        tn$sample_type ==
          "Primary Tumor",
        1,
        0
      )

    for (adj in adjustments) {

      rr <- run_model(
        tn,
        "All primary tumors vs normal",
        adj
      )

      if (!is.null(rr)) {

        rr$project <- p

        model_results[[
          length(model_results) + 1
        ]] <- rr
      }
    }
  }

  # ----------------------------------------------------------
  # STAGE I TUMORS VS NORMAL
  # ----------------------------------------------------------

  early <- d[
    d$sample_type ==
      "Solid Tissue Normal" |
    (
      d$sample_type ==
        "Primary Tumor" &
      !is.na(d$Stage_group) &
      d$Stage_group ==
        "Stage I"
    ),
  ]

  if (
    any(
      early$sample_type ==
        "Primary Tumor"
    ) &&
    any(
      early$sample_type ==
        "Solid Tissue Normal"
    )
  ) {

    early$Group <-
      ifelse(
        early$sample_type ==
          "Primary Tumor",
        1,
        0
      )

    for (adj in adjustments) {

      rr <- run_model(
        early,
        "Stage I tumors vs normal",
        adj
      )

      if (!is.null(rr)) {

        rr$project <- p

        model_results[[
          length(model_results) + 1
        ]] <- rr
      }
    }
  }
}

model_tab <- bind_rows(
  model_results
)

# ============================================================
# MULTIPLE-TESTING CORRECTION
#
# 1. Within-cancer:
#    preserves old KIRC six-model logic.
#
# 2. Model-family:
#    same comparison + adjustment across cancers.
#
# 3. Global:
#    all pan-cancer model tests together.
# ============================================================

model_tab <- model_tab %>%
  group_by(
    project
  ) %>%
  mutate(
    Within_cancer_FDR =
      p.adjust(
        P,
        method = "BH"
      )
  ) %>%
  ungroup() %>%
  group_by(
    Comparison,
    Adjustment
  ) %>%
  mutate(
    Model_family_FDR =
      p.adjust(
        P,
        method = "BH"
      )
  ) %>%
  ungroup()

model_tab$Global_FDR <-
  p.adjust(
    model_tab$P,
    method = "BH"
  )

model_tab <- model_tab %>%
  select(
    project,
    Comparison,
    Adjustment,
    N,
    Tumor_n,
    Normal_n,
    Coefficient,
    SE,
    T,
    P,
    Within_cancer_FDR,
    Model_family_FDR,
    Global_FDR,
    Low_N
  ) %>%
  arrange(
    Comparison,
    Adjustment,
    Model_family_FDR
  )

write.csv(
  model_tab,
  file.path(
    outdir,
    "JAK3_pan_cancer_immune_adjusted_models.csv"
  ),
  row.names = FALSE
)

# ============================================================
# CONSOLE SUMMARY
# ============================================================

tn_projects <-
  unique(
    model_tab$project[
      model_tab$Comparison ==
        "All primary tumors vs normal"
    ]
  )

early_projects <-
  unique(
    model_tab$project[
      model_tab$Comparison ==
        "Stage I tumors vs normal"
    ]
  )

cat("\n========================================\n")
cat("53b RESULTS SUMMARY\n")
cat("========================================\n")

cat(
  "\nTumor-normal cancers analyzed:",
  length(tn_projects),
  "\n"
)

cat(
  "Stage-I-normal cancers analyzed:",
  length(early_projects),
  "\n"
)

cat(
  "\nPTPRC-adjusted model-family FDR < 0.05:",
  sum(
    model_tab$Adjustment ==
      "Adjusted for PTPRC/CD45" &
    model_tab$Model_family_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "ImmuneScore-adjusted model-family FDR < 0.05:",
  sum(
    model_tab$Adjustment ==
      "Adjusted for multi-marker immune score" &
    model_tab$Model_family_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat("\n===== SIGNIFICANT IMMUNE-ADJUSTED MODELS =====\n")

sig <- model_tab[
  model_tab$Adjustment != "Unadjusted" &
  !is.na(model_tab$Model_family_FDR) &
  model_tab$Model_family_FDR < 0.05,
]

if (nrow(sig) > 0) {

  print(
    as.data.frame(
      sig[
        ,
        c(
          "project",
          "Comparison",
          "Adjustment",
          "Tumor_n",
          "Normal_n",
          "Coefficient",
          "P",
          "Model_family_FDR",
          "Global_FDR",
          "Low_N"
        )
      ]
    ),
    row.names = FALSE
  )

} else {

  cat("None\n")
}

cat("\nSaved:\n")
cat(
  " ",
  file.path(
    outdir,
    "TCGA_pan_cancer_immune_marker_patient_data.csv"
  ),
  "\n"
)

cat(
  " ",
  file.path(
    outdir,
    "JAK3_pan_cancer_immune_marker_correlations.csv"
  ),
  "\n"
)

cat(
  " ",
  file.path(
    outdir,
    "JAK3_pan_cancer_immune_adjusted_models.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("53b COMPLETE\n")
cat("========================================\n")
