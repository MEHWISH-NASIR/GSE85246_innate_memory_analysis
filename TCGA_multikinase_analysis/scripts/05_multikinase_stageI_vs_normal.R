
library(data.table)
setDTthreads(1)

root <- "TCGA_multikinase_analysis"

outdir <- file.path(root, "results/early_stage")
auditdir <- file.path(root, "audits")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(auditdir, recursive = TRUE, showWarnings = FALSE)

cat("1. Loading complete multikinase dataset...\n")

d <- fread(
  file.path(
    root,
    "results/expression",
    "03_multikinase_expression_clinical_complete.csv"
  ),
  nThread = 1
)

# ============================================================
# Helper functions
# ============================================================

first_nonmissing <- function(x) {

  x <- x[
    !is.na(x) &
    trimws(as.character(x)) != ""
  ]

  if (length(x) == 0) {
    NA_character_
  } else {
    as.character(x[1])
  }
}

safe_wilcox_unpaired <- function(x, y) {

  tryCatch(
    wilcox.test(
      x,
      y,
      paired = FALSE,
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

# ============================================================
# Canonical sample classes
# Same as completed JAK3 analysis
# ============================================================

d[, analysis_class := fifelse(
  project == "TCGA-LAML" &
    sample_type ==
    "Primary Blood Derived Cancer - Peripheral Blood",
  "Primary malignant",

  fifelse(
    project != "TCGA-LAML" &
      sample_type == "Primary Tumor",
    "Primary malignant",

    fifelse(
      sample_type == "Solid Tissue Normal",
      "Normal",
      "Other"
    )
  )
)]

cat("2. Sample classes assigned\n")

# ============================================================
# Collapse primary tumors to patient level
# Preserve AJCC stage
# ============================================================

tumor <- d[
  analysis_class == "Primary malignant",
  .(
    TPM = mean(TPM, na.rm = TRUE),

    stage = first_nonmissing(
      ajcc_pathologic_stage
    ),

    expression_samples = .N
  ),
  by = .(
    gene,
    project,
    patient_id
  )
]

tumor[, log2TPM1 := log2(TPM + 1)]

tumor[, stage_clean :=
  toupper(trimws(stage))
]

tumor[, is_stage_I :=
  !is.na(stage_clean) &
  grepl(
    "^STAGE I([ABC])?$",
    stage_clean
  )
]

# ============================================================
# Collapse normals to patient level
# ============================================================

normal <- d[
  analysis_class == "Normal",
  .(
    TPM = mean(TPM, na.rm = TRUE),
    expression_samples = .N
  ),
  by = .(
    gene,
    project,
    patient_id
  )
]

normal[, log2TPM1 := log2(TPM + 1)]

cat("3. Patient-level tumor/normal datasets prepared\n")

# Save patient-level stage table
fwrite(
  tumor,
  file.path(
    outdir,
    "05_primary_tumor_patient_level_with_stage.csv"
  )
)

# ============================================================
# Identify cancers containing normal tissue
# ============================================================

eligible <- merge(
  unique(tumor[, .(gene, project)]),
  unique(normal[, .(gene, project)]),
  by = c("gene", "project")
)

cat(
  "4. Eligible gene-cancer combinations:",
  nrow(eligible),
  "\n"
)

# ============================================================
# Stage I vs normal
# ============================================================

result_list <- vector(
  "list",
  nrow(eligible)
)

for (i in seq_len(nrow(eligible))) {

  g <- eligible$gene[i]
  cancer <- eligible$project[i]

  e <- tumor[
    gene == g &
    project == cancer &
    is_stage_I == TRUE
  ]

  n <- normal[
    gene == g &
    project == cancer
  ]

  if (nrow(e) == 0) {
    next
  }

  delta <-
    median(e$log2TPM1, na.rm = TRUE) -
    median(n$log2TPM1, na.rm = TRUE)

  result_list[[i]] <- data.table(

    gene = g,

    project = cancer,

    Early_stage_definition =
      "AJCC Stage I/IA/IB/IC only",

    Stage_I_n = nrow(e),

    normal_n = nrow(n),

    Stage_I_median_TPM =
      median(e$TPM, na.rm = TRUE),

    normal_median_TPM =
      median(n$TPM, na.rm = TRUE),

    Median_log2TPM1_difference =
      delta,

    Direction =
      ifelse(
        median(
          e$log2TPM1,
          na.rm = TRUE
        ) >
        median(
          n$log2TPM1,
          na.rm = TRUE
        ),
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
      nrow(n) < 10
  )

  if (i %% 25 == 0) {
    cat(
      "Stage-I groups processed:",
      i, "/", nrow(eligible), "\n"
    )
  }
}

early <- rbindlist(
  result_list,
  fill = TRUE
)

# ============================================================
# Same gene-specific BH framework used for JAK3
# ============================================================

early[, Wilcoxon_FDR :=
  p.adjust(
    Wilcoxon_P,
    method = "BH"
  ),
  by = gene
]

# Additional multikinase-wide correction
early[, Multikinase_Global_FDR :=
  p.adjust(
    Wilcoxon_P,
    method = "BH"
  )
]

setorder(
  early,
  gene,
  Wilcoxon_FDR,
  project
)

fwrite(
  early,
  file.path(
    outdir,
    "05_multikinase_stageI_vs_normal.csv"
  )
)

# ============================================================
# Summary comparison against JAK3
# ============================================================

summary_gene <- early[, .(

  Cancers_tested = .N,

  Significant_FDR_005 =
    sum(
      Wilcoxon_FDR < 0.05,
      na.rm = TRUE
    ),

  Significant_global_FDR_005 =
    sum(
      Multikinase_Global_FDR < 0.05,
      na.rm = TRUE
    ),

  Significant_UP =
    sum(
      Wilcoxon_FDR < 0.05 &
      Direction == "UP",
      na.rm = TRUE
    ),

  Significant_DOWN =
    sum(
      Wilcoxon_FDR < 0.05 &
      Direction == "DOWN",
      na.rm = TRUE
    ),

  Median_effect_across_cancers =
    median(
      Median_log2TPM1_difference,
      na.rm = TRUE
    )

), by = gene]

summary_gene[, Reference :=
  gene == "JAK3"
]

setorder(
  summary_gene,
  -Significant_FDR_005,
  -Significant_global_FDR_005
)

fwrite(
  summary_gene,
  file.path(
    outdir,
    "05_stageI_gene_comparison_summary.csv"
  )
)

# ============================================================
# JAK3 preview
# ============================================================

cat("\n========================================\n")
cat("MULTIKINASE STAGE-I VS NORMAL ANALYSIS\n")
cat("========================================\n")

cat(
  "Genes represented:",
  uniqueN(early$gene),
  "\n"
)

cat(
  "Stage-I comparisons:",
  nrow(early),
  "\n"
)

cat("\n===== JAK3 STAGE-I RESULTS =====\n")

print(
  early[
    gene == "JAK3",
    .(
      project,
      Stage_I_n,
      normal_n,
      Stage_I_median_TPM,
      normal_median_TPM,
      Median_log2TPM1_difference,
      Direction,
      Wilcoxon_P,
      Wilcoxon_FDR
    )
  ]
)

cat("\n===== GENE-LEVEL COMPARISON =====\n")

print(summary_gene)

cat("\n05 COMPLETE\n")

