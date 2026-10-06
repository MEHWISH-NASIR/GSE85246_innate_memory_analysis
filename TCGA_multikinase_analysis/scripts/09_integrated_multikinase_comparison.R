
library(data.table)
setDTthreads(1)

root <- "TCGA_multikinase_analysis"

outdir <- file.path(
  root,
  "results/integrated"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n========================================\n")
cat("INTEGRATED MULTIKINASE COMPARISON\n")
cat("========================================\n")

# ============================================================
# 1. ORIGINAL INNATE-MEMORY CANDIDATE TABLE
# ============================================================

memory <- fread(
  "results/15_multikinase_comparison/15_initial_TCGA_kinase_panel.csv"
)

setnames(
  memory,
  old = "SYMBOL",
  new = "gene"
)

keep_memory <- intersect(
  c(
    "gene",
    "strict_persistent",
    "log2FC_Day6",
    "FDR_Day6",
    "Downstream_status",
    "selection_group",
    "reference_gene"
  ),
  names(memory)
)

memory <- memory[
  ,
  ..keep_memory
]

cat(
  "Original candidate genes:",
  uniqueN(memory$gene),
  "\n"
)

# ============================================================
# 2. PRIMARY TUMOR-NORMAL
# ============================================================

tn <- fread(
  file.path(
    root,
    "results/primary",
    "04_multikinase_tumor_normal.csv"
  )
)

tn_summary <- tn[
  ,
  .(
    TN_cancers_tested = .N,

    TN_significant =
      sum(
        Wilcoxon_FDR < 0.05,
        na.rm = TRUE
      ),

    TN_UP_significant =
      sum(
        Wilcoxon_FDR < 0.05 &
        Direction == "UP",
        na.rm = TRUE
      ),

    TN_DOWN_significant =
      sum(
        Wilcoxon_FDR < 0.05 &
        Direction == "DOWN",
        na.rm = TRUE
      ),

    TN_median_effect =
      median(
        Median_log2TPM1_difference,
        na.rm = TRUE
      )
  ),
  by = gene
]

# ============================================================
# 3. MATCHED TUMOR-NORMAL
# ============================================================

matched <- fread(
  file.path(
    root,
    "results/primary",
    "04_multikinase_matched_tumor_normal.csv"
  )
)

matched_summary <- matched[
  ,
  .(
    Matched_cancers_tested = .N,

    Matched_significant =
      sum(
        Paired_Wilcoxon_FDR < 0.05,
        na.rm = TRUE
      ),

    Matched_positive_effect =
      sum(
        Paired_Wilcoxon_FDR < 0.05 &
        median_log2_difference > 0,
        na.rm = TRUE
      ),

    Matched_negative_effect =
      sum(
        Paired_Wilcoxon_FDR < 0.05 &
        median_log2_difference < 0,
        na.rm = TRUE
      ),

    Matched_median_effect =
      median(
        median_log2_difference,
        na.rm = TRUE
      )
  ),
  by = gene
]

# ============================================================
# 4. EARLY-STAGE / STAGE-I VS NORMAL
# ============================================================

early <- fread(
  file.path(
    root,
    "results/early_stage",
    "05_stageI_gene_comparison_summary.csv"
  )
)

setnames(
  early,
  old = c(
    "Cancers_tested",
    "Significant_FDR_005",
    "Significant_global_FDR_005",
    "Significant_UP",
    "Significant_DOWN",
    "Median_effect_across_cancers"
  ),
  new = c(
    "StageI_cancers_tested",
    "StageI_significant",
    "StageI_global_significant",
    "StageI_UP_significant",
    "StageI_DOWN_significant",
    "StageI_median_effect"
  ),
  skip_absent = TRUE
)

if ("Reference" %in% names(early)) {
  early[, Reference := NULL]
}

# ============================================================
# 5. STAGE / GRADE PROGRESSION
# ============================================================

progression <- fread(
  file.path(
    root,
    "results/stage_grade",
    "06_stage_grade_gene_comparison_summary.csv"
  )
)

if ("Reference" %in% names(progression)) {
  progression[, Reference := NULL]
}

# ============================================================
# 6. SURVIVAL
# ============================================================

survival <- fread(
  file.path(
    root,
    "results/survival",
    "07_survival_gene_comparison_summary.csv"
  )
)

if ("Reference" %in% names(survival)) {
  survival[, Reference := NULL]
}

# ============================================================
# 7. PANIMMUNE-ADJUSTED RESULTS
# ============================================================

immune <- fread(
  file.path(
    root,
    "results/immune",
    "08_PanImmune_gene_comparison_summary.csv"
  )
)

if ("Reference" %in% names(immune)) {
  immune[, Reference := NULL]
}

# ============================================================
# 8. MERGE ALL EVIDENCE
# ============================================================

tables <- list(
  memory,
  tn_summary,
  matched_summary,
  early,
  progression,
  survival,
  immune
)

integrated <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "gene",
      all = TRUE
    )
  },
  tables
)

integrated[, Reference :=
  gene == "JAK3"
]

# ============================================================
# 9. CALCULATE SUPPORT PROPORTIONS
# ============================================================

safe_prop <- function(num, den) {

  out <- rep(
    NA_real_,
    length(num)
  )

  ok <-
    !is.na(num) &
    !is.na(den) &
    den > 0

  out[ok] <-
    num[ok] / den[ok]

  out
}

# Primary tumor-normal
integrated[, TN_support :=
  safe_prop(
    TN_significant,
    TN_cancers_tested
  )
]

# Matched tumor-normal
integrated[, Matched_support :=
  safe_prop(
    Matched_significant,
    Matched_cancers_tested
  )
]

# Combine primary and matched evidence
integrated[, Primary_RNA_support :=
  rowMeans(
    cbind(
      TN_support,
      Matched_support
    ),
    na.rm = TRUE
  )
]

# Stage-I
integrated[, Early_stage_support :=
  safe_prop(
    StageI_significant,
    StageI_cancers_tested
  )
]

# Stage progression
integrated[, Stage_progression_support :=
  safe_prop(
    Stage_trend_significant,
    Stage_cancers_tested
  )
]

# Grade progression
integrated[, Grade_progression_support :=
  safe_prop(
    Grade_trend_significant,
    Grade_cancers_tested
  )
]

integrated[, Progression_support :=
  rowMeans(
    cbind(
      Stage_progression_support,
      Grade_progression_support
    ),
    na.rm = TRUE
  )
]

# Survival
integrated[, Survival_support :=
  safe_prop(
    Cox_FDR_significant,
    Projects_analysed
  )
]

# Leukocyte-adjusted stage
integrated[, Adjusted_stage_support :=
  safe_prop(
    Adjusted_stage_trend_sig,
    Adjusted_stage_cancers
  )
]

# Leukocyte-adjusted grade
integrated[, Adjusted_grade_support :=
  safe_prop(
    Adjusted_grade_trend_sig,
    Adjusted_grade_cancers
  )
]

integrated[, Immune_adjusted_support :=
  rowMeans(
    cbind(
      Adjusted_stage_support,
      Adjusted_grade_support
    ),
    na.rm = TRUE
  )
]

# Leukocyte association:
# recorded separately, NOT penalized
integrated[, Leukocyte_association_fraction :=
  safe_prop(
    Significant_LF_correlations,
    Cancers_with_LF
  )
]

# ============================================================
# 10. UNWEIGHTED RNA EVIDENCE INDEX
#
# Five independent RNA-level domains:
#   primary tumour behaviour
#   early-stage behaviour
#   stage/grade progression
#   survival
#   immune-adjusted progression
#
# Not a biological probability.
# Used only as a transparent triage index.
# ============================================================

integrated[, RNA_evidence_index :=
  rowMeans(
    cbind(
      Primary_RNA_support,
      Early_stage_support,
      Progression_support,
      Survival_support,
      Immune_adjusted_support
    ),
    na.rm = TRUE
  )
]

# ============================================================
# 11. JAK3 BENCHMARK COMPARISON
#
# For each domain:
# candidate gets 1 if support >= JAK3.
#
# This answers:
# "How many RNA evidence domains are at least as strong
#  as the completed JAK3 reference?"
# ============================================================

jak3 <- integrated[
  gene == "JAK3"
]

stopifnot(
  nrow(jak3) == 1
)

integrated[, At_least_JAK3_primary :=
  Primary_RNA_support >=
  jak3$Primary_RNA_support
]

integrated[, At_least_JAK3_early :=
  Early_stage_support >=
  jak3$Early_stage_support
]

integrated[, At_least_JAK3_progression :=
  Progression_support >=
  jak3$Progression_support
]

integrated[, At_least_JAK3_survival :=
  Survival_support >=
  jak3$Survival_support
]

integrated[, At_least_JAK3_immune_adjusted :=
  Immune_adjusted_support >=
  jak3$Immune_adjusted_support
]

integrated[, JAK3_benchmark_domains :=
  rowSums(
    cbind(
      At_least_JAK3_primary,
      At_least_JAK3_early,
      At_least_JAK3_progression,
      At_least_JAK3_survival,
      At_least_JAK3_immune_adjusted
    ),
    na.rm = TRUE
  )
]

# ============================================================
# 12. ORIGINAL MEMORY CLASS
# ============================================================

integrated[, Memory_class :=
  fifelse(
    reference_gene %in% TRUE,
    "JAK3 reference",

    fifelse(
      strict_persistent %in% TRUE,
      "Strict persistent",

      fifelse(
        selection_group ==
        "Washout_emergent",
        "Washout emergent",
        "Other"
      )
    )
  )
]

# Keep FJX1 scope explicitly visible
integrated[, Scope_note :=
  fifelse(
    gene == "FJX1",
    paste0(
      "Retained from original candidate panel; ",
      "confirm kinase-family scope before final mechanistic ranking"
    ),
    ""
  )
]

# ============================================================
# 13. RNA-BASED ATAC TRIAGE
#
# IMPORTANT:
# This is NOT the final biological ranking.
# It simply identifies candidates with broad RNA support
# while retaining original memory status.
# ============================================================

integrated[, RNA_ATAC_priority :=
  fifelse(
    gene == "JAK3",
    "REFERENCE",

    fifelse(
      JAK3_benchmark_domains >= 4 &
      Memory_class %in%
      c(
        "Strict persistent",
        "Washout emergent"
      ),
      "HIGH",

      fifelse(
        JAK3_benchmark_domains >= 2 &
        Memory_class %in%
        c(
          "Strict persistent",
          "Washout emergent"
        ),
        "INTERMEDIATE",
        "EXPLORATORY"
      )
    )
  )
]

# ============================================================
# 14. ORDER FOR REVIEW
# ============================================================

priority_order <- c(
  "REFERENCE",
  "HIGH",
  "INTERMEDIATE",
  "EXPLORATORY"
)

integrated[
  ,
  Priority_order :=
    match(
      RNA_ATAC_priority,
      priority_order
    )
]

setorder(
  integrated,
  Priority_order,
  -JAK3_benchmark_domains,
  -RNA_evidence_index,
  gene
)

# ============================================================
# 15. SAVE COMPLETE TABLE
# ============================================================

fwrite(
  integrated,
  file.path(
    outdir,
    "09_integrated_multikinase_evidence.csv"
  )
)

# ============================================================
# 16. CREATE COMPACT SUPERVISOR TABLE
# ============================================================

compact <- integrated[
  ,
  .(
    gene,

    Memory_class,

    Day6_log2FC =
      log2FC_Day6,

    Day6_FDR =
      FDR_Day6,

    TN_sig =
      TN_significant,

    Matched_sig =
      Matched_significant,

    StageI_sig =
      StageI_significant,

    Stage_trend_sig =
      Stage_trend_significant,

    Grade_trend_sig =
      Grade_trend_significant,

    Survival_sig =
      Cox_FDR_significant,

    LF_correlations_sig =
      Significant_LF_correlations,

    Adjusted_stage_sig =
      Adjusted_stage_trend_sig,

    Adjusted_grade_sig =
      Adjusted_grade_trend_sig,

    RNA_evidence_index =
      round(
        RNA_evidence_index,
        4
      ),

    JAK3_benchmark_domains,

    RNA_ATAC_priority,

    Scope_note
  )
]

fwrite(
  compact,
  file.path(
    outdir,
    "09_supervisor_multikinase_comparison.csv"
  )
)

# ============================================================
# 17. JAK3 BENCHMARK VALUES
# ============================================================

benchmark <- jak3[
  ,
  .(
    Primary_RNA_support,
    Early_stage_support,
    Progression_support,
    Survival_support,
    Immune_adjusted_support,
    Leukocyte_association_fraction,
    RNA_evidence_index
  )
]

fwrite(
  benchmark,
  file.path(
    outdir,
    "09_JAK3_benchmark_values.csv"
  )
)

# ============================================================
# 18. CONSOLE REPORT
# ============================================================

cat("\n===== JAK3 BENCHMARK =====\n")

print(
  benchmark
)

cat("\n===== INTEGRATED MULTIKINASE RANKING =====\n")

print(
  compact[
    ,
    .(
      gene,
      Memory_class,
      TN_sig,
      Matched_sig,
      StageI_sig,
      Stage_trend_sig,
      Grade_trend_sig,
      Survival_sig,
      LF_correlations_sig,
      Adjusted_stage_sig,
      Adjusted_grade_sig,
      RNA_evidence_index,
      JAK3_benchmark_domains,
      RNA_ATAC_priority
    )
  ]
)

cat("\n===== ATAC TRIAGE COUNTS =====\n")

print(
  compact[
    ,
    .N,
    by = RNA_ATAC_priority
  ][order(
    match(
      RNA_ATAC_priority,
      priority_order
    )
  )]
)

cat("\n09 COMPLETE\n")

