
library(data.table)
setDTthreads(1)

root <- "TCGA_multikinase_analysis"

outdir <- file.path(
  root,
  "results/final_integration"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n========================================\n")
cat("FINAL MULTIKINASE RNA + ATAC INTEGRATION\n")
cat("========================================\n")

# ============================================================
# Helper
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
    num[ok] /
    den[ok]

  out
}

# ============================================================
# 1. Load integrated RNA evidence from Step 09
# ============================================================

cat("\n1. Loading integrated RNA evidence...\n")

rna <- fread(
  file.path(
    root,
    "results/integrated",
    "09_integrated_multikinase_evidence.csv"
  ),
  nThread = 1
)

cat(
  "RNA genes:",
  uniqueN(rna$gene),
  "\n"
)

# ============================================================
# 2. Load complete ATAC-RNA correlations
# ============================================================

cat("2. Loading ATAC-RNA results...\n")

atac <- fread(
  file.path(
    root,
    "results/ATAC",
    "12_multikinase_ATAC_RNA_correlations.csv"
  ),
  nThread = 1
)

cat(
  "ATAC-RNA genes:",
  uniqueN(atac$gene),
  "\n"
)

# ============================================================
# 3. Primary-eligible ATAC tests only
# ============================================================

primary <- atac[
  Analysis_status == "Primary" &
  !is.na(P)
]

cat(
  "Primary ATAC-RNA tests:",
  nrow(primary),
  "\n"
)

# ============================================================
# 4. Gene-level normalized ATAC evidence
# ============================================================

atac_summary <- primary[
  ,
  .(
    Official_elements =
      uniqueN(
        Enhancer_ID
      ),

    Primary_ATAC_tests =
      .N,

    Primary_ATAC_cancers =
      uniqueN(
        Cancer
      ),

    Gene_FDR_significant_tests =
      sum(
        Gene_global_FDR < 0.05,
        na.rm = TRUE
      ),

    Multikinase_FDR_significant_tests =
      sum(
        Multikinase_Global_FDR < 0.05,
        na.rm = TRUE
      ),

    Significant_ATAC_cancers =
      uniqueN(
        Cancer[
          Gene_global_FDR < 0.05
        ]
      ),

    Multikinase_significant_cancers =
      uniqueN(
        Cancer[
          Multikinase_Global_FDR < 0.05
        ]
      ),

    Positive_significant_tests =
      sum(
        Gene_global_FDR < 0.05 &
        Spearman_rho > 0,
        na.rm = TRUE
      ),

    Negative_significant_tests =
      sum(
        Gene_global_FDR < 0.05 &
        Spearman_rho < 0,
        na.rm = TRUE
      ),

    Direction_concordant_significant =
      sum(
        Gene_global_FDR < 0.05 &
        Direction_concordant_with_DataS7 %in% TRUE,
        na.rm = TRUE
      ),

    Median_abs_rho_significant =
      if (
        any(
          Gene_global_FDR < 0.05,
          na.rm = TRUE
        )
      ) {
        median(
          abs(
            Spearman_rho[
              Gene_global_FDR < 0.05
            ]
          ),
          na.rm = TRUE
        )
      } else {
        NA_real_
      },

    Max_abs_rho_significant =
      if (
        any(
          Gene_global_FDR < 0.05,
          na.rm = TRUE
        )
      ) {
        max(
          abs(
            Spearman_rho[
              Gene_global_FDR < 0.05
            ]
          ),
          na.rm = TRUE
        )
      } else {
        NA_real_
      }
  ),
  by = gene
]

# ============================================================
# 5. Normalize ATAC support
#
# Important:
# Counts alone would favor genes with many linked enhancers.
#
# Therefore calculate:
#   A) significant-test fraction
#   B) stringent multikinase-significant fraction
#   C) significant-cancer coverage
#
# These three are combined into one ATAC evidence domain.
# ============================================================

atac_summary[, ATAC_test_support :=
  safe_prop(
    Gene_FDR_significant_tests,
    Primary_ATAC_tests
  )
]

atac_summary[, ATAC_global_test_support :=
  safe_prop(
    Multikinase_FDR_significant_tests,
    Primary_ATAC_tests
  )
]

atac_summary[, ATAC_cancer_support :=
  safe_prop(
    Significant_ATAC_cancers,
    Primary_ATAC_cancers
  )
]

atac_summary[, ATAC_global_cancer_support :=
  safe_prop(
    Multikinase_significant_cancers,
    Primary_ATAC_cancers
  )
]

# Primary ATAC evidence domain:
#
# 1/3 gene-level significant-test fraction
# 1/3 stringent multikinase-wide significant-test fraction
# 1/3 cancer-level coverage
#
# ATAC remains ONE evidence domain, preventing it from
# overwhelming the five RNA domains.
atac_summary[, ATAC_evidence_support :=
  rowMeans(
    cbind(
      ATAC_test_support,
      ATAC_global_test_support,
      ATAC_cancer_support
    ),
    na.rm = TRUE
  )
]

# Directional consistency metric
atac_summary[, DataS7_concordance_fraction :=
  safe_prop(
    Direction_concordant_significant,
    Gene_FDR_significant_tests
  )
]

# ============================================================
# 6. Merge RNA + ATAC
# ============================================================

final <- merge(
  rna,
  atac_summary,
  by = "gene",
  all.x = TRUE
)

if (nrow(final) != 10) {
  stop(
    paste(
      "Expected 10 genes after integration; found",
      nrow(final)
    )
  )
}

# ============================================================
# 7. Six-domain final evidence index
#
# Original five RNA domains:
#   1 Primary RNA
#   2 Early-stage RNA
#   3 Stage/grade progression
#   4 Survival
#   5 Immune-adjusted progression
#
# New sixth domain:
#   6 ATAC-RNA regulatory coupling
#
# Equal domain weighting.
# ============================================================

final[, Final_multiomic_evidence_index :=
  rowMeans(
    cbind(
      Primary_RNA_support,
      Early_stage_support,
      Progression_support,
      Survival_support,
      Immune_adjusted_support,
      ATAC_evidence_support
    ),
    na.rm = TRUE
  )
]

# ============================================================
# 8. JAK3 benchmark
# ============================================================

jak3 <- final[
  gene == "JAK3"
]

if (nrow(jak3) != 1) {
  stop(
    "JAK3 reference row not uniquely identified."
  )
}

final[, At_least_JAK3_ATAC :=
  ATAC_evidence_support >=
  jak3$ATAC_evidence_support
]

# Existing five RNA benchmark domains from Step 09
# + new ATAC benchmark domain

final[, Final_JAK3_benchmark_domains :=
  JAK3_benchmark_domains +
  as.integer(
    At_least_JAK3_ATAC
  )
]

# Maximum = 6
# JAK3 itself should equal 6.

# ============================================================
# 9. Original innate-memory evidence retained separately
#
# We DO NOT blend Day-6 persistence into the TCGA score.
# It represents the discovery biology, not another TCGA test.
# ============================================================

final[, Discovery_support :=
  fifelse(
    gene == "JAK3",
    "Reference persistent kinase",

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

# ============================================================
# 10. Evidence tier
#
# Transparent triage relative to JAK3.
#
# Not a biological probability and not a formal statistical
# ranking.
# ============================================================

final[, Final_evidence_tier :=
  fifelse(
    gene == "JAK3",
    "REFERENCE",

    fifelse(
      Final_JAK3_benchmark_domains >= 5,
      "VERY_STRONG",

      fifelse(
        Final_JAK3_benchmark_domains >= 4,
        "STRONG",

        fifelse(
          Final_JAK3_benchmark_domains >= 2,
          "SUPPORTED",
          "EXPLORATORY"
        )
      )
    )
  )
]

# ============================================================
# 11. FJX1 scope warning
# ============================================================

final[, Final_scope_note :=
  fifelse(
    gene == "FJX1",

    paste0(
      "Strong computational candidate retained from original panel; ",
      "confirm kinase-family classification before presenting ",
      "as a kinase-level mechanistic candidate."
    ),

    ""
  )
]

# ============================================================
# 12. Ranking
#
# First by number of domains >= JAK3,
# then by final multiomic evidence index.
# ============================================================

tier_order <- c(
  "REFERENCE",
  "VERY_STRONG",
  "STRONG",
  "SUPPORTED",
  "EXPLORATORY"
)

final[, Tier_order :=
  match(
    Final_evidence_tier,
    tier_order
  )
]

setorder(
  final,
  Tier_order,
  -Final_JAK3_benchmark_domains,
  -Final_multiomic_evidence_index,
  gene
)

# ============================================================
# 13. Save complete integrated table
# ============================================================

fwrite(
  final,
  file.path(
    outdir,
    "13_final_multikinase_RNA_ATAC_evidence.csv"
  )
)

# ============================================================
# 14. Supervisor-ready compact table
# ============================================================

compact <- final[
  ,
  .(
    gene,

    Discovery_support,

    Day6_log2FC =
      log2FC_Day6,

    Day6_FDR =
      FDR_Day6,

    # Primary expression
    TN_sig =
      TN_significant,

    Matched_sig =
      Matched_significant,

    # Early disease
    StageI_sig =
      StageI_significant,

    # Progression
    Stage_trend_sig =
      Stage_trend_significant,

    Grade_trend_sig =
      Grade_trend_significant,

    # Clinical outcome
    Survival_sig =
      Cox_FDR_significant,

    # Immune biology
    LF_correlations =
      Significant_LF_correlations,

    Immune_adjusted_stage_sig =
      Adjusted_stage_trend_sig,

    Immune_adjusted_grade_sig =
      Adjusted_grade_trend_sig,

    # ATAC
    Official_ATAC_elements =
      Official_elements,

    Primary_ATAC_tests,

    ATAC_FDR_sig =
      Gene_FDR_significant_tests,

    ATAC_global_FDR_sig =
      Multikinase_FDR_significant_tests,

    ATAC_sig_cancers =
      Significant_ATAC_cancers,

    ATAC_test_support =
      round(
        ATAC_test_support,
        4
      ),

    ATAC_cancer_support =
      round(
        ATAC_cancer_support,
        4
      ),

    ATAC_evidence_support =
      round(
        ATAC_evidence_support,
        4
      ),

    ATAC_DataS7_concordance =
      round(
        DataS7_concordance_fraction,
        4
      ),

    Median_abs_ATAC_RNA_rho =
      round(
        Median_abs_rho_significant,
        4
      ),

    # Final integration
    RNA_evidence_index =
      round(
        RNA_evidence_index,
        4
      ),

    Final_multiomic_index =
      round(
        Final_multiomic_evidence_index,
        4
      ),

    JAK3_domains_RNA =
      JAK3_benchmark_domains,

    JAK3_domains_final =
      Final_JAK3_benchmark_domains,

    Final_evidence_tier,

    Final_scope_note
  )
]

fwrite(
  compact,
  file.path(
    outdir,
    "13_supervisor_final_multikinase_comparison.csv"
  )
)

# ============================================================
# 15. Domain-level table
#
# Useful later for heatmaps.
# ============================================================

domains <- final[
  ,
  .(
    gene,

    Primary_RNA_support,
    Early_stage_support,
    Progression_support,
    Survival_support,
    Immune_adjusted_support,
    ATAC_evidence_support,

    Leukocyte_association_fraction,

    Final_multiomic_evidence_index,

    Reference =
      gene == "JAK3"
  )
]

fwrite(
  domains,
  file.path(
    outdir,
    "13_final_domain_support_matrix.csv"
  )
)

# ============================================================
# 16. ATAC details table
# ============================================================

atac_details <- primary[
  Gene_global_FDR < 0.05,
  .(
    gene,
    Cancer,
    Enhancer_ID,
    Matched_N,
    Spearman_rho,
    P,
    Gene_global_FDR,
    Multikinase_Global_FDR,
    DataS7_panCancer_correlation,
    DataS7_panCancer_FDR,
    Direction_concordant_with_DataS7,
    Constituent_peak_count
  )
]

setorder(
  atac_details,
  gene,
  Gene_global_FDR,
  Cancer
)

fwrite(
  atac_details,
  file.path(
    outdir,
    "13_significant_ATAC_RNA_associations.csv"
  )
)

# ============================================================
# 17. JAK3 benchmark record
# ============================================================

jak3_benchmark <- final[
  gene == "JAK3",
  .(
    gene,

    Primary_RNA_support,
    Early_stage_support,
    Progression_support,
    Survival_support,
    Immune_adjusted_support,
    ATAC_evidence_support,

    RNA_evidence_index,
    Final_multiomic_evidence_index,

    Final_JAK3_benchmark_domains
  )
]

fwrite(
  jak3_benchmark,
  file.path(
    outdir,
    "13_JAK3_final_benchmark.csv"
  )
)

# ============================================================
# 18. Console output
# ============================================================

cat("\n===== JAK3 FINAL BENCHMARK =====\n")

print(
  jak3_benchmark
)

cat("\n===== FINAL MULTIKINASE COMPARISON =====\n")

print(
  compact[
    ,
    .(
      gene,
      Discovery_support,
      TN_sig,
      Matched_sig,
      StageI_sig,
      Stage_trend_sig,
      Grade_trend_sig,
      Survival_sig,
      Immune_adjusted_stage_sig,
      Immune_adjusted_grade_sig,
      Official_ATAC_elements,
      ATAC_FDR_sig,
      ATAC_sig_cancers,
      ATAC_evidence_support,
      RNA_evidence_index,
      Final_multiomic_index,
      JAK3_domains_final,
      Final_evidence_tier
    )
  ]
)

cat("\n===== ATAC NORMALIZED SUPPORT =====\n")

print(
  compact[
    ,
    .(
      gene,
      Official_ATAC_elements,
      Primary_ATAC_tests,
      ATAC_FDR_sig,
      ATAC_global_FDR_sig,
      ATAC_sig_cancers,
      ATAC_test_support,
      ATAC_cancer_support,
      ATAC_evidence_support,
      ATAC_DataS7_concordance,
      Median_abs_ATAC_RNA_rho
    )
  ]
)

cat("\n===== FINAL EVIDENCE TIERS =====\n")

print(
  compact[
    ,
    .(
      gene,
      JAK3_domains_final,
      Final_multiomic_index,
      Final_evidence_tier
    )
  ]
)

cat("\n13 COMPLETE\n")

