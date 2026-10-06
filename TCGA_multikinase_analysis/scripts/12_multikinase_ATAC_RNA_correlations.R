
library(data.table)
setDTthreads(1)

root <- "TCGA_multikinase_analysis"

outdir <- file.path(
  root,
  "results/ATAC"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n========================================\n")
cat("MULTIKINASE PATIENT-MATCHED ATAC-RNA\n")
cat("========================================\n")

# ============================================================
# 1. Load patient-level ATAC
# ============================================================

cat("\n1. Loading patient-level ATAC values...\n")

atac <- fread(
  file.path(
    outdir,
    "11_multikinase_linked_ATAC_patient_values.csv"
  ),
  nThread = 1
)

cat(
  "ATAC rows:",
  nrow(atac),
  "\n"
)

cat(
  "ATAC genes:",
  uniqueN(atac$gene),
  "\n"
)

cat(
  "ATAC enhancer units:",
  uniqueN(atac$Enhancer_ID),
  "\n"
)

cat(
  "ATAC patients:",
  uniqueN(atac$patient_id),
  "\n"
)

# ============================================================
# 2. Load RNA
# ============================================================

cat("\n2. Loading multikinase RNA expression...\n")

rna <- fread(
  file.path(
    root,
    "results/expression",
    "03_multikinase_expression_clinical_complete.csv"
  ),
  nThread = 1
)

rna[, Cancer :=
  sub(
    "^TCGA-",
    "",
    project
  )
]

# canonical primary malignant
rna[, primary_malignant :=
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

rna <- rna[
  primary_malignant == TRUE
]

# ============================================================
# 3. Patient-level RNA
#
# Same handling as main pipeline:
# mean TPM first, then log2(TPM + 1)
# ============================================================

rna_patient <- rna[
  ,
  .(
    TPM =
      mean(
        TPM,
        na.rm = TRUE
      )
  ),
  by = .(
    gene,
    Cancer,
    patient_id
  )
]

rna_patient[, log2TPM1 :=
  log2(
    TPM + 1
  )
]

cat(
  "RNA patient-level rows:",
  nrow(rna_patient),
  "\n"
)

# ============================================================
# 4. ATAC-RNA match
# ============================================================

cat("\n3. Matching ATAC and RNA by gene + cancer + patient...\n")

matched <- merge(
  atac,
  rna_patient[
    ,
    .(
      gene,
      Cancer,
      patient_id,
      TPM,
      log2TPM1
    )
  ],
  by = c(
    "gene",
    "Cancer",
    "patient_id"
  ),
  all = FALSE
)

cat(
  "Matched patient-enhancer observations:",
  nrow(matched),
  "\n"
)

cat(
  "Matched unique patients:",
  uniqueN(matched$patient_id),
  "\n"
)

cat(
  "Matched genes:",
  uniqueN(matched$gene),
  "\n"
)

cat(
  "Matched cancers:",
  uniqueN(matched$Cancer),
  "\n"
)

fwrite(
  matched,
  file.path(
    outdir,
    "12_multikinase_ATAC_RNA_matched_patient_values.csv"
  )
)

# ============================================================
# 5. Match audit by gene / cancer
# ============================================================

atac_patients <- atac[
  ,
  .(
    ATAC_patients =
      uniqueN(patient_id)
  ),
  by = .(
    gene,
    Cancer
  )
]

rna_patients <- rna_patient[
  ,
  .(
    RNA_primary_patients =
      uniqueN(patient_id)
  ),
  by = .(
    gene,
    Cancer
  )
]

matched_patients <- matched[
  ,
  .(
    Matched_ATAC_RNA_patients =
      uniqueN(patient_id)
  ),
  by = .(
    gene,
    Cancer
  )
]

match_audit <- merge(
  atac_patients,
  rna_patients,
  by = c(
    "gene",
    "Cancer"
  ),
  all = TRUE
)

match_audit <- merge(
  match_audit,
  matched_patients,
  by = c(
    "gene",
    "Cancer"
  ),
  all = TRUE
)

match_audit[
  is.na(Matched_ATAC_RNA_patients),
  Matched_ATAC_RNA_patients := 0L
]

match_audit[, Match_percent :=
  100 *
  Matched_ATAC_RNA_patients /
  ATAC_patients
]

fwrite(
  match_audit,
  file.path(
    outdir,
    "12_multikinase_ATAC_RNA_match_audit.csv"
  )
)

# ============================================================
# 6. Safe Spearman
# ============================================================

safe_spearman <- function(x, y) {

  keep <-
    is.finite(x) &
    is.finite(y)

  x <- x[keep]
  y <- y[keep]

  n <- length(x)

  if (
    n < 5 ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {

    return(
      c(
        N = n,
        rho = NA_real_,
        P = NA_real_
      )
    )
  }

  out <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )

  c(
    N =
      n,

    rho =
      unname(
        out$estimate
      ),

    P =
      out$p.value
  )
}

# ============================================================
# 7. All gene × cancer × enhancer combinations
# ============================================================

all_combos <- unique(
  atac[
    ,
    .(
      gene,
      Cancer,
      Enhancer_ID,
      Enhancer_Distance,
      Enhancer_Correlation,
      Enhancer_FDR,
      Constituent_peak_count
    )
  ]
)

results_list <- vector(
  "list",
  nrow(all_combos)
)

cat(
  "\n4. Testing",
  nrow(all_combos),
  "gene-cancer-enhancer combinations...\n"
)

for (i in seq_len(nrow(all_combos))) {

  g <-
    all_combos$gene[i]

  ca <-
    all_combos$Cancer[i]

  enh <-
    all_combos$Enhancer_ID[i]

  d <- matched[
    gene == g &
    Cancer == ca &
    Enhancer_ID == enh
  ]

  stat <- safe_spearman(
    d$accessibility,
    d$log2TPM1
  )

  results_list[[i]] <-
    data.table(

      gene =
        g,

      Cancer =
        ca,

      Enhancer_ID =
        enh,

      Matched_N =
        as.integer(
          stat["N"]
        ),

      Spearman_rho =
        as.numeric(
          stat["rho"]
        ),

      P =
        as.numeric(
          stat["P"]
        ),

      DataS7_panCancer_correlation =
        all_combos$Enhancer_Correlation[i],

      DataS7_panCancer_FDR =
        all_combos$Enhancer_FDR[i],

      Enhancer_Distance =
        all_combos$Enhancer_Distance[i],

      Constituent_peak_count =
        all_combos$Constituent_peak_count[i]
    )

  if (i %% 250 == 0) {

    cat(
      "Processed:",
      i,
      "/",
      nrow(all_combos),
      "\n"
    )
  }
}

res <- rbindlist(
  results_list,
  fill = TRUE
)

# ============================================================
# 8. Sample-size classification
# ============================================================

res[, Analysis_status :=
  fifelse(
    Matched_N >= 10,
    "Primary",

    fifelse(
      Matched_N >= 5,
      "Exploratory_low_N",
      "Not_tested_N_lt_5"
    )
  )
]

# ============================================================
# 9. Gene-specific global BH-FDR
#
# Equivalent concept to original JAK3 global FDR,
# but applied separately to each gene.
# ============================================================

res[, Gene_global_FDR := NA_real_]

res[
  Analysis_status == "Primary" &
  !is.na(P),
  Gene_global_FDR :=
    p.adjust(
      P,
      method = "BH"
    ),
  by = gene
]

# ============================================================
# 10. Within gene + cancer BH-FDR
# ============================================================

res[, Within_gene_cancer_FDR := NA_real_]

res[
  Analysis_status == "Primary" &
  !is.na(P),
  Within_gene_cancer_FDR :=
    p.adjust(
      P,
      method = "BH"
    ),
  by = .(
    gene,
    Cancer
  )
]

# ============================================================
# 11. Multikinase-wide global FDR
# ============================================================

res[, Multikinase_Global_FDR := NA_real_]

idx <- which(
  res$Analysis_status == "Primary" &
  !is.na(res$P)
)

res$Multikinase_Global_FDR[idx] <-
  p.adjust(
    res$P[idx],
    method = "BH"
  )

# ============================================================
# 12. Direction concordance with Data S7
# ============================================================

res[, Direction_concordant_with_DataS7 :=
  fifelse(
    is.na(Spearman_rho) |
    is.na(DataS7_panCancer_correlation),

    NA,

    sign(Spearman_rho) ==
    sign(DataS7_panCancer_correlation)
  )
]

res[, Gene_global_significant :=
  !is.na(Gene_global_FDR) &
  Gene_global_FDR < 0.05
]

res[, Multikinase_global_significant :=
  !is.na(Multikinase_Global_FDR) &
  Multikinase_Global_FDR < 0.05
]

# Temporary absolute-correlation column for sorting
res[, Abs_Spearman_rho := abs(Spearman_rho)]

setorder(
  res,
  gene,
  Gene_global_FDR,
  -Abs_Spearman_rho,
  Cancer,
  Enhancer_ID
)

res[, Abs_Spearman_rho := NULL]

fwrite(
  res,
  file.path(
    outdir,
    "12_multikinase_ATAC_RNA_correlations.csv"
  )
)

# ============================================================
# 13. Gene-level summary
# ============================================================

gene_summary <- res[
  ,
  .(
    Official_elements =
      uniqueN(
        Enhancer_ID
      ),

    Cancer_enhancer_tests =
      .N,

    Primary_tests =
      sum(
        Analysis_status == "Primary"
      ),

    Exploratory_tests =
      sum(
        Analysis_status ==
        "Exploratory_low_N"
      ),

    Primary_gene_FDR_sig =
      sum(
        Gene_global_FDR < 0.05,
        na.rm = TRUE
      ),

    Primary_multikinase_FDR_sig =
      sum(
        Multikinase_Global_FDR < 0.05,
        na.rm = TRUE
      ),

    Positive_gene_FDR_sig =
      sum(
        Gene_global_FDR < 0.05 &
        Spearman_rho > 0,
        na.rm = TRUE
      ),

    Negative_gene_FDR_sig =
      sum(
        Gene_global_FDR < 0.05 &
        Spearman_rho < 0,
        na.rm = TRUE
      ),

    DataS7_direction_concordant_sig =
      sum(
        Gene_global_FDR < 0.05 &
        Direction_concordant_with_DataS7 %in% TRUE,
        na.rm = TRUE
      ),

    Significant_cancers =
      uniqueN(
        Cancer[
          Gene_global_FDR < 0.05
        ]
      )
  ),
  by = gene
]

gene_summary[, Reference :=
  gene == "JAK3"
]

setorder(
  gene_summary,
  -Primary_gene_FDR_sig,
  -Significant_cancers,
  gene
)

fwrite(
  gene_summary,
  file.path(
    outdir,
    "12_ATAC_RNA_gene_comparison_summary.csv"
  )
)

# ============================================================
# 14. Historical nine-cancer JAK3 subset
#
# Directly comparable to completed JAK3 analysis.
# ============================================================

jak3_priority_cancers <- c(
  "CHOL",
  "KIRC",
  "KIRP",
  "THCA",
  "HNSC",
  "STAD",
  "LUAD",
  "COAD",
  "LIHC"
)

jak3_validation <- res[
  gene == "JAK3" &
  Cancer %in% jak3_priority_cancers
]

# Recalculate historical JAK3-style global FDR
# only across primary tests in these nine cancers.
jak3_validation[
  ,
  Historical_9cancer_Global_FDR :=
    NA_real_
]

jidx <- which(
  jak3_validation$Analysis_status ==
  "Primary" &
  !is.na(jak3_validation$P)
)

jak3_validation$Historical_9cancer_Global_FDR[jidx] <-
  p.adjust(
    jak3_validation$P[jidx],
    method = "BH"
  )

fwrite(
  jak3_validation,
  file.path(
    outdir,
    "12_JAK3_historical_9cancer_validation.csv"
  )
)

# ============================================================
# 15. Console output
# ============================================================

cat("\n========================================\n")
cat("ATAC-RNA GENE-LEVEL COMPARISON\n")
cat("========================================\n")

print(
  gene_summary
)

cat("\n===== JAK3 HISTORICAL NINE-CANCER CHECK =====\n")

print(
  jak3_validation[
    ,
    .(
      Cancer,
      Enhancer_ID,
      Matched_N,
      Spearman_rho,
      P,
      Historical_9cancer_Global_FDR,
      Analysis_status
    )
  ][
    order(
      Historical_9cancer_Global_FDR
    )
  ]
)

cat("\n===== TOP MULTIKINASE ATAC-RNA ASSOCIATIONS =====\n")

print(
  head(
    res[
      Analysis_status == "Primary"
    ][
      order(
        Multikinase_Global_FDR
      )
    ],
    30
  )[
    ,
    .(
      gene,
      Cancer,
      Enhancer_ID,
      Matched_N,
      Spearman_rho,
      Gene_global_FDR,
      Multikinase_Global_FDR,
      Direction_concordant_with_DataS7
    )
  ]
)

cat("\n12 COMPLETE\n")

