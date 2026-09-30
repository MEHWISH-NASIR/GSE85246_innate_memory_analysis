# ============================================================
# 60b — Patient-matched JAK3 ATAC-RNA correlations
#
# Primary inference:
#   cancer × enhancer combinations with N >= 10 matched patients
#
# Exploratory:
#   N = 5-9
#
# Not tested:
#   N < 5
#
# Spearman correlation is used within each cancer.
# BH-FDR is calculated:
#   1) globally across primary-eligible tests
#   2) within each cancer across the three JAK3 elements
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

atac_file <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results/JAK3_linked_ATAC_tissue_values.csv"

rna_file <-
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_expression_clinical_complete_33projects.csv"

outdir <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

priority_cancers <- c(
  "CHOL","KIRC","KIRP","THCA","HNSC",
  "STAD","LUAD","COAD","LIHC"
)

enhancer_order <- c(
  "JAK3_m1",
  "JAK3_p1",
  "JAK3_p2"
)

cat("\n========================================\n")
cat("60b JAK3 PATIENT-MATCHED ATAC-RNA\n")
cat("========================================\n\n")

# ------------------------------------------------------------
# Read ATAC tissue-level values
# ------------------------------------------------------------

atac <- read.csv(
  atac_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Collapse multiple tissue fragments from the same patient
# to one patient-level ATAC value per enhancer.

atac_patient <- aggregate(
  accessibility ~ Cancer + Enhancer_ID + patient_id,
  data = atac,
  FUN = mean,
  na.rm = TRUE
)

# Preserve official Data S7 metadata
link_meta <- unique(
  atac[
    ,
    c(
      "Enhancer_ID",
      "Enhancer_Distance",
      "Enhancer_Correlation",
      "Enhancer_FDR"
    )
  ]
)

# ------------------------------------------------------------
# Read RNA
# ------------------------------------------------------------

rna <- read.csv(
  rna_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rna$Cancer <- sub(
  "^TCGA-",
  "",
  rna$project
)

rna$patient_id <-
  as.character(
    rna$cases.submitter_id
  )

rna <- rna[
  rna$Cancer %in% priority_cancers &
  rna$sample_type == "Primary Tumor",
]

# Keep patient-level RNA handling consistent with the main
# JAK3 pipeline: mean TPM first, then log2(TPM + 1).

if ("TPM" %in% names(rna)) {

  rna_patient <- aggregate(
    TPM ~ Cancer + patient_id,
    data = rna,
    FUN = mean,
    na.rm = TRUE
  )

  rna_patient$JAK3_log2TPM1 <-
    log2(
      rna_patient$TPM + 1
    )

} else {

  rna_patient <- aggregate(
    log2TPM1 ~ Cancer + patient_id,
    data = rna,
    FUN = mean,
    na.rm = TRUE
  )

  names(rna_patient)[
    names(rna_patient) == "log2TPM1"
  ] <- "JAK3_log2TPM1"
}

# ------------------------------------------------------------
# Merge matched ATAC + RNA
# ------------------------------------------------------------

matched <- merge(
  atac_patient,
  rna_patient[
    ,
    c(
      "Cancer",
      "patient_id",
      "JAK3_log2TPM1"
    )
  ],
  by = c(
    "Cancer",
    "patient_id"
  ),
  all = FALSE
)

matched <- merge(
  matched,
  link_meta,
  by = "Enhancer_ID",
  all.x = TRUE
)

write.csv(
  matched,
  file.path(
    outdir,
    "JAK3_ATAC_RNA_matched_patient_values.csv"
  ),
  row.names = FALSE
)

cat(
  "Matched patient-enhancer observations:",
  nrow(matched),
  "\n"
)

cat(
  "Matched unique patients:",
  length(unique(matched$patient_id)),
  "\n\n"
)

# ------------------------------------------------------------
# Safe Spearman correlation
# ------------------------------------------------------------

safe_spearman <- function(x, y) {

  keep <- is.finite(x) & is.finite(y)

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

  z <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )

  c(
    N = n,
    rho = unname(z$estimate),
    P = z$p.value
  )
}

# ------------------------------------------------------------
# Cancer × enhancer correlations
# ------------------------------------------------------------

rows <- list()

for (ca in priority_cancers) {

  for (enh in enhancer_order) {

    d <- matched[
      matched$Cancer == ca &
      matched$Enhancer_ID == enh,
    ]

    stat <- safe_spearman(
      d$accessibility,
      d$JAK3_log2TPM1
    )

    meta <- link_meta[
      link_meta$Enhancer_ID == enh,
    ][1, ]

    rows[[
      length(rows) + 1
    ]] <- data.frame(

      Cancer = ca,
      Enhancer_ID = enh,

      Matched_N =
        as.integer(stat["N"]),

      Spearman_rho =
        as.numeric(stat["rho"]),

      P =
        as.numeric(stat["P"]),

      DataS7_panCancer_correlation =
        meta$Enhancer_Correlation,

      DataS7_panCancer_FDR =
        meta$Enhancer_FDR,

      stringsAsFactors = FALSE
    )
  }
}

res <- do.call(
  rbind,
  rows
)

# ------------------------------------------------------------
# Sample-size classification
# ------------------------------------------------------------

res$Analysis_status <- ifelse(
  res$Matched_N >= 10,
  "Primary",
  ifelse(
    res$Matched_N >= 5,
    "Exploratory_low_N",
    "Not_tested_N_lt_5"
  )
)

# ------------------------------------------------------------
# Global BH-FDR across PRIMARY tests only
# ------------------------------------------------------------

res$Global_FDR <- NA_real_

primary_idx <-
  res$Analysis_status == "Primary" &
  !is.na(res$P)

res$Global_FDR[primary_idx] <-
  p.adjust(
    res$P[primary_idx],
    method = "BH"
  )

# ------------------------------------------------------------
# Within-cancer BH-FDR
# ------------------------------------------------------------

res$Within_cancer_FDR <- NA_real_

for (ca in priority_cancers) {

  idx <-
    res$Cancer == ca &
    res$Analysis_status == "Primary" &
    !is.na(res$P)

  if (sum(idx) > 0) {

    res$Within_cancer_FDR[idx] <-
      p.adjust(
        res$P[idx],
        method = "BH"
      )
  }
}

# ------------------------------------------------------------
# Direction concordance with official Data S7 link
# ------------------------------------------------------------

res$Direction_concordant_with_DataS7 <-
  ifelse(
    is.na(res$Spearman_rho),
    NA,
    sign(res$Spearman_rho) ==
      sign(res$DataS7_panCancer_correlation)
  )

res$Primary_global_significant <-
  !is.na(res$Global_FDR) &
  res$Global_FDR < 0.05

# ------------------------------------------------------------
# Save full results
# ------------------------------------------------------------

write.csv(
  res,
  file.path(
    outdir,
    "JAK3_ATAC_RNA_within_cancer_correlations.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console summaries
# ------------------------------------------------------------

cat("========================================\n")
cat("PRIMARY TESTS — N >= 10\n")
cat("========================================\n\n")

primary <- res[
  res$Analysis_status == "Primary",
]

primary <- primary[
  order(
    primary$Global_FDR,
    -abs(primary$Spearman_rho),
    na.last = TRUE
  ),
]

print(
  primary[
    ,
    c(
      "Cancer",
      "Enhancer_ID",
      "Matched_N",
      "Spearman_rho",
      "P",
      "Global_FDR",
      "Within_cancer_FDR",
      "Direction_concordant_with_DataS7",
      "Primary_global_significant"
    )
  ],
  row.names = FALSE
)

cat("\n========================================\n")
cat("LOW-N / NOT-TESTED CANCERS\n")
cat("========================================\n\n")

print(
  res[
    res$Analysis_status != "Primary",
    c(
      "Cancer",
      "Enhancer_ID",
      "Matched_N",
      "Spearman_rho",
      "P",
      "Analysis_status"
    )
  ],
  row.names = FALSE
)

cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

cat(
  "Primary tests:",
  sum(res$Analysis_status == "Primary"),
  "\n"
)

cat(
  "Global FDR < 0.05:",
  sum(
    res$Primary_global_significant,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Positive direction among primary tests:",
  sum(
    res$Spearman_rho[
      res$Analysis_status == "Primary"
    ] > 0,
    na.rm = TRUE
  ),
  "of",
  sum(res$Analysis_status == "Primary"),
  "\n"
)

cat(
  "Direction concordant with Data S7:",
  sum(
    res$Direction_concordant_with_DataS7[
      res$Analysis_status == "Primary"
    ],
    na.rm = TRUE
  ),
  "of",
  sum(res$Analysis_status == "Primary"),
  "\n"
)

cat("\nSaved:\n")
cat(
  " ",
  file.path(
    outdir,
    "JAK3_ATAC_RNA_within_cancer_correlations.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("60b COMPLETE\n")
cat("========================================\n")
