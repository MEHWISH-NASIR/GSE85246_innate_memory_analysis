
library(data.table)

# Avoid data.table multithreading instability on this system
setDTthreads(1)

root <- "TCGA_multikinase_analysis"
outdir <- file.path(root, "results/primary")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("1. Loading merged multikinase dataset...\n")

d <- fread(
  file.path(
    root,
    "results/expression",
    "03_multikinase_expression_clinical_complete.csv"
  ),
  nThread = 1
)

cat("2. Dataset loaded:", nrow(d), "rows\n")

# ============================================================
# Analysis class
# ============================================================

d[, analysis_class := fifelse(
  project == "TCGA-LAML" &
    sample_type == "Primary Blood Derived Cancer - Peripheral Blood",
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

cat("3. Analysis classes assigned\n")

# ============================================================
# Patient-level aggregation
# Same approach as original JAK3 workflow
# ============================================================

p <- d[
  analysis_class %in% c("Primary malignant", "Normal"),
  .(
    TPM = mean(TPM),
    expression_samples = .N
  ),
  by = .(
    gene,
    project,
    patient_id,
    analysis_class
  )
]

p[, log2TPM1 := log2(TPM + 1)]

stopifnot(
  !anyDuplicated(
    p[, .(gene, project, patient_id, analysis_class)]
  ),
  !anyNA(p$TPM)
)

fwrite(
  p,
  file.path(outdir, "04_patient_level_expression.csv")
)

cat("4. Patient-level aggregation complete:", nrow(p), "rows\n")

tumor <- p[analysis_class == "Primary malignant"]
normal <- p[analysis_class == "Normal"]

# ============================================================
# Safe Wilcoxon function
# ============================================================

safe_wilcox <- function(x, y, paired = FALSE) {

  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }

  tryCatch(
    wilcox.test(
      x,
      y,
      paired = paired,
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

# ============================================================
# Tumor vs normal
# ============================================================

groups <- unique(
  tumor[, .(gene, project)]
)

tn_list <- vector("list", nrow(groups))

for (i in seq_len(nrow(groups))) {

  g <- groups$gene[i]
  cancer <- groups$project[i]

  t <- tumor[
    gene == g &
      project == cancer
  ]

  n <- normal[
    gene == g &
      project == cancer
  ]

  if (nrow(n) == 0) {
    next
  }

  delta <-
    median(t$log2TPM1) -
    median(n$log2TPM1)

  tn_list[[i]] <- data.table(
    gene = g,
    project = cancer,

    tumor_n = nrow(t),
    normal_n = nrow(n),

    tumor_median_TPM =
      median(t$TPM),

    normal_median_TPM =
      median(n$TPM),

    Median_log2TPM1_difference =
      delta,

    Direction =
      fifelse(
        delta > 0,
        "UP",
        fifelse(
          delta < 0,
          "DOWN",
          "UNCHANGED"
        )
      ),

    Wilcoxon_P =
      safe_wilcox(
        t$log2TPM1,
        n$log2TPM1
      ),

    Small_normal_n =
      nrow(n) < 10
  )

  if (i %% 25 == 0) {
    cat(
      "Tumor-normal groups processed:",
      i, "/", nrow(groups), "\n"
    )
  }
}

tn <- rbindlist(
  tn_list,
  fill = TRUE
)

tn[, Wilcoxon_FDR :=
  p.adjust(
    Wilcoxon_P,
    method = "BH"
  ),
  by = gene
]

tn[, Global_FDR :=
  p.adjust(
    Wilcoxon_P,
    method = "BH"
  )
]

setorder(
  tn,
  gene,
  project
)

fwrite(
  tn,
  file.path(
    outdir,
    "04_multikinase_tumor_normal.csv"
  )
)

cat("5. Tumor-normal analysis complete:", nrow(tn), "comparisons\n")

# ============================================================
# Matched tumor-normal
# ============================================================

matched <- merge(
  tumor[, .(
    gene,
    project,
    patient_id,
    tumor_TPM = TPM,
    tumor_log2TPM1 = log2TPM1
  )],

  normal[, .(
    gene,
    project,
    patient_id,
    normal_TPM = TPM,
    normal_log2TPM1 = log2TPM1
  )],

  by = c(
    "gene",
    "project",
    "patient_id"
  )
)

matched[, log2_difference :=
  tumor_log2TPM1 -
  normal_log2TPM1
]

matched[, tumor_higher :=
  tumor_log2TPM1 >
  normal_log2TPM1
]

fwrite(
  matched,
  file.path(
    outdir,
    "04_multikinase_matched_patient_values.csv"
  )
)

mp <- matched[, .(
  matched_pairs = .N,

  median_log2_difference =
    median(log2_difference),

  tumor_higher_n =
    sum(tumor_higher),

  tumor_lower_or_equal_n =
    sum(!tumor_higher),

  Paired_Wilcoxon_P =
    safe_wilcox(
      tumor_log2TPM1,
      normal_log2TPM1,
      paired = TRUE
    ),

  Small_matched_n =
    .N < 10

), by = .(
  gene,
  project
)]

mp[, Paired_Wilcoxon_FDR :=
  p.adjust(
    Paired_Wilcoxon_P,
    method = "BH"
  ),
  by = gene
]

mp[, Global_FDR :=
  p.adjust(
    Paired_Wilcoxon_P,
    method = "BH"
  )
]

setorder(
  mp,
  gene,
  project
)

fwrite(
  mp,
  file.path(
    outdir,
    "04_multikinase_matched_tumor_normal.csv"
  )
)

cat("6. Matched tumor-normal analysis complete:", nrow(mp), "comparisons\n")

# ============================================================
# Final console summary
# ============================================================

cat("\n====================================\n")
cat("MULTIKINASE PRIMARY RNA ANALYSIS\n")
cat("====================================\n")

cat("Genes:", uniqueN(p$gene), "\n")
cat("Cancer projects:", uniqueN(p$project), "\n")
cat("Tumor-normal comparisons:", nrow(tn), "\n")
cat("Matched comparisons:", nrow(mp), "\n")

cat("\n===== JAK3 VALIDATION PREVIEW =====\n")

print(
  tn[
    gene == "JAK3",
    .(
      project,
      tumor_n,
      normal_n,
      Median_log2TPM1_difference,
      Wilcoxon_P,
      Wilcoxon_FDR
    )
  ]
)

cat("\n04 COMPLETE\n")

