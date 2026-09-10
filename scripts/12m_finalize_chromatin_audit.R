# ============================================================
# 12m — Finalize corrected Day-6 chromatin audit
#
# PURPOSE
#   Integrate:
#     1. locally reproduced 12j limma results
#     2. independently reported corrected Tobias benchmarks
#
# IMPORTANT
#   Tobias benchmark FDRs are NOT claimed as locally reproduced.
#   EPHB2 and MET benchmark emails did not explicitly specify
#   which Day-6 contrast carried the reported FDR, so that
#   information is kept separate from local contrast rows.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

outdir <- "results/12_kinase_reevaluation"
figdir <- "figures/12_kinase_reevaluation"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

local_file <- file.path(outdir, "12j_target_results.csv")

if (!file.exists(local_file)) {
  stop("Missing 12j_target_results.csv")
}

local <- read.csv(
  local_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

targets <- c(
  "JAK3",
  "EPHB2",
  "MET",
  "MAP3K8",
  "BMPR1A"
)

# ------------------------------------------------------------
# 1. Local reproduction: formal LPS vs RPMI promoter results
# ------------------------------------------------------------

local_lps <- local[
  local$SYMBOL %in% targets &
    local$Contrast == "LPS_vs_RPMI" &
    local$Mark %in% c("H3K27ac", "H3K4me3"),
  c(
    "SYMBOL",
    "Mark",
    "Contrast",
    "Log_signal_difference",
    "P_value",
    "FDR",
    "Direction"
  ),
  drop = FALSE
]

local_lps$Local_significant_FDR05 <- local_lps$FDR < 0.05

write.csv(
  local_lps,
  file.path(
    outdir,
    "12m_local_formal_promoter_results.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 2. Independent corrected Tobias benchmark
#
# These values come from supervisor correspondence.
# They are NOT outputs of the local model.
#
# Do not assign EPHB2/MET values to both BG and LPS contrasts:
# the correspondence quoted the corrected Day-6 promoter FDR
# but did not explicitly identify the contrast in those lines.
# ------------------------------------------------------------

benchmark <- data.frame(
  SYMBOL = c(
    "JAK3",
    "JAK3",
    "EPHB2",
    "MET",
    "MAP3K8",
    "BMPR1A"
  ),

  Mark = c(
    "H3K27ac",
    "H3K4me3",
    "H3K27ac",
    "H3K27ac",
    "multiple",
    "multiple"
  ),

  Reported_FDR = c(
    0.0065,
    0.0141,
    0.030,
    0.051,
    NA,
    NA
  ),

  Reported_contrast = c(
    "LPS_vs_RPMI",
    "LPS_vs_RPMI",
    "Day6 promoter comparison; exact contrast not specified in quoted email",
    "Day6 promoter comparison; exact contrast not specified in quoted email",
    "Day6 promoter tests",
    "Day6 promoter tests"
  ),

  Independent_interpretation = c(
    "Significant corrected promoter support",
    "Significant corrected promoter support",
    "Significant corrected H3K27ac support",
    "Not significant after normalization correction",
    "No significant corrected promoter support",
    "No significant corrected promoter support"
  ),

  stringsAsFactors = FALSE
)

write.csv(
  benchmark,
  file.path(
    outdir,
    "12m_independent_corrected_benchmark.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Candidate-level final audit status
# ------------------------------------------------------------

candidate_status <- data.frame(
  SYMBOL = targets,

  RNA_status = c(
    "Persistent RNA-supported candidate",
    "Washout-emergent RNA candidate",
    "RNA-supported candidate",
    "RNA-supported candidate",
    "RNA-supported candidate"
  ),

  Local_chromatin_reproduction = c(
    "Strong nominal H3K27ac/H3K4me3 signal; exact Tobias FDRs not reproduced",
    "Independent significant H3K27ac result not reproduced locally",
    "Strong nominal H3K27ac signal; corrected independent FDR is non-significant",
    "No local formal promoter support",
    "No local formal promoter support"
  ),

  Independent_corrected_chromatin = c(
    "H3K27ac FDR 0.0065; H3K4me3 FDR 0.0141",
    "H3K27ac FDR 0.030",
    "H3K27ac FDR 0.051 after correction",
    "No significant corrected promoter result",
    "No significant corrected promoter result"
  ),

  Downstream_status = c(
    "PRIORITIZE",
    "PRIORITIZE",
    "RNA_ONLY",
    "RNA_ONLY",
    "RNA_ONLY"
  ),

  Priority_basis = c(
    "Our RNA evidence + independent corrected promoter chromatin support",
    "Our washout RNA evidence + independent corrected promoter H3K27ac support",
    "RNA evidence; corrected promoter chromatin does not cross FDR 0.05",
    "RNA evidence without confirmed formal promoter chromatin support",
    "RNA evidence without confirmed formal promoter chromatin support"
  ),

  stringsAsFactors = FALSE
)

write.csv(
  candidate_status,
  file.path(
    outdir,
    "12m_FINAL_candidate_status.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. Reproducibility summary
# ------------------------------------------------------------

repro <- data.frame(
  Component = c(
    "Promoter definition",
    "Signal quantification",
    "Transformation",
    "Design",
    "Testing scope",
    "RPMI_d6_rep2 normalization",
    "Exact Tobias FDR reproduction",
    "Downstream candidate decision"
  ),

  Status = c(
    "TSS +/- 2 kb",
    "Mean BigWig promoter signal",
    "log2(signal + 1) local implementation",
    "Unpaired",
    "KinHub kinase genes",
    "Three NotNormalized tracks rescaled locally using normalized-peer global-mean matching",
    "NOT reproduced exactly",
    "JAK3 + EPHB2"
  ),

  Notes = c(
    "",
    "",
    "Tobias confirmed log transformation but did not specify exact pseudocount",
    "",
    "Tobias confirmed testing was restricted to kinase genes",
    "Local scaling is an approximation of the correction; exact Tobias scaling factors were not supplied",
    "Independent corrected results are retained separately from local statistics",
    "Priority incorporates our RNA results and independent corrected chromatin review"
  ),

  stringsAsFactors = FALSE
)

write.csv(
  repro,
  file.path(
    outdir,
    "12m_reproducibility_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Simple final figure
# ------------------------------------------------------------

png(
  file.path(
    figdir,
    "12m_FINAL_candidate_summary.png"
  ),
  width = 1800,
  height = 1050,
  res = 180
)

par(mar = c(7, 8, 4, 2))

genes <- rev(candidate_status$SYMBOL)

status_value <- ifelse(
  rev(candidate_status$Downstream_status) == "PRIORITIZE",
  2,
  1
)

barplot(
  status_value,
  names.arg = genes,
  horiz = TRUE,
  las = 1,
  xlim = c(0, 2.6),
  axes = FALSE,
  xlab = "",
  main = "Corrected candidate integration"
)

axis(
  1,
  at = c(1, 2),
  labels = c(
    "RNA-supported",
    "Downstream priority"
  )
)

mtext(
  "Priority integrates local RNA evidence with independently corrected promoter review",
  side = 1,
  line = 5,
  cex = 0.8
)

dev.off()

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("12m FINAL CHROMATIN AUDIT\n")
cat("========================================\n\n")

print(candidate_status, row.names = FALSE)

cat("\nIMPORTANT:\n")
cat(
  "Tobias FDR values are independent corrected benchmarks ",
  "and are NOT claimed as locally reproduced statistics.\n",
  sep = ""
)

cat("\nDownstream priority:\n")
cat(
  paste(
    candidate_status$SYMBOL[
      candidate_status$Downstream_status == "PRIORITIZE"
    ],
    collapse = ", "
  ),
  "\n"
)

cat("\nFiles written:\n")
cat("  12m_local_formal_promoter_results.csv\n")
cat("  12m_independent_corrected_benchmark.csv\n")
cat("  12m_FINAL_candidate_status.csv\n")
cat("  12m_reproducibility_summary.csv\n")
cat("  figures/12_kinase_reevaluation/12m_FINAL_candidate_summary.png\n")

cat("\n========================================\n")
cat("12m COMPLETE\n")
cat("========================================\n")
