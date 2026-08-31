# ============================================================
# 06_make_final_kinase_figure.R
#
# GSE85246 / GSE85243
#
# FINAL INTEGRATED RNA KINASE FIGURES
#
# Input:
#   results/05_integrated_trajectory/
#   05_integrated_kinase_trajectory.csv
#
# No new statistical testing is performed here.
# ============================================================

rm(list = ls())

# ============================================================
# 1. INPUT / OUTPUT
# ============================================================

input_file <- paste0(
  "results/05_integrated_trajectory/",
  "05_integrated_kinase_trajectory.csv"
)

outdir <- "results/06_final_kinase_figure"
figdir <- "figures/06_final_kinase"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop(
    paste(
      "Missing Step-05 input:",
      input_file
    )
  )
}

# ============================================================
# 2. LOAD INTEGRATED DATA
# ============================================================

dat <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("06 — FINAL INTEGRATED KINASE FIGURES\n")
cat("========================================\n\n")

cat(
  "Integrated kinase rows:",
  nrow(dat),
  "\n"
)

# ============================================================
# 3. REQUIRED COLUMNS
# ============================================================

required_cols <- c(
  "ENSEMBL",
  "SYMBOL",
  "Day1_log2FC",
  "Day1_P",
  "Day6_memory_log2FC",
  "Day6_memory_FDR_kinase",
  "Response_pattern",
  "Interaction_P",
  "Interaction_FDR_genome",
  "BG_baseline_reversal_supported",
  "BG_response_restoration_supported",
  "Directional_support_score",
  "RNA_chromatin_priority"
)

missing_cols <- setdiff(
  required_cols,
  names(dat)
)

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# ============================================================
# 4. SAFE LOGICAL CONVERSION
# ============================================================

as_logical_safe <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  y <- toupper(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA,
    length(y)
  )

  out[y == "TRUE"] <- TRUE
  out[y == "FALSE"] <- FALSE

  out
}

dat$BG_baseline_reversal_supported <-
  as_logical_safe(
    dat$BG_baseline_reversal_supported
  )

dat$BG_response_restoration_supported <-
  as_logical_safe(
    dat$BG_response_restoration_supported
  )

# ============================================================
# 5. ADD DISPLAY FLAGS
# ============================================================

dat$Day6_memory_FDR05 <-
  dat$Day6_memory_FDR_kinase < 0.05

dat$Day1_nominal_P05 <-
  dat$Day1_P < 0.05

dat$Blunted_restimulation <-
  dat$Response_pattern ==
  "blunted_same_direction"

dat$Interaction_nominal_P05 <-
  dat$Interaction_P < 0.05

dat$Interaction_genome_FDR05 <-
  dat$Interaction_FDR_genome < 0.05

# ============================================================
# 6. ORDER BY RNA PRIORITY
# ============================================================

priority_rank <- c(
  HIGH = 1,
  INTERMEDIATE = 2,
  EXPLORATORY = 3
)

dat$priority_rank <-
  unname(
    priority_rank[
      dat$RNA_chromatin_priority
    ]
  )

dat <- dat[
  order(
    dat$priority_rank,
    -dat$Directional_support_score,
    dat$Day6_memory_FDR_kinase
  ),
  ,
  drop = FALSE
]

dat$priority_rank <- NULL

write.csv(
  dat,
  file.path(
    outdir,
    "06_final_kinase_figure_data.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 7. BUILD EVIDENCE MATRIX
# ============================================================

evidence_matrix <- cbind(

  Day6_memory =
    dat$Day6_memory_FDR05,

  Day1_nominal =
    dat$Day1_nominal_P05,

  Blunted_restim =
    dat$Blunted_restimulation,

  BG_baseline =
    dat$BG_baseline_reversal_supported,

  BG_response =
    dat$BG_response_restoration_supported,

  Interaction_nominal =
    dat$Interaction_nominal_P05
)

evidence_matrix[
  is.na(evidence_matrix)
] <- FALSE

storage.mode(
  evidence_matrix
) <- "numeric"

rownames(
  evidence_matrix
) <- dat$SYMBOL

# ============================================================
# 8. FIGURE 06A — EVIDENCE MATRIX
# ============================================================

png(
  file.path(
    figdir,
    "06A_integrated_kinase_evidence_matrix.png"
  ),
  width = 2500,
  height = 2400,
  res = 220
)

par(
  mar = c(
    11,
    10,
    5,
    3
  )
)

image(
  x = seq_len(
    ncol(evidence_matrix)
  ),
  y = seq_len(
    nrow(evidence_matrix)
  ),
  z = t(
    evidence_matrix[
      nrow(evidence_matrix):1,
      ,
      drop = FALSE
    ]
  ),
  col = c(
    "white",
    "black"
  ),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main =
    "Integrated RNA evidence for Day-6 memory kinases"
)

axis(
  1,
  at = seq_len(
    ncol(evidence_matrix)
  ),
  labels = c(
    "Day6\nmemory",
    "Day1\nnominal",
    "Blunted\nrestim",
    "BG\nbaseline",
    "BG\nresponse",
    "Interaction\nnominal"
  ),
  las = 2,
  cex.axis = 0.8
)

axis(
  2,
  at = seq_len(
    nrow(evidence_matrix)
  ),
  labels = rev(
    rownames(
      evidence_matrix
    )
  ),
  las = 2,
  cex.axis = 0.8
)

abline(
  v = seq(
    1.5,
    ncol(evidence_matrix) - 0.5,
    by = 1
  )
)

abline(
  h = seq(
    1.5,
    nrow(evidence_matrix) - 0.5,
    by = 1
  )
)

mtext(
  paste(
    "Black = supporting evidence;",
    "interaction evidence is nominal only"
  ),
  side = 1,
  line = 9,
  cex = 0.7
)

dev.off()

# ============================================================
# 9. FIGURE 06B — DAY-6 MEMORY EFFECTS
# ============================================================

plot_dat <- dat[
  nrow(dat):1,
  ,
  drop = FALSE
]

ypos <- seq_len(
  nrow(plot_dat)
)

xmax <- max(
  abs(
    plot_dat$Day6_memory_log2FC
  ),
  na.rm = TRUE
)

png(
  file.path(
    figdir,
    "06B_Day6_memory_kinase_effects.png"
  ),
  width = 2100,
  height = 2300,
  res = 220
)

par(
  mar = c(
    6,
    10,
    5,
    3
  )
)

plot(
  plot_dat$Day6_memory_log2FC,
  ypos,
  pch = 16,
  cex = 1.1,
  xlim = c(
    -xmax,
    xmax
  ),
  yaxt = "n",
  xlab =
    "Day-6 LPS vs RPMI log2 fold change",
  ylab = "",
  main =
    "Persistent Day-6 memory-associated kinase effects"
)

axis(
  2,
  at = ypos,
  labels =
    plot_dat$SYMBOL,
  las = 2,
  cex.axis = 0.8
)

abline(
  v = 0,
  lty = 2
)

segments(
  x0 = 0,
  y0 = ypos,
  x1 =
    plot_dat$Day6_memory_log2FC,
  y1 = ypos
)

text(
  x =
    plot_dat$Day6_memory_log2FC,
  y = ypos,
  labels =
    plot_dat$RNA_chromatin_priority,
  pos = ifelse(
    plot_dat$Day6_memory_log2FC >= 0,
    4,
    2
  ),
  cex = 0.6
)

dev.off()

# ============================================================
# 10. FIGURE 06C — DAY-1 VS DAY-6
# ============================================================

png(
  file.path(
    figdir,
    "06C_Day1_vs_Day6_kinase_effects.png"
  ),
  width = 2100,
  height = 1800,
  res = 220
)

plot(
  dat$Day1_log2FC,
  dat$Day6_memory_log2FC,
  pch = 16,
  cex = 1.1,
  xlab =
    "Day-1 LPS vs RPMI log2 fold change",
  ylab =
    "Day-6 LPS vs RPMI log2 fold change",
  main =
    "Initial response versus persistent Day-6 memory effect"
)

abline(
  h = 0,
  lty = 2
)

abline(
  v = 0,
  lty = 2
)

text(
  dat$Day1_log2FC,
  dat$Day6_memory_log2FC,
  labels =
    dat$SYMBOL,
  pos = 3,
  cex = 0.65
)

dev.off()

# ============================================================
# 11. FIGURE 06D — FINAL PROFESSOR-FACING SUMMARY
# ============================================================

png(
  file.path(
    figdir,
    "06D_FINAL_integrated_kinase_summary.png"
  ),
  width = 3600,
  height = 2400,
  res = 240
)

layout(
  matrix(
    c(
      1,
      2
    ),
    nrow = 1
  ),
  widths = c(
    1.1,
    1.6
  )
)

# ---------------- PANEL A ----------------

par(
  mar = c(
    7,
    10,
    6,
    2
  )
)

plot(
  plot_dat$Day6_memory_log2FC,
  ypos,
  pch = 16,
  cex = 1,
  xlim = c(
    -xmax,
    xmax
  ),
  yaxt = "n",
  xlab =
    "Day-6 LPS vs RPMI log2FC",
  ylab = "",
  main =
    "A. Persistent memory effect"
)

axis(
  2,
  at = ypos,
  labels =
    plot_dat$SYMBOL,
  las = 2,
  cex.axis = 0.72
)

abline(
  v = 0,
  lty = 2
)

segments(
  0,
  ypos,
  plot_dat$Day6_memory_log2FC,
  ypos
)

# ---------------- PANEL B ----------------

par(
  mar = c(
    11,
    4,
    6,
    3
  )
)

image(
  x = seq_len(
    ncol(evidence_matrix)
  ),
  y = seq_len(
    nrow(evidence_matrix)
  ),
  z = t(
    evidence_matrix[
      nrow(evidence_matrix):1,
      ,
      drop = FALSE
    ]
  ),
  col = c(
    "white",
    "black"
  ),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main =
    "B. Integrated supporting evidence"
)

axis(
  1,
  at = seq_len(
    ncol(evidence_matrix)
  ),
  labels = c(
    "Day6\nmemory",
    "Day1\nnominal",
    "Blunted\nrestim",
    "BG\nbaseline",
    "BG\nresponse",
    "Interaction\nnominal"
  ),
  las = 2,
  cex.axis = 0.72
)

abline(
  v = seq(
    1.5,
    ncol(evidence_matrix) - 0.5,
    by = 1
  )
)

abline(
  h = seq(
    1.5,
    nrow(evidence_matrix) - 0.5,
    by = 1
  )
)

mtext(
  paste(
    "Black = supporting evidence;",
    "interaction evidence is nominal only"
  ),
  side = 1,
  line = 9,
  cex = 0.65
)

dev.off()

# ============================================================
# 12. HIGH-PRIORITY TABLE
# ============================================================

high <- dat[
  dat$RNA_chromatin_priority ==
    "HIGH",
  ,
  drop = FALSE
]

write.csv(
  high,
  file.path(
    outdir,
    "06_HIGH_priority_kinases.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 13. SUMMARY COUNTS
# ============================================================

n_high <- sum(
  dat$RNA_chromatin_priority ==
    "HIGH",
  na.rm = TRUE
)

n_intermediate <- sum(
  dat$RNA_chromatin_priority ==
    "INTERMEDIATE",
  na.rm = TRUE
)

n_exploratory <- sum(
  dat$RNA_chromatin_priority ==
    "EXPLORATORY",
  na.rm = TRUE
)

n_interaction_fdr <- sum(
  dat$Interaction_genome_FDR05,
  na.rm = TRUE
)

summary_lines <- c(

  "GSE85246 — Step 06 Final Integrated Kinase Figures",

  "",

  paste(
    "Total memory kinases:",
    nrow(dat)
  ),

  paste(
    "HIGH priority:",
    n_high
  ),

  paste(
    "INTERMEDIATE priority:",
    n_intermediate
  ),

  paste(
    "EXPLORATORY priority:",
    n_exploratory
  ),

  "",

  paste(
    "Genome-wide significant restimulation interactions:",
    n_interaction_fdr
  ),

  "",

  "HIGH-priority kinase symbols:",

  paste(
    high$SYMBOL,
    collapse = ", "
  ),

  "",

  paste(
    "IMPORTANT:",
    "priority categories represent integrated directional",
    "RNA evidence and are not new statistical significance tests."
  )
)

writeLines(
  summary_lines,
  file.path(
    outdir,
    "06_final_kinase_figure_summary.txt"
  )
)

# ============================================================
# 14. FINAL CONSOLE OUTPUT
# ============================================================

cat("\n========================================\n")
cat("06 FINAL KINASE FIGURES COMPLETE\n")
cat("========================================\n\n")

cat(
  "Total memory kinases:",
  nrow(dat),
  "\n"
)

cat(
  "HIGH:",
  n_high,
  "\n"
)

cat(
  "INTERMEDIATE:",
  n_intermediate,
  "\n"
)

cat(
  "EXPLORATORY:",
  n_exploratory,
  "\n"
)

cat(
  "Genome-wide significant interaction:",
  n_interaction_fdr,
  "\n"
)

cat(
  "\nHIGH-priority kinases:\n"
)

cat(
  paste(
    high$SYMBOL,
    collapse = ", "
  ),
  "\n"
)

cat(
  "\nFigures created:\n"
)

cat(
  "06A_integrated_kinase_evidence_matrix.png\n"
)

cat(
  "06B_Day6_memory_kinase_effects.png\n"
)

cat(
  "06C_Day1_vs_Day6_kinase_effects.png\n"
)

cat(
  "06D_FINAL_integrated_kinase_summary.png\n"
)

cat(
  "\nNext canonical step:\n"
)

cat(
  "scripts/07_epigenetic_persistence.R\n"
)

cat("\n========================================\n")
