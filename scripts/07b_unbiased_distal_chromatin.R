# ============================================================
# 07b_unbiased_distal_chromatin.R
#
# GSE85246 / GSE85245
#
# CORRECTED DISTAL CHROMATIN ANALYSIS
#
# Purpose:
#   Remove the selection bias created by choosing the "best"
#   distal window AFTER looking at Day1/Day6/BG directional
#   effects.
#
# Strategy:
#
#   1. Read the already-quantified Step-07 BigWig matrix.
#   2. Consider distal windows only.
#   3. Select ONE distal window per gene using CONDITION-BLIND
#      chromatin signal:
#         - mean H3K27ac across all H3K27ac tracks
#         - mean H3K4me1 across all H3K4me1 tracks
#         - convert each mark to a global distal percentile
#         - average the two percentiles
#   4. Freeze that one window.
#   5. Only then calculate:
#         Day1 LPS - RPMI
#         Day6 LPS - RPMI
#         Day6 RESCUE - LPS
#   6. Integrate the fixed distal locus with the already-valid
#      single predefined promoter result.
#
# IMPORTANT:
#   - No BigWig files are reread.
#   - No peak calling is performed.
#   - This remains directional/descriptive chromatin evidence.
#   - Day1 and Day6 ChIP donor sets are not longitudinally paired.
# ============================================================

rm(list = ls())

# ============================================================
# 1. INPUT FILES
# ============================================================

signal_file <- paste0(
  "results/07_epigenetic_persistence/",
  "07_all_region_BigWig_signal.csv"
)

promoter_file <- paste0(
  "results/07_epigenetic_persistence/",
  "07_promoter_integrated_marks.csv"
)

rna_file <- paste0(
  "results/05_integrated_trajectory/",
  "05_integrated_kinase_trajectory.csv"
)

required_files <- c(
  signal_file,
  promoter_file,
  rna_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required file(s):",
      paste(missing_files, collapse = "\n")
    )
  )
}

# ============================================================
# 2. OUTPUT DIRECTORIES
# ============================================================

outdir <- "results/07b_unbiased_distal_chromatin"
figdir <- "figures/07b_unbiased_distal_chromatin"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 3. LOAD DATA
# ============================================================

sig <- read.csv(
  signal_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

promoter <- read.csv(
  promoter_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rna <- read.csv(
  rna_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("07b — UNBIASED DISTAL CHROMATIN\n")
cat("========================================\n\n")

cat("Rows in quantified Step-07 matrix:", nrow(sig), "\n")
cat("RNA memory kinases:", nrow(rna), "\n")

# ============================================================
# 4. REQUIRED METADATA COLUMNS
# ============================================================

required_meta <- c(
  "region_id",
  "ENSEMBL",
  "SYMBOL",
  "RNA_priority",
  "RNA_support_score",
  "region_type",
  "chr",
  "start",
  "end",
  "distance_to_TSS"
)

missing_meta <- setdiff(
  required_meta,
  names(sig)
)

if (length(missing_meta) > 0) {
  stop(
    paste(
      "Missing Step-07 metadata columns:",
      paste(missing_meta, collapse = ", ")
    )
  )
}

# ============================================================
# 5. EXPECTED TRACK NAMES
# ============================================================

h3k27ac_tracks <- c(
  "D1_RPMI_r1_H3K27ac",
  "D1_RPMI_r2_H3K27ac",
  "D1_LPS_r1_H3K27ac",
  "D1_LPS_r2_H3K27ac",
  "D6_RPMI_d73_H3K27ac",
  "D6_RPMI_d74_H3K27ac",
  "D6_LPS_d73_H3K27ac",
  "D6_LPS_d74_H3K27ac",
  "D6_RESCUE_d73_H3K27ac",
  "D6_RESCUE_d74_H3K27ac"
)

h3k4me1_tracks <- c(
  "D1_RPMI_r1_H3K4me1",
  "D1_RPMI_r2_H3K4me1",
  "D1_LPS_r1_H3K4me1",
  "D1_LPS_r2_H3K4me1",
  "D6_RPMI_d73_H3K4me1",
  "D6_RPMI_d74_H3K4me1",
  "D6_LPS_d73_H3K4me1",
  "D6_LPS_d74_H3K4me1",
  "D6_RESCUE_d73_H3K4me1",
  "D6_RESCUE_d74_H3K4me1"
)

required_tracks <- c(
  h3k27ac_tracks,
  h3k4me1_tracks
)

missing_tracks <- setdiff(
  required_tracks,
  names(sig)
)

if (length(missing_tracks) > 0) {
  stop(
    paste(
      "Missing expected BigWig signal columns:",
      paste(missing_tracks, collapse = ", ")
    )
  )
}

# ============================================================
# 6. DISTAL WINDOWS ONLY
# ============================================================

distal <- sig[
  sig$region_type == "distal",
  ,
  drop = FALSE
]

if (nrow(distal) == 0) {
  stop("No distal windows found.")
}

cat("Distal windows available:", nrow(distal), "\n")

# ============================================================
# 7. CONDITION-BLIND DISTAL SIGNAL
#
# Selection does NOT use:
#   - LPS-RPMI direction
#   - Day1/Day6 persistence
#   - BG reversal
#
# It only uses overall chromatin abundance across ALL tracks.
# ============================================================

distal$H3K27ac_condition_blind_mean <- rowMeans(
  distal[
    ,
    h3k27ac_tracks,
    drop = FALSE
  ],
  na.rm = TRUE
)

distal$H3K4me1_condition_blind_mean <- rowMeans(
  distal[
    ,
    h3k4me1_tracks,
    drop = FALSE
  ],
  na.rm = TRUE
)

percentile_rank <- function(x) {

  ok <- is.finite(x)

  out <- rep(
    NA_real_,
    length(x)
  )

  if (sum(ok) == 0) {
    return(out)
  }

  out[ok] <- rank(
    x[ok],
    ties.method = "average"
  ) / sum(ok)

  out
}

distal$H3K27ac_condition_blind_percentile <-
  percentile_rank(
    distal$H3K27ac_condition_blind_mean
  )

distal$H3K4me1_condition_blind_percentile <-
  percentile_rank(
    distal$H3K4me1_condition_blind_mean
  )

distal$condition_blind_combined_score <- rowMeans(
  cbind(
    distal$H3K27ac_condition_blind_percentile,
    distal$H3K4me1_condition_blind_percentile
  ),
  na.rm = TRUE
)

# ============================================================
# 8. SELECT EXACTLY ONE DISTAL WINDOW PER GENE
#
# Primary rule:
#   highest condition-blind combined score
#
# Tie-breaker:
#   closest absolute distance to TSS
# ============================================================

distal$abs_distance_to_TSS <- abs(
  distal$distance_to_TSS
)

distal_ordered <- distal[
  order(
    distal$SYMBOL,
    -distal$condition_blind_combined_score,
    distal$abs_distance_to_TSS,
    distal$chr,
    distal$start
  ),
  ,
  drop = FALSE
]

selected <- distal_ordered[
  !duplicated(
    distal_ordered$SYMBOL
  ),
  ,
  drop = FALSE
]

if (nrow(selected) != nrow(rna)) {

  cat("\nSelected distal symbols:\n")
  print(selected$SYMBOL)

  stop(
    paste(
      "Expected one fixed distal window for each RNA kinase.",
      "Selected:",
      nrow(selected),
      "RNA kinases:",
      nrow(rna)
    )
  )
}

write.csv(
  selected,
  file.path(
    outdir,
    "07b_condition_blind_selected_distal_windows.csv"
  ),
  row.names = FALSE
)

cat(
  "Condition-blind fixed distal windows selected:",
  nrow(selected),
  "\n"
)

# ============================================================
# 9. EFFECT-CALCULATION HELPER
# ============================================================

calc_mark <- function(
    d,
    mark
) {

  if (mark == "H3K27ac") {

    d1_rpmi <- rowMeans(
      d[
        ,
        c(
          "D1_RPMI_r1_H3K27ac",
          "D1_RPMI_r2_H3K27ac"
        ),
        drop = FALSE
      ]
    )

    d1_lps <- rowMeans(
      d[
        ,
        c(
          "D1_LPS_r1_H3K27ac",
          "D1_LPS_r2_H3K27ac"
        ),
        drop = FALSE
      ]
    )

    d6_rpmi <- rowMeans(
      d[
        ,
        c(
          "D6_RPMI_d73_H3K27ac",
          "D6_RPMI_d74_H3K27ac"
        ),
        drop = FALSE
      ]
    )

    d6_lps <- rowMeans(
      d[
        ,
        c(
          "D6_LPS_d73_H3K27ac",
          "D6_LPS_d74_H3K27ac"
        ),
        drop = FALSE
      ]
    )

    d6_rescue <- rowMeans(
      d[
        ,
        c(
          "D6_RESCUE_d73_H3K27ac",
          "D6_RESCUE_d74_H3K27ac"
        ),
        drop = FALSE
      ]
    )

  } else if (mark == "H3K4me1") {

    d1_rpmi <- rowMeans(
      d[
        ,
        c(
          "D1_RPMI_r1_H3K4me1",
          "D1_RPMI_r2_H3K4me1"
        ),
        drop = FALSE
      ]
    )

    d1_lps <- rowMeans(
      d[
        ,
        c(
          "D1_LPS_r1_H3K4me1",
          "D1_LPS_r2_H3K4me1"
        ),
        drop = FALSE
      ]
    )

    d6_rpmi <- rowMeans(
      d[
        ,
        c(
          "D6_RPMI_d73_H3K4me1",
          "D6_RPMI_d74_H3K4me1"
        ),
        drop = FALSE
      ]
    )

    d6_lps <- rowMeans(
      d[
        ,
        c(
          "D6_LPS_d73_H3K4me1",
          "D6_LPS_d74_H3K4me1"
        ),
        drop = FALSE
      ]
    )

    d6_rescue <- rowMeans(
      d[
        ,
        c(
          "D6_RESCUE_d73_H3K4me1",
          "D6_RESCUE_d74_H3K4me1"
        ),
        drop = FALSE
      ]
    )

  } else {

    stop("Unknown mark.")
  }

  d1_delta <- d1_lps - d1_rpmi
  d6_delta <- d6_lps - d6_rpmi
  bg_shift <- d6_rescue - d6_lps

  same_direction <-
    is.finite(d1_delta) &
    is.finite(d6_delta) &
    (d1_delta * d6_delta > 0)

  bg_opposes_day6 <-
    is.finite(bg_shift) &
    is.finite(d6_delta) &
    (bg_shift * d6_delta < 0)

  rescue_closer <-
    abs(d6_rescue - d6_rpmi) <
    abs(d6_lps - d6_rpmi)

  persistent_bg_reversible <-
    same_direction &
    bg_opposes_day6 &
    rescue_closer

  data.frame(
    Day1_RPMI = d1_rpmi,
    Day1_LPS = d1_lps,
    Day1_delta = d1_delta,
    Day6_RPMI = d6_rpmi,
    Day6_LPS = d6_lps,
    Day6_RESCUE = d6_rescue,
    Day6_delta = d6_delta,
    BG_shift = bg_shift,
    cross_time_same_direction = same_direction,
    BG_shift_opposes_Day6 = bg_opposes_day6,
    BG_closer_to_RPMI = rescue_closer,
    persistent_BG_reversible = persistent_bg_reversible,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 10. TEST THE FROZEN DISTAL WINDOWS
# ============================================================

ac <- calc_mark(
  selected,
  "H3K27ac"
)

me1 <- calc_mark(
  selected,
  "H3K4me1"
)

fixed_distal <- selected[
  ,
  c(
    "region_id",
    "ENSEMBL",
    "SYMBOL",
    "RNA_priority",
    "RNA_support_score",
    "chr",
    "start",
    "end",
    "distance_to_TSS",
    "H3K27ac_condition_blind_mean",
    "H3K4me1_condition_blind_mean",
    "H3K27ac_condition_blind_percentile",
    "H3K4me1_condition_blind_percentile",
    "condition_blind_combined_score"
  )
]

fixed_distal$H3K27ac_Day1_delta <-
  ac$Day1_delta

fixed_distal$H3K27ac_Day6_delta <-
  ac$Day6_delta

fixed_distal$H3K27ac_BG_shift <-
  ac$BG_shift

fixed_distal$H3K27ac_cross_time_same_direction <-
  ac$cross_time_same_direction

fixed_distal$H3K27ac_BG_shift_opposes_Day6 <-
  ac$BG_shift_opposes_Day6

fixed_distal$H3K27ac_BG_closer_to_RPMI <-
  ac$BG_closer_to_RPMI

fixed_distal$H3K27ac_persistent_BG_reversible <-
  ac$persistent_BG_reversible

fixed_distal$H3K4me1_Day1_delta <-
  me1$Day1_delta

fixed_distal$H3K4me1_Day6_delta <-
  me1$Day6_delta

fixed_distal$H3K4me1_BG_shift <-
  me1$BG_shift

fixed_distal$H3K4me1_cross_time_same_direction <-
  me1$cross_time_same_direction

fixed_distal$H3K4me1_BG_shift_opposes_Day6 <-
  me1$BG_shift_opposes_Day6

fixed_distal$H3K4me1_BG_closer_to_RPMI <-
  me1$BG_closer_to_RPMI

fixed_distal$H3K4me1_persistent_BG_reversible <-
  me1$persistent_BG_reversible

fixed_distal$Distal_mark_support_count <-
  as.integer(
    fixed_distal$H3K27ac_persistent_BG_reversible
  ) +
  as.integer(
    fixed_distal$H3K4me1_persistent_BG_reversible
  )

fixed_distal$Distal_any_mark_support <-
  fixed_distal$Distal_mark_support_count >= 1

fixed_distal$Distal_dual_mark_support <-
  fixed_distal$Distal_mark_support_count == 2

write.csv(
  fixed_distal,
  file.path(
    outdir,
    "07b_fixed_distal_chromatin_results.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 11. VALIDATE PROMOTER TABLE
# ============================================================

required_promoter_cols <- c(
  "ENSEMBL",
  "SYMBOL",
  "mark_support_count",
  "dual_mark_support",
  "any_mark_support"
)

missing_promoter_cols <- setdiff(
  required_promoter_cols,
  names(promoter)
)

if (length(missing_promoter_cols) > 0) {
  stop(
    paste(
      "Promoter table missing columns:",
      paste(missing_promoter_cols, collapse = ", ")
    )
  )
}

if (nrow(promoter) != 24) {
  stop(
    paste(
      "Expected 24 promoter rows, found",
      nrow(promoter)
    )
  )
}

# ============================================================
# 12. BUILD CORRECTED GENE-LEVEL TABLE
# ============================================================

gene <- rna[
  ,
  c(
    "ENSEMBL",
    "SYMBOL",
    "RNA_chromatin_priority",
    "Directional_support_score",
    "Day6_memory_log2FC",
    "Day6_memory_FDR_kinase"
  )
]

pidx <- match(
  gene$ENSEMBL,
  promoter$ENSEMBL
)

didx <- match(
  gene$ENSEMBL,
  fixed_distal$ENSEMBL
)

if (anyNA(pidx)) {
  stop("One or more RNA kinases are missing from promoter results.")
}

if (anyNA(didx)) {
  stop("One or more RNA kinases are missing from fixed distal results.")
}

gene$Promoter_mark_support_count <-
  promoter$mark_support_count[pidx]

gene$Promoter_any_mark_support <-
  promoter$any_mark_support[pidx]

gene$Promoter_dual_mark_support <-
  promoter$dual_mark_support[pidx]

gene$Fixed_distal_region_id <-
  fixed_distal$region_id[didx]

gene$Fixed_distal_chr <-
  fixed_distal$chr[didx]

gene$Fixed_distal_start <-
  fixed_distal$start[didx]

gene$Fixed_distal_end <-
  fixed_distal$end[didx]

gene$Fixed_distal_distance_to_TSS <-
  fixed_distal$distance_to_TSS[didx]

gene$Fixed_distal_selection_score <-
  fixed_distal$condition_blind_combined_score[didx]

gene$Distal_mark_support_count <-
  fixed_distal$Distal_mark_support_count[didx]

gene$Distal_any_mark_support <-
  fixed_distal$Distal_any_mark_support[didx]

gene$Distal_dual_mark_support <-
  fixed_distal$Distal_dual_mark_support[didx]

gene$Any_corrected_chromatin_support <-
  gene$Promoter_any_mark_support |
  gene$Distal_any_mark_support

gene$Any_corrected_dual_mark_support <-
  gene$Promoter_dual_mark_support |
  gene$Distal_dual_mark_support

gene$Corrected_chromatin_support_class <- ifelse(
  gene$Any_corrected_dual_mark_support,
  "DUAL_MARK_PERSISTENT_BG_REVERSIBLE",
  ifelse(
    gene$Any_corrected_chromatin_support,
    "SINGLE_MARK_PERSISTENT_BG_REVERSIBLE",
    "NO_DIRECTIONAL_CHROMATIN_SUPPORT"
  )
)

gene$Corrected_RNA_plus_chromatin_class <- ifelse(
  gene$RNA_chromatin_priority == "HIGH" &
    gene$Any_corrected_dual_mark_support,
  "TOP_RNA_HIGH_PLUS_DUAL_CHROMATIN",
  ifelse(
    gene$RNA_chromatin_priority == "HIGH" &
      gene$Any_corrected_chromatin_support,
    "RNA_HIGH_PLUS_SINGLE_MARK_CHROMATIN",
    ifelse(
      gene$RNA_chromatin_priority == "HIGH",
      "RNA_HIGH_NO_DIRECTIONAL_CHROMATIN",
      ifelse(
        gene$Any_corrected_dual_mark_support,
        "NONHIGH_RNA_PLUS_DUAL_CHROMATIN",
        ifelse(
          gene$Any_corrected_chromatin_support,
          "NONHIGH_RNA_PLUS_SINGLE_MARK_CHROMATIN",
          "EXPLORATORY"
        )
      )
    )
  )
)

priority_rank <- c(
  "TOP_RNA_HIGH_PLUS_DUAL_CHROMATIN" = 1,
  "RNA_HIGH_PLUS_SINGLE_MARK_CHROMATIN" = 2,
  "RNA_HIGH_NO_DIRECTIONAL_CHROMATIN" = 3,
  "NONHIGH_RNA_PLUS_DUAL_CHROMATIN" = 4,
  "NONHIGH_RNA_PLUS_SINGLE_MARK_CHROMATIN" = 5,
  "EXPLORATORY" = 6
)

gene$priority_rank <- unname(
  priority_rank[
    gene$Corrected_RNA_plus_chromatin_class
  ]
)

gene <- gene[
  order(
    gene$priority_rank,
    gene$Day6_memory_FDR_kinase
  ),
  ,
  drop = FALSE
]

gene$priority_rank <- NULL

write.csv(
  gene,
  file.path(
    outdir,
    "07b_corrected_gene_level_chromatin_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 13. HIGH RNA SUBSET
# ============================================================

high <- gene[
  gene$RNA_chromatin_priority == "HIGH",
  ,
  drop = FALSE
]

write.csv(
  high,
  file.path(
    outdir,
    "07b_HIGH_RNA_chromatin_summary.csv"
  ),
  row.names = FALSE
)

top <- gene[
  gene$Corrected_RNA_plus_chromatin_class ==
    "TOP_RNA_HIGH_PLUS_DUAL_CHROMATIN",
  ,
  drop = FALSE
]

write.csv(
  top,
  file.path(
    outdir,
    "07b_TOP_RNA_HIGH_plus_dual_chromatin.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 14. SIMPLE SUPPORT FIGURE
# ============================================================

mat <- cbind(
  Promoter_H3K27ac =
    promoter$H3K27ac_persistent_BG_reversible[
      match(
        gene$ENSEMBL,
        promoter$ENSEMBL
      )
    ],

  Promoter_H3K4me1 =
    promoter$H3K4me1_persistent_BG_reversible[
      match(
        gene$ENSEMBL,
        promoter$ENSEMBL
      )
    ],

  Fixed_Distal_H3K27ac =
    fixed_distal$H3K27ac_persistent_BG_reversible[
      match(
        gene$ENSEMBL,
        fixed_distal$ENSEMBL
      )
    ],

  Fixed_Distal_H3K4me1 =
    fixed_distal$H3K4me1_persistent_BG_reversible[
      match(
        gene$ENSEMBL,
        fixed_distal$ENSEMBL
      )
    ]
)

mat[is.na(mat)] <- FALSE
storage.mode(mat) <- "numeric"
rownames(mat) <- gene$SYMBOL

png(
  file.path(
    figdir,
    "07b_corrected_chromatin_support_matrix.png"
  ),
  width = 2200,
  height = 2400,
  res = 220
)

par(
  mar = c(11, 10, 5, 3)
)

image(
  x = seq_len(ncol(mat)),
  y = seq_len(nrow(mat)),
  z = t(
    mat[
      nrow(mat):1,
      ,
      drop = FALSE
    ]
  ),
  col = c("white", "black"),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main =
    "Corrected chromatin support at fixed memory-kinase loci"
)

axis(
  1,
  at = seq_len(ncol(mat)),
  labels = colnames(mat),
  las = 2,
  cex.axis = 0.75
)

axis(
  2,
  at = seq_len(nrow(mat)),
  labels = rev(rownames(mat)),
  las = 2,
  cex.axis = 0.75
)

abline(
  v = seq(
    1.5,
    ncol(mat) - 0.5,
    by = 1
  )
)

abline(
  h = seq(
    1.5,
    nrow(mat) - 0.5,
    by = 1
  )
)

dev.off()

# ============================================================
# 15. SUMMARY COUNTS
# ============================================================

n_prom_any <- sum(
  gene$Promoter_any_mark_support,
  na.rm = TRUE
)

n_prom_dual <- sum(
  gene$Promoter_dual_mark_support,
  na.rm = TRUE
)

n_dist_any <- sum(
  gene$Distal_any_mark_support,
  na.rm = TRUE
)

n_dist_dual <- sum(
  gene$Distal_dual_mark_support,
  na.rm = TRUE
)

n_any <- sum(
  gene$Any_corrected_chromatin_support,
  na.rm = TRUE
)

n_dual <- sum(
  gene$Any_corrected_dual_mark_support,
  na.rm = TRUE
)

n_high_any <- sum(
  high$Any_corrected_chromatin_support,
  na.rm = TRUE
)

n_high_dual <- sum(
  high$Any_corrected_dual_mark_support,
  na.rm = TRUE
)

summary_lines <- c(
  "GSE85246 — Step 07b Corrected Distal Chromatin",
  "",
  paste("RNA memory kinases:", nrow(gene)),
  paste("RNA HIGH-priority kinases:", nrow(high)),
  "",
  paste("Promoter any-mark support:", n_prom_any),
  paste("Promoter dual-mark support:", n_prom_dual),
  "",
  paste("Fixed distal any-mark support:", n_dist_any),
  paste("Fixed distal dual-mark support:", n_dist_dual),
  "",
  paste("Any corrected chromatin support:", n_any),
  paste("Any corrected dual-mark support:", n_dual),
  "",
  paste("RNA HIGH with any corrected chromatin support:", n_high_any),
  paste("RNA HIGH with corrected dual-mark support:", n_high_dual),
  "",
  paste("TOP RNA-HIGH + dual-chromatin candidates:", nrow(top)),
  if (nrow(top) > 0) {
    paste(top$SYMBOL, collapse = ", ")
  } else {
    "None"
  },
  "",
  paste(
    "IMPORTANT:",
    "The distal locus was selected using condition-blind total",
    "chromatin signal before testing Day1/Day6/BG direction."
  ),
  paste(
    "This removes the previous best-window directional",
    "selection bias."
  )
)

writeLines(
  summary_lines,
  file.path(
    outdir,
    "07b_corrected_chromatin_summary.txt"
  )
)

# ============================================================
# 16. FINAL CONSOLE REPORT
# ============================================================

cat("\n========================================\n")
cat("07b CORRECTED DISTAL CHROMATIN COMPLETE\n")
cat("========================================\n\n")

cat(
  "RNA memory kinases:",
  nrow(gene),
  "\n"
)

cat(
  "RNA HIGH-priority kinases:",
  nrow(high),
  "\n"
)

cat(
  "\nPromoter any-mark support:",
  n_prom_any,
  "/",
  nrow(gene),
  "\n"
)

cat(
  "Promoter dual-mark support:",
  n_prom_dual,
  "/",
  nrow(gene),
  "\n"
)

cat(
  "\nFixed distal any-mark support:",
  n_dist_any,
  "/",
  nrow(gene),
  "\n"
)

cat(
  "Fixed distal dual-mark support:",
  n_dist_dual,
  "/",
  nrow(gene),
  "\n"
)

cat(
  "\nAny corrected chromatin support:",
  n_any,
  "/",
  nrow(gene),
  "\n"
)

cat(
  "Any corrected dual-mark support:",
  n_dual,
  "/",
  nrow(gene),
  "\n"
)

cat(
  "\nRNA HIGH with any corrected chromatin support:",
  n_high_any,
  "/",
  nrow(high),
  "\n"
)

cat(
  "RNA HIGH with corrected dual-mark support:",
  n_high_dual,
  "/",
  nrow(high),
  "\n"
)

cat(
  "\nTOP RNA-HIGH + dual-chromatin candidates:",
  nrow(top),
  "\n"
)

if (nrow(top) > 0) {
  cat(
    paste(
      top$SYMBOL,
      collapse = ", "
    ),
    "\n"
  )
}

cat(
  "\nIMPORTANT:\n",
  "Distal loci were selected condition-blind before ",
  "testing Day1/Day6/BG direction.\n",
  sep = ""
)

cat(
  "\nNext step after review:\n",
  "scripts/08_make_final_chromatin_figures.R\n",
  sep = ""
)

cat("\n========================================\n")
