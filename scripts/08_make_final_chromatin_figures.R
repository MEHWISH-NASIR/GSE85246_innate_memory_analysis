# ============================================================
# 08_make_final_chromatin_figures.R
#
# GSE85246 / GSE85245
#
# FINAL CORRECTED CHROMATIN FIGURES
#
# Uses ONLY the corrected Step-07b distal analysis plus the
# predefined promoter analysis from Step 07.
#
# No new statistical testing is performed here.
#
# Main goals:
#   1. Show corrected chromatin support across all 24 memory kinases
#   2. Show chromatin status of the 10 HIGH RNA-priority kinases
#   3. Visualize directional chromatin behavior of the strongest
#      integrated candidates
#   4. Produce a professor-facing final RNA + chromatin summary
#
# IMPORTANT:
#   Chromatin evidence remains directional/descriptive.
#   Day1 and Day6 ChIP donor sets are not longitudinally paired.
# ============================================================

rm(list = ls())

# ============================================================
# 1. INPUTS
# ============================================================

gene_file <- paste0(
  "results/07b_unbiased_distal_chromatin/",
  "07b_corrected_gene_level_chromatin_summary.csv"
)

distal_file <- paste0(
  "results/07b_unbiased_distal_chromatin/",
  "07b_fixed_distal_chromatin_results.csv"
)

promoter_file <- paste0(
  "results/07_epigenetic_persistence/",
  "07_promoter_integrated_marks.csv"
)

required_files <- c(
  gene_file,
  distal_file,
  promoter_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required input file(s):",
      paste(missing_files, collapse = "\n")
    )
  )
}

# ============================================================
# 2. OUTPUTS
# ============================================================

outdir <- "results/08_final_chromatin"
figdir <- "figures/08_final_chromatin"

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
# 3. LOAD TABLES
# ============================================================

gene <- read.csv(
  gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

distal <- read.csv(
  distal_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

promoter <- read.csv(
  promoter_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("08 — FINAL CHROMATIN FIGURES\n")
cat("========================================\n\n")

cat("Corrected memory kinases:", nrow(gene), "\n")

# ============================================================
# 4. VALIDATE EXPECTED COLUMNS
# ============================================================

required_gene_cols <- c(
  "ENSEMBL",
  "SYMBOL",
  "RNA_chromatin_priority",
  "Directional_support_score",
  "Day6_memory_log2FC",
  "Day6_memory_FDR_kinase",
  "Promoter_mark_support_count",
  "Promoter_any_mark_support",
  "Promoter_dual_mark_support",
  "Fixed_distal_region_id",
  "Fixed_distal_chr",
  "Fixed_distal_start",
  "Fixed_distal_end",
  "Fixed_distal_distance_to_TSS",
  "Distal_mark_support_count",
  "Distal_any_mark_support",
  "Distal_dual_mark_support",
  "Any_corrected_chromatin_support",
  "Any_corrected_dual_mark_support",
  "Corrected_RNA_plus_chromatin_class"
)

required_distal_cols <- c(
  "ENSEMBL",
  "SYMBOL",
  "region_id",
  "chr",
  "start",
  "end",
  "distance_to_TSS",
  "H3K27ac_Day1_delta",
  "H3K27ac_Day6_delta",
  "H3K27ac_BG_shift",
  "H3K27ac_persistent_BG_reversible",
  "H3K4me1_Day1_delta",
  "H3K4me1_Day6_delta",
  "H3K4me1_BG_shift",
  "H3K4me1_persistent_BG_reversible"
)

required_promoter_cols <- c(
  "ENSEMBL",
  "SYMBOL",
  "region_id",
  "chr",
  "start",
  "end",
  "H3K27ac_Day1_delta",
  "H3K27ac_Day6_delta",
  "H3K27ac_BG_shift",
  "H3K27ac_persistent_BG_reversible",
  "H3K4me1_Day1_delta",
  "H3K4me1_Day6_delta",
  "H3K4me1_BG_shift",
  "H3K4me1_persistent_BG_reversible"
)

check_columns <- function(x, required, label) {

  missing <- setdiff(
    required,
    names(x)
  )

  if (length(missing) > 0) {
    stop(
      paste(
        label,
        "is missing columns:",
        paste(missing, collapse = ", ")
      )
    )
  }
}

check_columns(
  gene,
  required_gene_cols,
  "Corrected gene table"
)

check_columns(
  distal,
  required_distal_cols,
  "Fixed distal table"
)

check_columns(
  promoter,
  required_promoter_cols,
  "Promoter table"
)

if (nrow(gene) != 24) {
  warning(
    paste(
      "Expected 24 memory kinases; found",
      nrow(gene)
    )
  )
}

# ============================================================
# 5. MATCH PROMOTER / DISTAL TO GENE TABLE
# ============================================================

pidx <- match(
  gene$ENSEMBL,
  promoter$ENSEMBL
)

didx <- match(
  gene$ENSEMBL,
  distal$ENSEMBL
)

if (anyNA(pidx)) {
  stop("One or more genes are missing from promoter results.")
}

if (anyNA(didx)) {
  stop("One or more genes are missing from fixed distal results.")
}

promoter_g <- promoter[pidx, , drop = FALSE]
distal_g <- distal[didx, , drop = FALSE]

if (!all(gene$ENSEMBL == promoter_g$ENSEMBL)) {
  stop("Promoter alignment failed.")
}

if (!all(gene$ENSEMBL == distal_g$ENSEMBL)) {
  stop("Distal alignment failed.")
}

# ============================================================
# 6. DEFINE HIGH AND TOP INTEGRATED SETS
# ============================================================

high <- gene[
  gene$RNA_chromatin_priority == "HIGH",
  ,
  drop = FALSE
]

top <- gene[
  gene$Corrected_RNA_plus_chromatin_class ==
    "TOP_RNA_HIGH_PLUS_DUAL_CHROMATIN",
  ,
  drop = FALSE
]

cat("RNA HIGH-priority kinases:", nrow(high), "\n")
cat("TOP RNA-HIGH + dual-chromatin:", nrow(top), "\n")

if (nrow(top) > 0) {
  cat(
    "TOP symbols:",
    paste(top$SYMBOL, collapse = ", "),
    "\n"
  )
}

# ============================================================
# 7. SAVE FINAL SUMMARY TABLES
# ============================================================

write.csv(
  high,
  file.path(
    outdir,
    "08_HIGH_RNA_final_chromatin_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  top,
  file.path(
    outdir,
    "08_TOP_integrated_candidates.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 8. PREPARE IGV LOCI FOR TOP CANDIDATES
# ============================================================

if (nrow(top) > 0) {

  top_pidx <- match(
    top$ENSEMBL,
    promoter$ENSEMBL
  )

  top_didx <- match(
    top$ENSEMBL,
    distal$ENSEMBL
  )

  igv_prom <- data.frame(
    SYMBOL = top$SYMBOL,
    locus_type = "promoter",
    genome = "hg19",
    chr = promoter$chr[top_pidx],
    start = promoter$start[top_pidx],
    end = promoter$end[top_pidx],
    stringsAsFactors = FALSE
  )

  igv_dist <- data.frame(
    SYMBOL = top$SYMBOL,
    locus_type = "fixed_distal",
    genome = "hg19",
    chr = distal$chr[top_didx],
    start = distal$start[top_didx],
    end = distal$end[top_didx],
    stringsAsFactors = FALSE
  )

  igv_loci <- rbind(
    igv_prom,
    igv_dist
  )

  igv_loci$IGV_search <- paste0(
    igv_loci$chr,
    ":",
    igv_loci$start,
    "-",
    igv_loci$end
  )

  write.csv(
    igv_loci,
    file.path(
      outdir,
      "08_TOP_candidate_IGV_loci_hg19.csv"
    ),
    row.names = FALSE
  )

  # BED is zero-based, half-open.
  igv_bed <- data.frame(
    chr = igv_loci$chr,
    start0 = pmax(
      0,
      igv_loci$start - 1
    ),
    end = igv_loci$end,
    name = paste0(
      igv_loci$SYMBOL,
      "_",
      igv_loci$locus_type
    ),
    stringsAsFactors = FALSE
  )

  write.table(
    igv_bed,
    file.path(
      outdir,
      "08_TOP_candidate_IGV_loci_hg19.bed"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

# ============================================================
# 9. FIGURE 08A
#
# Corrected support matrix across all 24 kinases.
# "*" before gene symbol = RNA HIGH priority.
# ============================================================

support_all <- cbind(
  Promoter_H3K27ac =
    as.integer(
      promoter_g$H3K27ac_persistent_BG_reversible
    ),

  Promoter_H3K4me1 =
    as.integer(
      promoter_g$H3K4me1_persistent_BG_reversible
    ),

  Fixed_Distal_H3K27ac =
    as.integer(
      distal_g$H3K27ac_persistent_BG_reversible
    ),

  Fixed_Distal_H3K4me1 =
    as.integer(
      distal_g$H3K4me1_persistent_BG_reversible
    )
)

support_all[is.na(support_all)] <- 0

labels_all <- ifelse(
  gene$RNA_chromatin_priority == "HIGH",
  paste0("* ", gene$SYMBOL),
  gene$SYMBOL
)

png(
  file.path(
    figdir,
    "08A_corrected_chromatin_support_all24.png"
  ),
  width = 2200,
  height = 2500,
  res = 220
)

par(
  mar = c(13, 11, 6, 4)
)

image(
  x = seq_len(ncol(support_all)),
  y = seq_len(nrow(support_all)),
  z = t(
    support_all[
      nrow(support_all):1,
      ,
      drop = FALSE
    ]
  ),
  col = c(
    "#F2F2F2",
    "#25364A"
  ),
  zlim = c(0, 1),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main =
    "Corrected chromatin persistence and BG-reversal support"
)

axis(
  1,
  at = seq_len(ncol(support_all)),
  labels = c(
    "Promoter\nH3K27ac",
    "Promoter\nH3K4me1",
    "Fixed distal\nH3K27ac",
    "Fixed distal\nH3K4me1"
  ),
  las = 2,
  cex.axis = 0.85
)

axis(
  2,
  at = seq_len(nrow(support_all)),
  labels = rev(labels_all),
  las = 2,
  cex.axis = 0.8
)

abline(
  v = seq(
    1.5,
    ncol(support_all) - 0.5,
    by = 1
  ),
  col = "white"
)

abline(
  h = seq(
    1.5,
    nrow(support_all) - 0.5,
    by = 1
  ),
  col = "white"
)

mtext(
  "* = RNA HIGH priority; filled cell = same-direction Day1/Day6 plus BG reversal toward RPMI",
  side = 1,
  line = 10.5,
  cex = 0.7
)

dev.off()

# ============================================================
# 10. FIGURE 08B
#
# Day6 RNA memory effect among the 10 HIGH candidates,
# colored by corrected chromatin class.
# ============================================================

class_short <- ifelse(
  high$Any_corrected_dual_mark_support,
  "Dual-mark chromatin",
  ifelse(
    high$Any_corrected_chromatin_support,
    "Single-mark chromatin",
    "No directional chromatin"
  )
)

class_colors <- c(
  "Dual-mark chromatin" = "#173F5F",
  "Single-mark chromatin" = "#6C8EAD",
  "No directional chromatin" = "#D9D9D9"
)

bar_colors <- unname(
  class_colors[class_short]
)

png(
  file.path(
    figdir,
    "08B_HIGH_RNA_Day6_effect_with_chromatin_class.png"
  ),
  width = 2400,
  height = 1800,
  res = 220
)

par(
  mar = c(6, 10, 5, 3)
)

mids <- barplot(
  high$Day6_memory_log2FC,
  names.arg = high$SYMBOL,
  horiz = TRUE,
  las = 1,
  col = bar_colors,
  border = NA,
  xlab = "Day-6 LPS vs RPMI RNA log2FC",
  main =
    "HIGH-priority memory kinases: RNA effect and corrected chromatin support",
  cex.names = 0.85
)

abline(
  v = 0,
  lwd = 1
)

legend(
  "bottomright",
  legend = names(class_colors),
  fill = class_colors,
  border = NA,
  bty = "n",
  cex = 0.8
)

mtext(
  "Chromatin classes are directional/descriptive; no ChIP significance test is implied.",
  side = 1,
  line = 4,
  cex = 0.7
)

dev.off()

# ============================================================
# 11. FIGURE 08C
#
# For TOP candidates only:
# normalized directional effects for four chromatin units:
#   promoter H3K27ac
#   promoter H3K4me1
#   fixed distal H3K27ac
#   fixed distal H3K4me1
#
# Within each mark-region unit, Day1, Day6 and BG shift are
# scaled by the maximum absolute value for that unit.
# Therefore SIGNS/PATTERN should be interpreted, not magnitude
# across different marks.
# ============================================================

scale_triplet <- function(x) {

  x <- as.numeric(x)

  den <- max(
    abs(x),
    na.rm = TRUE
  )

  if (!is.finite(den) || den == 0) {
    return(
      rep(
        0,
        length(x)
      )
    )
  }

  x / den
}

if (nrow(top) > 0) {

  top_pidx <- match(
    top$ENSEMBL,
    promoter$ENSEMBL
  )

  top_didx <- match(
    top$ENSEMBL,
    distal$ENSEMBL
  )

  png(
    file.path(
      figdir,
      "08C_TOP_candidates_chromatin_direction.png"
    ),
    width = 2500,
    height = max(
      1800,
      650 * nrow(top)
    ),
    res = 220
  )

  par(
    mfrow = c(
      nrow(top),
      1
    ),
    mar = c(6, 5, 4, 2),
    oma = c(3, 0, 2, 0)
  )

  phase_colors <- c(
    "#4C78A8",
    "#1F3B5B",
    "#D97706"
  )

  for (i in seq_len(nrow(top))) {

    p <- promoter[top_pidx[i], , drop = FALSE]
    d <- distal[top_didx[i], , drop = FALSE]

    m <- cbind(
      Promoter_H3K27ac = scale_triplet(
        c(
          p$H3K27ac_Day1_delta,
          p$H3K27ac_Day6_delta,
          p$H3K27ac_BG_shift
        )
      ),

      Promoter_H3K4me1 = scale_triplet(
        c(
          p$H3K4me1_Day1_delta,
          p$H3K4me1_Day6_delta,
          p$H3K4me1_BG_shift
        )
      ),

      Fixed_Distal_H3K27ac = scale_triplet(
        c(
          d$H3K27ac_Day1_delta,
          d$H3K27ac_Day6_delta,
          d$H3K27ac_BG_shift
        )
      ),

      Fixed_Distal_H3K4me1 = scale_triplet(
        c(
          d$H3K4me1_Day1_delta,
          d$H3K4me1_Day6_delta,
          d$H3K4me1_BG_shift
        )
      )
    )

    rownames(m) <- c(
      "Day1 LPS-RPMI",
      "Day6 LPS-RPMI",
      "BG shift"
    )

    barplot(
      m,
      beside = TRUE,
      col = phase_colors,
      border = NA,
      ylim = c(-1.15, 1.15),
      names.arg = c(
        "Prom\nH3K27ac",
        "Prom\nH3K4me1",
        "Distal\nH3K27ac",
        "Distal\nH3K4me1"
      ),
      ylab = "Within-unit scaled effect",
      main = top$SYMBOL[i],
      cex.names = 0.78
    )

    abline(
      h = 0,
      lwd = 1
    )

    if (i == 1) {
      legend(
        "bottomleft",
        legend = rownames(m),
        fill = phase_colors,
        border = NA,
        bty = "n",
        cex = 0.72
      )
    }
  }

  mtext(
    "Interpret direction/pattern only; each mark-region unit is internally scaled.",
    side = 1,
    outer = TRUE,
    line = 1,
    cex = 0.75
  )

  dev.off()
}

# ============================================================
# 12. FIGURE 08D — FINAL PROFESSOR-FACING SUMMARY
#
# Table-style figure for the 10 HIGH RNA candidates.
# TOP integrated rows are highlighted.
# ============================================================

high_table <- high

high_table$Chromatin_label <- ifelse(
  high_table$Any_corrected_dual_mark_support,
  "DUAL",
  ifelse(
    high_table$Any_corrected_chromatin_support,
    "SINGLE",
    "NONE"
  )
)

high_table$Top_integrated <-
  high_table$Corrected_RNA_plus_chromatin_class ==
  "TOP_RNA_HIGH_PLUS_DUAL_CHROMATIN"

write.csv(
  high_table,
  file.path(
    outdir,
    "08_FINAL_professor_table.csv"
  ),
  row.names = FALSE
)

png(
  file.path(
    figdir,
    "08D_FINAL_RNA_chromatin_integration.png"
  ),
  width = 3000,
  height = 1900,
  res = 220
)

par(
  mar = c(1, 1, 4, 1)
)

plot.new()

plot.window(
  xlim = c(0, 1),
  ylim = c(
    0,
    nrow(high_table) + 2.4
  )
)

header_y <- nrow(high_table) + 1.5

x_gene <- 0.05
x_fc <- 0.25
x_rna <- 0.42
x_prom <- 0.58
x_dist <- 0.72
x_chr <- 0.88

text(
  0.5,
  nrow(high_table) + 2.15,
  "Final RNA + corrected chromatin integration",
  cex = 1.35,
  font = 2
)

text(
  x_gene,
  header_y,
  "Kinase",
  adj = 0,
  font = 2,
  cex = 0.9
)

text(
  x_fc,
  header_y,
  "Day6 RNA\nlog2FC",
  font = 2,
  cex = 0.82
)

text(
  x_rna,
  header_y,
  "RNA support\nscore",
  font = 2,
  cex = 0.82
)

text(
  x_prom,
  header_y,
  "Promoter\nmarks (0-2)",
  font = 2,
  cex = 0.82
)

text(
  x_dist,
  header_y,
  "Fixed distal\nmarks (0-2)",
  font = 2,
  cex = 0.82
)

text(
  x_chr,
  header_y,
  "Corrected\nchromatin",
  font = 2,
  cex = 0.82
)

abline(
  h = header_y - 0.55,
  col = "#777777"
)

for (i in seq_len(nrow(high_table))) {

  y <- nrow(high_table) - i + 1

  if (high_table$Top_integrated[i]) {

    rect(
      0.025,
      y - 0.42,
      0.975,
      y + 0.42,
      col = "#E8F0F7",
      border = NA
    )
  }

  this_font <- ifelse(
    high_table$Top_integrated[i],
    2,
    1
  )

  text(
    x_gene,
    y,
    high_table$SYMBOL[i],
    adj = 0,
    font = this_font,
    cex = 0.88
  )

  text(
    x_fc,
    y,
    sprintf(
      "%.2f",
      high_table$Day6_memory_log2FC[i]
    ),
    font = this_font,
    cex = 0.84
  )

  text(
    x_rna,
    y,
    high_table$Directional_support_score[i],
    font = this_font,
    cex = 0.84
  )

  text(
    x_prom,
    y,
    high_table$Promoter_mark_support_count[i],
    font = this_font,
    cex = 0.84
  )

  text(
    x_dist,
    y,
    high_table$Distal_mark_support_count[i],
    font = this_font,
    cex = 0.84
  )

  text(
    x_chr,
    y,
    high_table$Chromatin_label[i],
    font = this_font,
    cex = 0.82
  )

  segments(
    0.03,
    y - 0.48,
    0.97,
    y - 0.48,
    col = "#E3E3E3"
  )
}

mtext(
  "Highlighted rows = RNA HIGH plus dual-mark corrected chromatin support.",
  side = 1,
  line = -1.2,
  cex = 0.78
)

mtext(
  "Chromatin evidence is directional/descriptive and does not establish causal epigenetic inheritance.",
  side = 1,
  line = -2.2,
  cex = 0.68
)

dev.off()

# ============================================================
# 13. TEXT SUMMARY
# ============================================================

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
  "GSE85246 — Step 08 Final Corrected Chromatin Figures",
  "",
  paste("Memory kinases:", nrow(gene)),
  paste("RNA HIGH-priority kinases:", nrow(high)),
  paste("Any corrected chromatin support:", n_any),
  paste("Any corrected dual-mark support:", n_dual),
  paste("RNA HIGH with any corrected chromatin support:", n_high_any),
  paste("RNA HIGH with corrected dual-mark support:", n_high_dual),
  "",
  paste(
    "TOP RNA-HIGH + dual-chromatin candidates:",
    nrow(top)
  ),
  if (nrow(top) > 0) {
    paste(top$SYMBOL, collapse = ", ")
  } else {
    "None"
  },
  "",
  "Figure 08A: all 24 corrected chromatin support matrix",
  "Figure 08B: HIGH RNA Day6 effect with corrected chromatin class",
  "Figure 08C: TOP-candidate directional chromatin patterns",
  "Figure 08D: final professor-facing RNA + chromatin integration",
  "",
  paste(
    "Interpretation:",
    "The strongest integrated candidates combine HIGH RNA-level",
    "memory/rescue evidence with corrected dual-mark chromatin",
    "support at a predefined promoter or condition-blind selected",
    "distal locus."
  ),
  "",
  paste(
    "Caution:",
    "Chromatin evidence is directional/descriptive; Day1 and",
    "Day6 ChIP donor sets are not longitudinally paired, and",
    "these data do not establish causal epigenetic inheritance."
  )
)

writeLines(
  summary_lines,
  file.path(
    outdir,
    "08_final_chromatin_summary.txt"
  )
)

# ============================================================
# 14. FINAL CONSOLE REPORT
# ============================================================

cat("\n========================================\n")
cat("08 FINAL CHROMATIN FIGURES COMPLETE\n")
cat("========================================\n\n")

cat(
  "Memory kinases:",
  nrow(gene),
  "\n"
)

cat(
  "RNA HIGH-priority kinases:",
  nrow(high),
  "\n"
)

cat(
  "Any corrected chromatin support:",
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
  "RNA HIGH with any corrected chromatin support:",
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
  "\nTOP integrated candidates:",
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
  "\nMain professor-facing figure:\n",
  "figures/08_final_chromatin/",
  "08D_FINAL_RNA_chromatin_integration.png\n",
  sep = ""
)

cat(
  "\nIGV loci prepared:\n",
  "results/08_final_chromatin/",
  "08_TOP_candidate_IGV_loci_hg19.csv\n",
  sep = ""
)

cat(
  "\nNo new statistical testing was performed in Step 08.\n"
)

cat("\n========================================\n")
