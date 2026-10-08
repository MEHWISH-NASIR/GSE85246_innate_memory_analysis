
rm(list = ls())

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
})

setDTthreads(1)

cat("\n====================================================\n")
cat("21a — GSE262098 JAK3 BCG TRAJECTORY\n")
cat("PAIRED DONOR DESeq2 ANALYSIS\n")
cat("====================================================\n\n")

indir <- "macrophage_transcriptomics_validation/GSE262098/data/processed"
outdir <- "macrophage_transcriptomics_validation/GSE262098/results"
figdir <- "macrophage_transcriptomics_validation/GSE262098/figures"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

JAK3_ID <- "ENSG00000105639"

# ============================================================
# 1. EXACT SAMPLE MANIFEST
#
# exp1/exp2/exp3 are treated as the paired biological units.
# Only RPMI and BCG samples are included here.
# Later-added LPS samples are NOT mixed into this analysis.
# ============================================================

samples <- data.frame(
  sample = c(
    "exp1_d1_RPMI",
    "exp1_d1_BCG",
    "exp1_d6_RPMI",
    "exp1_d6_BCG",

    "exp2_d1_RPMI",
    "exp2_d1_BCG",
    "exp2_d6_RPMI",
    "exp2_d6_BCG",

    "exp3_d1_RPMI",
    "exp3_d1_BCG",
    "exp3_d6_RPMI",
    "exp3_d6_BCG"
  ),

  donor = factor(
    rep(
      c("exp1", "exp2", "exp3"),
      each = 4
    )
  ),

  time = factor(
    rep(
      c("D1", "D1", "D6", "D6"),
      times = 3
    ),
    levels = c("D1", "D6")
  ),

  condition = factor(
    rep(
      c("RPMI", "BCG", "RPMI", "BCG"),
      times = 3
    ),
    levels = c("RPMI", "BCG")
  ),

  file = c(
    "GSM8157504_exp1_d1_RPMI.counts.gz",
    "GSM8157505_exp1_d1_BCG.counts.gz",
    "GSM8157506_exp1_d6_RPMI.counts.gz",
    "GSM8157507_exp1_d6_BCG.counts.gz",

    "GSM8157509_exp2_d1_RPMI.counts.gz",
    "GSM8157510_exp2_d1_BCG.counts.gz",
    "GSM8157511_exp2_d6_RPMI.counts.gz",
    "GSM8157512_exp2_d6_BCG.counts.gz",

    "GSM8157514_exp3_d1_RPMI.counts.gz",
    "GSM8157515_exp3_d1_BCG.counts.gz",
    "GSM8157516_exp3_d6_RPMI.counts.gz",
    "GSM8157517_exp3_d6_BCG.counts.gz"
  ),

  stringsAsFactors = FALSE
)

samples$path <- file.path(
  indir,
  samples$file
)

if (!all(file.exists(samples$path))) {
  stop(
    paste(
      "Missing files:",
      paste(
        samples$path[!file.exists(samples$path)],
        collapse = "\n"
      )
    )
  )
}

rownames(samples) <- samples$sample

cat("Samples:", nrow(samples), "\n")
cat("Donors:", length(unique(samples$donor)), "\n")

# ============================================================
# 2. READ HTSEQ COUNTS
# ============================================================

count_list <- vector(
  "list",
  nrow(samples)
)

for (i in seq_len(nrow(samples))) {

  x <- fread(
    samples$path[i],
    header = FALSE,
    nThread = 1
  )

  if (ncol(x) != 2) {
    stop(
      "Expected two columns in ",
      samples$file[i]
    )
  }

  setnames(
    x,
    c("gene_id_version", "count")
  )

  # Remove Ensembl version suffix.
  x[
    ,
    gene_id :=
      sub("\\.[0-9]+$", "", gene_id_version)
  ]

  # Remove HTSeq technical summary rows if present.
  x <- x[
    !grepl("^__", gene_id)
  ]

  if (anyDuplicated(x$gene_id)) {
    stop(
      "Duplicated gene IDs after version stripping in ",
      samples$file[i]
    )
  }

  if (
    any(!is.finite(x$count)) ||
    any(x$count < 0) ||
    any(x$count != round(x$count))
  ) {
    stop(
      "Non-integer or invalid counts in ",
      samples$file[i]
    )
  }

  x <- x[
    ,
    .(
      gene_id,
      count = as.integer(count)
    )
  ]

  count_list[[i]] <- x
}

# ============================================================
# 3. VERIFY IDENTICAL GENE UNIVERSE
# ============================================================

reference_genes <- count_list[[1]]$gene_id

for (i in seq_along(count_list)) {

  if (!identical(
    reference_genes,
    count_list[[i]]$gene_id
  )) {
    stop(
      "Gene ordering/universe differs in file: ",
      samples$file[i]
    )
  }
}

counts <- do.call(
  cbind,
  lapply(
    count_list,
    function(z) z$count
  )
)

rownames(counts) <- reference_genes
colnames(counts) <- samples$sample

storage.mode(counts) <- "integer"

cat(
  "Count matrix:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "samples\n"
)

if (!JAK3_ID %in% rownames(counts)) {
  stop("JAK3 not found in count matrix.")
}

# ============================================================
# 4. SAVE RAW JAK3 COUNTS
# ============================================================

jak3_raw <- data.frame(
  samples[
    ,
    c(
      "sample",
      "donor",
      "time",
      "condition"
    )
  ],
  JAK3_raw_count =
    as.integer(
      counts[JAK3_ID, ]
    ),
  stringsAsFactors = FALSE
)

write.csv(
  jak3_raw,
  file.path(
    outdir,
    "21a_GSE262098_JAK3_raw_counts.csv"
  ),
  row.names = FALSE
)

cat("\n===== RAW JAK3 COUNTS =====\n")
print(jak3_raw, row.names = FALSE)

# ============================================================
# 5. FILTER LOW-COUNT GENES
#
# Require >=10 counts in >=3 samples.
# The filter is defined independently of JAK3 significance.
# ============================================================

keep <- rowSums(
  counts >= 10
) >= 3

cat("\nGenes before filtering:", nrow(counts), "\n")
cat("Genes retained:", sum(keep), "\n")

if (!keep[JAK3_ID]) {
  stop("JAK3 failed the predefined count filter.")
}

counts_f <- counts[
  keep,
  ,
  drop = FALSE
]

# ============================================================
# 6. DESEQ2 INTERACTION MODEL
#
# Reference:
# time      = D1
# condition = RPMI
#
# condition_BCG_vs_RPMI
#     = BCG effect at Day1
#
# timeD6.conditionBCG
#     = difference in BCG effect at Day6 vs Day1
# ============================================================

coldata <- samples[
  ,
  c(
    "donor",
    "time",
    "condition"
  )
]

dds <- DESeqDataSetFromMatrix(
  countData = counts_f,
  colData = coldata,
  design =
    ~ donor +
      time +
      condition +
      time:condition
)

dds <- DESeq(
  dds,
  quiet = FALSE
)

cat("\n===== DESEQ2 COEFFICIENTS =====\n")
print(resultsNames(dds))

rn <- resultsNames(dds)

# Identify coefficient names robustly.
bcg_main <- grep(
  "^condition.*BCG.*RPMI",
  rn,
  value = TRUE
)

interaction <- grep(
  "timeD6.*conditionBCG|conditionBCG.*timeD6",
  rn,
  value = TRUE
)

if (length(bcg_main) != 1) {
  stop(
    "Could not uniquely identify Day1 BCG coefficient.\n",
    paste(rn, collapse = "\n")
  )
}

if (length(interaction) != 1) {
  stop(
    "Could not uniquely identify Day6 interaction coefficient.\n",
    paste(rn, collapse = "\n")
  )
}

cat("\nDay1 BCG coefficient:", bcg_main, "\n")
cat("Trajectory interaction:", interaction, "\n")

# ============================================================
# 7. THREE CORE TESTS
# ============================================================

# Day1 BCG vs Day1 RPMI
res_D1 <- results(
  dds,
  name = bcg_main,
  alpha = 0.05
)

# Direct trajectory:
# (D6 BCG - D6 RPMI) -
# (D1 BCG - D1 RPMI)
res_interaction <- results(
  dds,
  name = interaction,
  alpha = 0.05
)

# Day6 BCG vs Day6 RPMI =
# Day1 main BCG effect + interaction.
res_D6 <- results(
  dds,
  contrast = list(
    c(
      bcg_main,
      interaction
    )
  ),
  alpha = 0.05
)

# ============================================================
# 8. SAVE GENOME-WIDE TABLES
# ============================================================

result_to_dt <- function(
  res,
  label
) {

  z <- as.data.table(
    as.data.frame(res),
    keep.rownames = "gene_id"
  )

  z[, contrast := label]

  z
}

all_D1 <- result_to_dt(
  res_D1,
  "Day1_BCG_vs_RPMI"
)

all_D6 <- result_to_dt(
  res_D6,
  "Day6_BCG_vs_RPMI"
)

all_int <- result_to_dt(
  res_interaction,
  "Trajectory_D6_vs_D1"
)

fwrite(
  all_D1,
  file.path(
    outdir,
    "21a_GSE262098_Day1_BCG_vs_RPMI_DESeq2.csv"
  )
)

fwrite(
  all_D6,
  file.path(
    outdir,
    "21a_GSE262098_Day6_BCG_vs_RPMI_DESeq2.csv"
  )
)

fwrite(
  all_int,
  file.path(
    outdir,
    "21a_GSE262098_D6_vs_D1_interaction_DESeq2.csv"
  )
)

# ============================================================
# 9. EXTRACT JAK3
# ============================================================

extract_jak3 <- function(
  tab
) {

  tab[
    gene_id == JAK3_ID,
    .(
      contrast,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj
    )
  ]
}

jak3_stats <- rbind(
  extract_jak3(all_D1),
  extract_jak3(all_D6),
  extract_jak3(all_int)
)

if (nrow(jak3_stats) != 3) {
  stop("Expected exactly three JAK3 result rows.")
}

fwrite(
  jak3_stats,
  file.path(
    outdir,
    "21a_GSE262098_JAK3_DESeq2_summary.csv"
  )
)

cat("\n========================================\n")
cat("JAK3 DESEQ2 RESULTS\n")
cat("========================================\n")
print(jak3_stats)

# ============================================================
# 10. NORMALIZED JAK3 COUNTS
# ============================================================

norm <- counts(
  dds,
  normalized = TRUE
)

jak3_norm <- data.table(
  sample = colnames(norm),
  donor = samples[colnames(norm), "donor"],
  time = samples[colnames(norm), "time"],
  condition =
    samples[colnames(norm), "condition"],
  normalized_count =
    as.numeric(
      norm[JAK3_ID, ]
    )
)

fwrite(
  jak3_norm,
  file.path(
    outdir,
    "21a_GSE262098_JAK3_normalized_counts.csv"
  )
)

cat("\n===== NORMALIZED JAK3 COUNTS =====\n")
print(jak3_norm)

# ============================================================
# 11. DONOR-WISE EFFECTS FOR VISUALIZATION
#
# Use log2(normalized_count + 1)
# only for visualization / descriptive paired differences.
# Formal inference remains DESeq2 on raw counts.
# ============================================================

jak3_norm[
  ,
  log2_norm :=
    log2(normalized_count + 1)
]

wide <- dcast(
  jak3_norm,
  donor + time ~ condition,
  value.var = "log2_norm"
)

wide[
  ,
  BCG_minus_RPMI :=
    BCG - RPMI
]

fwrite(
  wide,
  file.path(
    outdir,
    "21a_GSE262098_JAK3_donor_effects.csv"
  )
)

cat("\n===== DONOR-WISE JAK3 EFFECTS =====\n")
print(wide)

# ============================================================
# 12. DIRECTION CONSISTENCY
# ============================================================

direction <- wide[
  ,
  .(
    mean_descriptive_effect =
      mean(BCG_minus_RPMI),

    median_descriptive_effect =
      median(BCG_minus_RPMI),

    positive_donors =
      sum(BCG_minus_RPMI > 0),

    negative_donors =
      sum(BCG_minus_RPMI < 0)
  ),
  by = time
]

direction[
  ,
  direction_consistency :=
    pmax(
      positive_donors,
      negative_donors
    ) / 3
]

fwrite(
  direction,
  file.path(
    outdir,
    "21a_GSE262098_JAK3_direction_consistency.csv"
  )
)

cat("\n===== DIRECTION CONSISTENCY =====\n")
print(direction)

# ============================================================
# 13. FIGURE — PAIRED JAK3 EXPRESSION
# ============================================================

plotdat <- copy(jak3_norm)

plotdat$time <- factor(
  plotdat$time,
  levels = c("D1", "D6"),
  labels = c(
    "24 h",
    "Day 6 post-washout"
  )
)

plotdat$condition <- factor(
  plotdat$condition,
  levels = c("RPMI", "BCG")
)

p1 <- ggplot(
  plotdat,
  aes(
    x = condition,
    y = log2_norm,
    group = donor
  )
) +
  geom_line(
    alpha = 0.65
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ time,
    nrow = 1
  ) +
  labs(
    title =
      "GSE262098: JAK3 expression during BCG training and post-washout state",
    subtitle =
      "Paired measurements across three independent donors",
    x = NULL,
    y = "log2(DESeq2 normalized count + 1)"
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "21a_GSE262098_JAK3_paired_expression.png"
  ),
  p1,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "21a_GSE262098_JAK3_paired_expression.pdf"
  ),
  p1,
  width = 8,
  height = 5
)

# ============================================================
# 14. FIGURE — EFFECT SIZES + 95% CI
# ============================================================

forest <- copy(jak3_stats)

forest <- forest[
  contrast %in%
    c(
      "Day1_BCG_vs_RPMI",
      "Day6_BCG_vs_RPMI"
    )
]

forest[
  ,
  CI_low :=
    log2FoldChange -
    1.96 * lfcSE
]

forest[
  ,
  CI_high :=
    log2FoldChange +
    1.96 * lfcSE
]

forest[
  ,
  time_label :=
    fifelse(
      contrast == "Day1_BCG_vs_RPMI",
      "24 h",
      "Day 6 post-washout"
    )
]

forest$time_label <- factor(
  forest$time_label,
  levels = c(
    "24 h",
    "Day 6 post-washout"
  )
)

p2 <- ggplot(
  forest,
  aes(
    x = log2FoldChange,
    y = time_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_segment(
    aes(
      x = CI_low,
      xend = CI_high,
      y = time_label,
      yend = time_label
    ),
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    title =
      "GSE262098: JAK3 BCG treatment effects",
    subtitle =
      "DESeq2 log2 fold changes with approximate 95% confidence intervals",
    x = "BCG vs RPMI log2 fold change",
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "21a_GSE262098_JAK3_effect_sizes.png"
  ),
  p2,
  width = 7,
  height = 4.5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "21a_GSE262098_JAK3_effect_sizes.pdf"
  ),
  p2,
  width = 7,
  height = 4.5
)

# ============================================================
# 15. ANALYSIS METADATA
# ============================================================

metadata <- data.frame(
  dataset = "GSE262098",
  species = "Homo sapiens",
  biological_units = 3,
  biological_unit_label = "exp1/exp2/exp3",
  assay = "RNA-seq",
  count_source =
    "Deposited unnormalized HTSeq gene counts",
  genome_annotation =
    "Ensembl GRCh37v70 / hg19",
  primary_stimulus = "BCG",
  initial_exposure = "24 h",
  late_state = "Day6 post-washout",
  design =
    "~ donor + time + condition + time:condition",
  primary_tests =
    paste(
      "Day1 BCG vs RPMI;",
      "Day6 BCG vs RPMI;",
      "Day6-vs-Day1 interaction"
    ),
  low_count_filter =
    ">=10 counts in >=3 of 12 samples",
  stringsAsFactors = FALSE
)

write.csv(
  metadata,
  file.path(
    outdir,
    "21a_GSE262098_analysis_metadata.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("21a COMPLETE\n")
cat("========================================\n")

