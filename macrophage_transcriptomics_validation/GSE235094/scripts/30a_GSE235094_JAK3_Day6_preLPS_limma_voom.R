
rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
})

setDTthreads(1)

cat("\n====================================================\n")
cat("30a — GSE235094 JAK3 DAY6 PRE-LPS ANALYSIS\n")
cat("iSARS-CoV-2 trained vs untrained monocytes\n")
cat("LIMMA-VOOM\n")
cat("====================================================\n\n")

infile <- "macrophage_transcriptomics_validation/GSE235094/data/GSE235094_Raw_counts_all_samples.txt.gz"
outdir <- "macrophage_transcriptomics_validation/GSE235094/results"
figdir <- "macrophage_transcriptomics_validation/GSE235094/figures"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

JAK3_ID <- "ENSG00000105639"

# ============================================================
# 1. READ MATRIX AS CHARACTER
# ============================================================

x <- fread(
  infile,
  sep = "\t",
  header = TRUE,
  colClasses = "character",
  nThread = 1
)

x[, gene_id := sub("\\..*$", "", Geneid)]

# Clean public Day6 pre-LPS samples only
untrained <- c(
  "Sample6",
  "Sample13",
  "Sample20",
  "Sample21",
  "Sample22"
)

trained <- c(
  "Sample3",
  "Sample7",
  "Sample10",
  "Sample14",
  "Sample16",
  "Sample26"
)

sample_cols <- c(
  untrained,
  trained
)

stopifnot(all(sample_cols %in% names(x)))

# ============================================================
# 2. NUMERIC PARSING
#
# Matrix contains both decimal commas and decimal points.
# Convert decimal commas to decimal points before numeric coercion.
# ============================================================

for (cc in sample_cols) {
  x[
    ,
    (cc) := as.numeric(
      gsub(",", ".", get(cc), fixed = TRUE)
    )
  ]
}

if (anyNA(x[, ..sample_cols])) {
  stop("NA introduced during numeric conversion.")
}

# ============================================================
# 3. COLLAPSE DUPLICATE STABLE ENSEMBL IDs
# ============================================================

count_dt <- x[
  ,
  lapply(.SD, sum),
  by = gene_id,
  .SDcols = sample_cols
]

counts <- as.matrix(
  count_dt[, ..sample_cols]
)

rownames(counts) <- count_dt$gene_id
storage.mode(counts) <- "numeric"

if (anyDuplicated(rownames(counts))) {
  stop("Duplicate gene IDs remain.")
}

cat(
  "Count matrix:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "samples\n"
)

# ============================================================
# 4. SAMPLE METADATA
# ============================================================

group <- factor(
  c(
    rep("Unt_Med", length(untrained)),
    rep("Trn_S_Med", length(trained))
  ),
  levels = c(
    "Unt_Med",
    "Trn_S_Med"
  )
)

meta <- data.frame(
  sample = sample_cols,
  group = group,
  row.names = sample_cols
)

cat("\n===== SAMPLE METADATA =====\n")
print(meta)

# ============================================================
# 5. FILTER + TMM NORMALIZATION
# ============================================================

dge <- DGEList(
  counts = counts,
  group = group
)

keep <- filterByExpr(
  dge,
  group = group
)

cat("\nGenes before filtering:", nrow(dge), "\n")
cat("Genes retained:", sum(keep), "\n")

dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

if (!JAK3_ID %in% rownames(dge)) {
  stop("JAK3 failed expression filter.")
}

# ============================================================
# 6. LIMMA-VOOM
# ============================================================

design <- model.matrix(
  ~ group
)

cat("\n===== DESIGN MATRIX =====\n")
print(design)

v <- voom(
  dge,
  design,
  plot = FALSE
)

fit <- lmFit(
  v,
  design
)

fit <- eBayes(
  fit,
  robust = TRUE
)

coef_name <- "groupTrn_S_Med"

tt <- topTable(
  fit,
  coef = coef_name,
  number = Inf,
  sort.by = "none",
  adjust.method = "BH",
  confint = TRUE
)

tt$gene_id <- rownames(tt)

fwrite(
  as.data.table(tt),
  file.path(
    outdir,
    "30a_GSE235094_Trn_S_Med_vs_Unt_Med_genomewide_limma_voom.csv"
  )
)

# ============================================================
# 7. JAK3 FORMAL RESULT
# ============================================================

j <- tt[
  tt$gene_id == JAK3_ID,
  ,
  drop = FALSE
]

if (nrow(j) != 1) {
  stop("Expected exactly one JAK3 row.")
}

jak3_result <- data.table(
  dataset = "GSE235094",
  gene = "JAK3",
  gene_id = JAK3_ID,
  contrast = "Day6_Trn_S_Med_vs_Unt_Med",
  logFC = j$logFC,
  CI_low = j$CI.L,
  CI_high = j$CI.R,
  t = j$t,
  P_value = j$P.Value,
  FDR = j$adj.P.Val
)

fwrite(
  jak3_result,
  file.path(
    outdir,
    "30a_GSE235094_JAK3_formal_result.csv"
  )
)

cat("\n========================================\n")
cat("JAK3 FORMAL RESULT\n")
cat("========================================\n")
print(jak3_result)

# ============================================================
# 8. NORMALIZED JAK3 EXPRESSION
# ============================================================

jak3_voom <- data.table(
  sample = sample_cols,
  group = group,
  voom_log2CPM = as.numeric(
    v$E[JAK3_ID, ]
  )
)

fwrite(
  jak3_voom,
  file.path(
    outdir,
    "30a_GSE235094_JAK3_voom_values.csv"
  )
)

cat("\n===== JAK3 NORMALIZED VALUES =====\n")
print(jak3_voom)

# ============================================================
# 9. GROUP SUMMARY
# ============================================================

summary <- jak3_voom[
  ,
  .(
    n = .N,
    mean = mean(voom_log2CPM),
    median = median(voom_log2CPM),
    sd = sd(voom_log2CPM)
  ),
  by = group
]

fwrite(
  summary,
  file.path(
    outdir,
    "30a_GSE235094_JAK3_group_summary.csv"
  )
)

cat("\n===== GROUP SUMMARY =====\n")
print(summary)

# ============================================================
# 10. FIGURE
# ============================================================

p <- ggplot(
  jak3_voom,
  aes(
    x = group,
    y = voom_log2CPM
  )
) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 3
  ) +
  labs(
    title =
      "GSE235094: JAK3 after iSARS-CoV-2 training",
    subtitle =
      "Day6 before secondary LPS stimulation",
    x = NULL,
    y = "JAK3 voom log2 CPM"
  ) +
  theme_bw(
    base_size = 11
  )

ggsave(
  file.path(
    figdir,
    "30a_GSE235094_JAK3_Day6_preLPS.png"
  ),
  p,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(
    figdir,
    "30a_GSE235094_JAK3_Day6_preLPS.pdf"
  ),
  p,
  width = 6,
  height = 5
)

# ============================================================
# 11. METADATA
# ============================================================

metadata <- data.frame(
  dataset = "GSE235094",
  species = "Homo sapiens",
  system = "Primary human monocytes",
  training = "Inactivated SARS-CoV-2",
  training_duration = "24 h",
  rest_period = "5 days",
  analysis_state = "Day6 before secondary LPS",
  untrained_samples = length(untrained),
  trained_samples = length(trained),
  excluded_samples =
    "Sample5 unmapped; Sample8 metadata inconsistent",
  expression_source =
    "Deposited fractional gene-level estimated counts",
  method =
    "TMM normalization + limma-voom",
  stringsAsFactors = FALSE
)

write.csv(
  metadata,
  file.path(
    outdir,
    "30a_GSE235094_analysis_metadata.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("30a COMPLETE\n")
cat("========================================\n")

