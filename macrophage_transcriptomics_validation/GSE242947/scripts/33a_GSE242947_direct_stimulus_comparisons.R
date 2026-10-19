suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

root <- "macrophage_transcriptomics_validation/GSE242947"
input <- file.path(root, "data/GSE242947_Cheng_RNAseqCounts.csv.gz")
outdir <- file.path(root, "results")

sample_cols <- c(
  "Untrained.0h.R1", "Untrained.0h.R2",
  "Zymo.0h.R1", "Zymo.0h.R2",
  "ETSU.0h.R1", "ETSU.0h.R2",
  "Depleted.0h.R1", "Depleted.0h.R2"
)

condition <- factor(
  rep(c("Untrained", "Zymosan", "ETSU_BG", "Depleted"), each=2),
  levels=c("Untrained", "Zymosan", "ETSU_BG", "Depleted")
)

x <- fread(input, nThread=1)

stopifnot("Geneid" %in% names(x))
stopifnot(all(sample_cols %in% names(x)))

# Match the original 28a preprocessing:
# remove Ensembl version suffixes and sum duplicate gene IDs.
x[, gene_id := sub("\\..*$", "", Geneid)]

count_dt <- x[
  ,
  lapply(.SD, sum),
  by=gene_id,
  .SDcols=sample_cols
]

mat <- as.matrix(count_dt[, ..sample_cols])
storage.mode(mat) <- "integer"
rownames(mat) <- count_dt$gene_id
colnames(mat) <- sample_cols

stopifnot(
  !anyNA(mat),
  all(mat >= 0),
  !anyDuplicated(rownames(mat))
)

keep <- rowSums(mat >= 10) >= 2
mat <- mat[keep, , drop=FALSE]

jak3 <- "ENSG00000105639"
stopifnot(jak3 %in% rownames(mat))

meta <- data.frame(condition=condition, row.names=sample_cols)

dds <- DESeqDataSetFromMatrix(
  countData=mat,
  colData=meta,
  design=~ condition
)

dds <- DESeq(dds)

comparisons <- list(
  Zymosan_vs_ETSU_BG=c("condition","Zymosan","ETSU_BG"),
  Zymosan_vs_Depleted=c("condition","Zymosan","Depleted"),
  Depleted_vs_ETSU_BG=c("condition","Depleted","ETSU_BG")
)

output <- list()

for (nm in names(comparisons)) {

  res <- results(
    dds,
    contrast=comparisons[[nm]],
    alpha=0.05
  )

  j <- as.data.frame(res[jak3, , drop=FALSE])

  output[[nm]] <- data.frame(
    Dataset="GSE242947",
    Gene="JAK3",
    Comparison=nm,
    log2FC=j$log2FoldChange,
    lfcSE=j$lfcSE,
    CI_low=j$log2FoldChange-1.96*j$lfcSE,
    CI_high=j$log2FoldChange+1.96*j$lfcSE,
    P_value=j$pvalue,
    FDR=j$padj
  )
}

final <- do.call(rbind, output)
rownames(final) <- NULL

outfile <- file.path(
  outdir,
  "33a_GSE242947_JAK3_direct_stimulus_comparisons.csv"
)

if (file.exists(outfile)) {
  stop("Output already exists; refusing to overwrite.")
}

fwrite(final, outfile)

cat("\n===== DIRECT JAK3 COMPARISONS =====\n")
print(final, row.names=FALSE)
cat("\nSaved:", outfile, "\n")
