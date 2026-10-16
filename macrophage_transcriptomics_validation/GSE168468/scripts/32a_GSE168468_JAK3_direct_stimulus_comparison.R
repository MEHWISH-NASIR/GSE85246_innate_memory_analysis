# 32a — GSE168468 direct Day6 stimulus contrasts (new outputs only)
# Run from the GSE85246_innate_memory_analysis project root.
# Reuses official MMSEQ normalization and the donor-adjusted limma model.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})
setDTthreads(1)

base <- 'macrophage_transcriptomics_validation/GSE168468'
indir <- file.path(base, 'data', 'processed')
outdir <- file.path(base, 'results')
helper <- file.path(base, 'scripts', 'mmseq.R')
stopifnot(dir.exists(indir), file.exists(helper))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
source(helper)

# Same 20 resting Day6 samples and order as the existing 20c analysis.
accessions <- c(5141114:5141118, 5141119:5141123,
                5141124:5141128, 5141129:5141133)
suffix <- c('A1_D6_RPMI', 'A2_d6R', 'A3_d6R', 'A4_d6R', 'A5_d6R',
            paste0('A',1:5,'_d6L'),
            'A1_D6_Den', paste0('A',2:5,'_d6D'),
            'A1_D6_Bul', paste0('A',2:5,'_d6B'))
stopifnot(length(accessions) == 20L, length(suffix) == 20L)
filenames <- paste0('GSM',accessions,'_',suffix,'.gene.mmseq.txt.gz')
samples <- data.frame(
  donor = rep(paste0('A', 1:5), 4),
  condition = factor(rep(c('RPMI','LPS','BCG_Denmark','BCG_Bulgaria'),each=5),
                     levels=c('RPMI','LPS','BCG_Denmark','BCG_Bulgaria')),
  file = filenames,
  stringsAsFactors=FALSE
)
samples$path <- file.path(indir,samples$file)
samples$sample_name <- paste(samples$donor,samples$condition,sep='_')
if (!all(file.exists(samples$path))) {
  stop('Missing processed samples:\n',paste(samples$path[!file.exists(samples$path)],collapse='\n'))
}
cat('Using',nrow(samples),'Day6 resting samples from',length(unique(samples$donor)),'matched donors.\n')

mm <- readmmseq(mmseq_files=samples$path,
               sample_names=samples$sample_name,
               normalize=TRUE)
expr <- mm$log_mu
if (ncol(expr) != nrow(samples)) stop('Expression/sample mismatch')
if (!identical(colnames(expr),samples$sample_name)) {
  stop('MMSEQ columns not identical to sample manifest ordering. Check before analysis.')
}
jak3_id <- 'ENSG00000105639'
keep <- rowSums(mm$observed == 1,na.rm=TRUE) >= 5
if (!isTRUE(keep[jak3_id])) stop('JAK3 did not pass original observed filter')
expr_f <- expr[keep,,drop=FALSE]
donor <- factor(samples$donor)
condition <- samples$condition
design <- model.matrix(~ donor + condition)
fit <- lmFit(expr_f,design)

# Explicit coefficient-weight matrix avoids model-name parsing issues.
coefs <- colnames(design)
needed <- c('conditionLPS','conditionBCG_Denmark','conditionBCG_Bulgaria')
if (!all(needed %in% coefs)) stop('Expected condition coefficients absent')
C <- matrix(0,nrow=ncol(design),ncol=4,
            dimnames=list(coefs,c('Denmark_vs_LPS','Bulgaria_vs_LPS',
                                  'BCG_average_vs_LPS','Denmark_vs_Bulgaria')))
C['conditionBCG_Denmark','Denmark_vs_LPS'] <- 1
C['conditionLPS','Denmark_vs_LPS'] <- -1
C['conditionBCG_Bulgaria','Bulgaria_vs_LPS'] <- 1
C['conditionLPS','Bulgaria_vs_LPS'] <- -1
C['conditionBCG_Denmark','BCG_average_vs_LPS'] <- 0.5
C['conditionBCG_Bulgaria','BCG_average_vs_LPS'] <- 0.5
C['conditionLPS','BCG_average_vs_LPS'] <- -1
C['conditionBCG_Denmark','Denmark_vs_Bulgaria'] <- 1
C['conditionBCG_Bulgaria','Denmark_vs_Bulgaria'] <- -1

fit2 <- contrasts.fit(fit,C)
fit2 <- eBayes(fit2,robust=TRUE)

# BH adjustment separately across all tested genes, within each contrast.
result_list <- lapply(colnames(C),function(nm) {
  tt <- topTable(fit2,coef=nm,number=Inf,sort.by='none',
                 adjust.method='BH',confint=TRUE)
  tt$feature_id <- rownames(tt)
  tt$contrast <- nm
  as.data.table(tt)
})
all_results <- rbindlist(result_list,fill=TRUE)
jak3 <- all_results[feature_id == jak3_id]
if (nrow(jak3) != 4L) stop('Expected four JAK3 contrast results')

# The full matrix preserves a reproducible genome-wide FDR calculation.
all_path <- file.path(outdir,'32a_GSE168468_Day6_direct_contrasts_genomewide.csv')
jak_path <- file.path(outdir,'32a_GSE168468_Day6_JAK3_direct_contrasts.csv')
fwrite(all_results,all_path)
fwrite(jak3,jak_path)
cat('\n=== JAK3 DIRECT CONTRASTS (Day6 resting, paired donor limma) ===\n')
print(jak3[,.(contrast,logFC,CI.L,CI.R,P.Value,adj.P.Val)])
cat('\nSaved:\n',all_path,'\n',jak_path,'\n',sep='')
cat('Notes: reported logFC follows the MMSEQ log_mu scale; check log base before labelling log2FC.\n')
cat('BH FDR is gene-wise within each contrast; four-contrast multiplicity is not separately adjusted.\n')
