# 32b — GSE168468 five-gene comparator panel
# Uses previously fitted genome-wide results; does not refit the original model.
# Run from repository root: Rscript macrophage_transcriptomics_validation/GSE168468/scripts/32b_GSE168468_five_gene_comparator_panel.R
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
root <- 'macrophage_transcriptomics_validation/GSE168468'
res <- file.path(root, 'results')
fig <- file.path(root, 'figures')
dir.create(fig, recursive=TRUE, showWarnings=FALSE)
original_path <- file.path(res,'20c_GSE168468_Day6_genomewide_limma.csv')
direct_path <- file.path(res,'32a_GSE168468_Day6_direct_contrasts_genomewide.csv')
if (!file.exists(original_path) || !file.exists(direct_path)) {
  stop('Missing an existing genome-wide CSV. Complete/check analyses 20c and 32a first.')
}
ids <- c(JAK3='ENSG00000105639', HCK='ENSG00000101336',
         LYN='ENSG00000254087', IRAK2='ENSG00000134070',
         EPHB2='ENSG00000133216')
required <- c('feature_id','contrast','logFC','CI.L','CI.R','P.Value','adj.P.Val')
extract_panel <- function(path, expected_contrasts) {
  x <- fread(path)
  if (!all(required %in% names(x))) stop('Missing required result columns: ',path)
  x <- x[feature_id %in% unname(ids) & contrast %in% expected_contrasts]
  x[, gene := names(ids)[match(feature_id, unname(ids))]]
  if (nrow(x) != length(ids)*length(expected_contrasts)) {
    stop('Some gene x contrast combinations are missing in: ',path,
         ' (observed ',nrow(x),', expected ',length(ids)*length(expected_contrasts),')')
  }
  if (anyDuplicated(x[,.(gene,contrast)])) stop('Duplicate gene/contrast combinations: ',path)
  setcolorder(x,c('gene','feature_id','contrast',setdiff(names(x),c('gene','feature_id','contrast'))))
  x[]
}
original <- extract_panel(original_path,c('Day6_LPS_vs_RPMI',
                                        'Day6_BCG_Denmark_vs_RPMI',
                                        'Day6_BCG_Bulgaria_vs_RPMI'))
direct <- extract_panel(direct_path,c('Denmark_vs_LPS','Bulgaria_vs_LPS',
                                    'BCG_average_vs_LPS','Denmark_vs_Bulgaria'))
orig_out <- file.path(res,'32b_GSE168468_five_gene_vs_RPMI.csv')
dir_out <- file.path(res,'32b_GSE168468_five_gene_direct_contrasts.csv')
fwrite(original,orig_out)
fwrite(direct,dir_out)
# Plot individual gene effects; do not standardize rows as it would hide effect magnitudes.
plot_data <- copy(original)
plot_data[,gene:=factor(gene,levels=rev(names(ids)))]
plot_data[,contrast:=factor(contrast,levels=c('Day6_LPS_vs_RPMI',
                                            'Day6_BCG_Denmark_vs_RPMI',
                                            'Day6_BCG_Bulgaria_vs_RPMI'),
                             labels=c('LPS','BCG Denmark','BCG Bulgaria'))]
p <- ggplot(plot_data,aes(x=contrast,y=gene,fill=logFC)) +
  geom_tile(color='white',linewidth=0.8) +
  geom_text(aes(label=sprintf('%.2f',logFC)),size=4) +
  scale_fill_gradient2(low='#3A78A5',mid='white',high='#B54D42',midpoint=0,
                       name='MMSEQ log_mu\neffect') +
  labs(title='GSE168468: Day 6 post-washout kinase responses',
       subtitle='5 donor-paired samples per condition; contrasts against RPMI',
       x=NULL,y=NULL,
       caption='Numbers show fitted effects, not statistical significance. FDR in accompanying CSV.') +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(),plot.title=element_text(face='bold'),
        axis.text.x=element_text(angle=15,hjust=1))
png <- file.path(fig,'32b_GSE168468_five_gene_effect_heatmap.png')
pdf <- file.path(fig,'32b_GSE168468_five_gene_effect_heatmap.pdf')
ggsave(png,p,width=8,height=4.5,dpi=300)
ggsave(pdf,p,width=8,height=4.5)
cat('\n=== FIVE GENES: DAY6 VS RPMI ===\n')
print(original[,.(gene,contrast,logFC,CI.L,CI.R,P.Value,adj.P.Val)])
cat('\n=== FIVE GENES: DIRECT STIMULUS CONTRASTS ===\n')
print(direct[,.(gene,contrast,logFC,CI.L,CI.R,P.Value,adj.P.Val)])
cat('\nSaved new results:\n',orig_out,'\n',dir_out,'\n',png,'\n',pdf,'\n',sep='')
cat('\nNote: MMSEQ log_mu effects have not been verified as base-2 log fold changes.\n')
cat('FDR is computed genome-wide per contrast. Cross-gene effect differences are descriptive, not formal gene-by-gene significance tests.\n')
