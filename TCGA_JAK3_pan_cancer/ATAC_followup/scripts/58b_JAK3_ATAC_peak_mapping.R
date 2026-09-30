# ============================================================
# 58b — Map TCGA pan-cancer ATAC peaks around JAK3
#
# Genome build: GRCh38 / hg38
#
# JAK3:
# chr19:17,824,782-17,847,982
# strand: minus
# TSS: 17,847,982
#
# Regions:
# promoter = TSS +/- 2 kb
# gene body = annotated JAK3 interval
# nearby = gene interval +/- 50 kb
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

peak_file <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/data/TCGA-ATAC_PanCancer_PeakSet.txt"

outdir <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# JAK3 GRCh38 coordinates
# ------------------------------------------------------------

chr <- "chr19"

gene_start <- 17824782
gene_end   <- 17847982

# JAK3 is on minus strand
tss <- gene_end

promoter_start <- tss - 2000
promoter_end   <- tss + 2000

nearby_start <- gene_start - 50000
nearby_end   <- gene_end + 50000

# ------------------------------------------------------------
# Read pan-cancer peak set
# ------------------------------------------------------------

cat("\n========================================\n")
cat("JAK3 ATAC PEAK MAPPING\n")
cat("========================================\n\n")

cat("Reading peak set...\n")

peaks <- read.delim(
  peak_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Total pan-cancer peaks:", nrow(peaks), "\n")

cat(
  "Columns:",
  paste(names(peaks), collapse = ", "),
  "\n\n"
)

# ------------------------------------------------------------
# Basic validation
# ------------------------------------------------------------

required <- c(
  "seqnames",
  "start",
  "end",
  "name",
  "score",
  "annotation",
  "percentGC"
)

missing_cols <- setdiff(
  required,
  names(peaks)
)

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# Restrict to chr19 first
# ------------------------------------------------------------

chr19 <- peaks[
  peaks$seqnames == chr,
]

cat("chr19 peaks:", nrow(chr19), "\n\n")

# ------------------------------------------------------------
# Overlap helper
# ------------------------------------------------------------

overlaps_region <- function(
  start,
  end,
  region_start,
  region_end
) {

  start <= region_end &
  end >= region_start
}

# ------------------------------------------------------------
# Identify JAK3-related peaks
# ------------------------------------------------------------

chr19$Overlap_promoter <- overlaps_region(
  chr19$start,
  chr19$end,
  promoter_start,
  promoter_end
)

chr19$Overlap_gene_body <- overlaps_region(
  chr19$start,
  chr19$end,
  gene_start,
  gene_end
)

chr19$Within_JAK3_50kb <- overlaps_region(
  chr19$start,
  chr19$end,
  nearby_start,
  nearby_end
)

jak3_peaks <- chr19[
  chr19$Within_JAK3_50kb,
]

# ------------------------------------------------------------
# Assign mutually exclusive region category
#
# Priority:
# promoter > gene body > nearby
# ------------------------------------------------------------

jak3_peaks$JAK3_region <- "Nearby_50kb"

jak3_peaks$JAK3_region[
  jak3_peaks$Overlap_gene_body
] <- "Gene_body"

jak3_peaks$JAK3_region[
  jak3_peaks$Overlap_promoter
] <- "Promoter_TSS"

# ------------------------------------------------------------
# Distance to TSS
#
# Use peak midpoint
# ------------------------------------------------------------

jak3_peaks$Peak_midpoint <-
  floor(
    (
      jak3_peaks$start +
      jak3_peaks$end
    ) / 2
  )

jak3_peaks$Distance_to_TSS_bp <-
  jak3_peaks$Peak_midpoint - tss

jak3_peaks$Abs_distance_to_TSS_bp <-
  abs(
    jak3_peaks$Distance_to_TSS_bp
  )

# ------------------------------------------------------------
# Sort by distance to TSS
# ------------------------------------------------------------

jak3_peaks <- jak3_peaks[
  order(
    jak3_peaks$Abs_distance_to_TSS_bp
  ),
]

# ------------------------------------------------------------
# Save all +/-50 kb peaks
# ------------------------------------------------------------

outfile_all <- file.path(
  outdir,
  "JAK3_ATAC_peaks_gene_plusminus50kb.csv"
)

write.csv(
  jak3_peaks,
  outfile_all,
  row.names = FALSE
)

# ------------------------------------------------------------
# Save promoter peaks separately
# ------------------------------------------------------------

promoter <- jak3_peaks[
  jak3_peaks$Overlap_promoter,
]

outfile_promoter <- file.path(
  outdir,
  "JAK3_ATAC_promoter_peaks.csv"
)

write.csv(
  promoter,
  outfile_promoter,
  row.names = FALSE
)

# ------------------------------------------------------------
# Save gene-body peaks separately
# ------------------------------------------------------------

gene_body <- jak3_peaks[
  jak3_peaks$Overlap_gene_body,
]

outfile_gene <- file.path(
  outdir,
  "JAK3_ATAC_gene_body_peaks.csv"
)

write.csv(
  gene_body,
  outfile_gene,
  row.names = FALSE
)

# ------------------------------------------------------------
# Console report
# ------------------------------------------------------------

cat("JAK3 coordinates:\n")
cat(
  " Gene:",
  chr,
  paste0(gene_start, "-", gene_end),
  "\n"
)

cat(
  " TSS:",
  tss,
  "(minus strand)\n"
)

cat(
  " Promoter:",
  paste0(
    chr,
    ":",
    promoter_start,
    "-",
    promoter_end
  ),
  "\n"
)

cat(
  " Gene +/-50 kb:",
  paste0(
    chr,
    ":",
    nearby_start,
    "-",
    nearby_end
  ),
  "\n\n"
)

cat(
  "Peaks within JAK3 +/-50 kb:",
  nrow(jak3_peaks),
  "\n"
)

cat(
  "Promoter/TSS peaks:",
  nrow(promoter),
  "\n"
)

cat(
  "Gene-body overlapping peaks:",
  nrow(gene_body),
  "\n\n"
)

cat("Region categories:\n")
print(
  table(
    jak3_peaks$JAK3_region
  )
)

cat("\n========================================\n")
cat("JAK3 PEAKS NEAREST TO TSS\n")
cat("========================================\n\n")

show <- jak3_peaks[
  ,
  c(
    "seqnames",
    "start",
    "end",
    "name",
    "score",
    "annotation",
    "JAK3_region",
    "Distance_to_TSS_bp"
  )
]

print(
  head(
    show,
    20
  ),
  row.names = FALSE
)

cat("\nSaved:\n")
cat(" ", outfile_all, "\n")
cat(" ", outfile_promoter, "\n")
cat(" ", outfile_gene, "\n")

cat("\n========================================\n")
cat("58b COMPLETE\n")
cat("========================================\n")
