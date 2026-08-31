# ============================================================
# 07_epigenetic_persistence.R
#
# GSE85246 / GSE85245
#
# EPIGENETIC PERSISTENCE AT MEMORY-KINASE LOCI
#
# Integrates normalized BigWig signal for H3K27ac and H3K4me1.
#
# Day 1: LPS vs RPMI (2 replicates each)
# Day 6: LPS vs RPMI vs beta-glucan RESCUE (don73, don74)
#
# Candidate universe:
#   24 Day-6 memory kinases from Step 05, retaining RNA priority.
#
# Regulatory regions:
#   promoter = +/-2 kb around hg19 gene TSS
#   distal   = 1-kb non-promoter windows from 2 kb to 100 kb
#              on either side of the TSS
#
# Directional chromatin support requires:
#   1) Day1 and Day6 LPS-RPMI effects have the same direction
#   2) Day6 RESCUE-LPS opposes the Day6 memory direction
#   3) Day6 RESCUE is closer to RPMI than Day6 LPS
#
# IMPORTANT:
#   - descriptive/directional chromatin evidence only
#   - Day1 and Day6 ChIP donor sets are not longitudinally paired
#   - distal windows are operational regions, not called enhancers
#
# Next: 08_make_final_chromatin_figures.R
# ============================================================

rm(list = ls())

# ============================================================
# 1. REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "GenomicRanges",
  "IRanges",
  "rtracklayer",
  "GenomicFeatures",
  "GenomeInfoDb",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste(
      "Missing required package(s):",
      paste(missing_packages, collapse = ", ")
    )
  )
}

# ============================================================
# 2. INPUT / OUTPUT
# ============================================================

rna_file <- paste0(
  "results/05_integrated_trajectory/",
  "05_integrated_kinase_trajectory.csv"
)

if (!file.exists(rna_file)) {
  stop(paste("Missing Step-05 integrated RNA file:", rna_file))
}

outdir <- "results/07_epigenetic_persistence"
figdir <- "figures/07_epigenetic_persistence"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3. BIGWIG MANIFEST
# ============================================================

bw <- data.frame(
  track_id = c(
    "D1_RPMI_r1_H3K27ac",
    "D1_RPMI_r2_H3K27ac",
    "D1_LPS_r1_H3K27ac",
    "D1_LPS_r2_H3K27ac",
    "D1_RPMI_r1_H3K4me1",
    "D1_RPMI_r2_H3K4me1",
    "D1_LPS_r1_H3K4me1",
    "D1_LPS_r2_H3K4me1",
    "D6_RPMI_d73_H3K27ac",
    "D6_RPMI_d74_H3K27ac",
    "D6_LPS_d73_H3K27ac",
    "D6_LPS_d74_H3K27ac",
    "D6_RESCUE_d73_H3K27ac",
    "D6_RESCUE_d74_H3K27ac",
    "D6_RPMI_d73_H3K4me1",
    "D6_RPMI_d74_H3K4me1",
    "D6_LPS_d73_H3K4me1",
    "D6_LPS_d74_H3K4me1",
    "D6_RESCUE_d73_H3K4me1",
    "D6_RESCUE_d74_H3K4me1"
  ),
  day = c(rep("Day1", 8), rep("Day6", 12)),
  condition = c(
    "RPMI", "RPMI", "LPS", "LPS",
    "RPMI", "RPMI", "LPS", "LPS",
    "RPMI", "RPMI", "LPS", "LPS", "RESCUE", "RESCUE",
    "RPMI", "RPMI", "LPS", "LPS", "RESCUE", "RESCUE"
  ),
  replicate = c(
    "rep1", "rep2", "rep1", "rep2",
    "rep1", "rep2", "rep1", "rep2",
    "don73", "don74", "don73", "don74", "don73", "don74",
    "don73", "don74", "don73", "don74", "don73", "don74"
  ),
  mark = c(
    rep("H3K27ac", 4),
    rep("H3K4me1", 4),
    rep("H3K27ac", 6),
    rep("H3K4me1", 6)
  ),
  file = c(
    "data/chip/day1/H3K27ac/RPMI_rep1_H3K27ac.bw",
    "data/chip/day1/H3K27ac/RPMI_rep2_H3K27ac.bw",
    "data/chip/day1/H3K27ac/LPS_rep1_H3K27ac.bw",
    "data/chip/day1/H3K27ac/LPS_rep2_H3K27ac.bw",
    "data/chip/day1/H3K4me1/RPMI_rep1_H3K4me1.bw",
    "data/chip/day1/H3K4me1/RPMI_rep2_H3K4me1.bw",
    "data/chip/day1/H3K4me1/LPS_rep1_H3K4me1.bw",
    "data/chip/day1/H3K4me1/LPS_rep2_H3K4me1.bw",
    "data/chip/H3K27ac/don73_RPMI_H3K27ac.bw",
    "data/chip/H3K27ac/don74_RPMI_H3K27ac.bw",
    "data/chip/H3K27ac/don73_LPS_H3K27ac.bw",
    "data/chip/H3K27ac/don74_LPS_H3K27ac.bw",
    "data/chip/H3K27ac/don73_RESCUE_H3K27ac.bw",
    "data/chip/H3K27ac/don74_RESCUE_H3K27ac.bw",
    "data/chip/H3K4me1/don73_RPMI_H3K4me1.bw",
    "data/chip/H3K4me1/don74_RPMI_H3K4me1.bw",
    "data/chip/H3K4me1/don73_LPS_H3K4me1.bw",
    "data/chip/H3K4me1/don74_LPS_H3K4me1.bw",
    "data/chip/H3K4me1/don73_RESCUE_H3K4me1.bw",
    "data/chip/H3K4me1/don74_RESCUE_H3K4me1.bw"
  ),
  stringsAsFactors = FALSE
)

if (anyDuplicated(bw$track_id)) {
  stop("Duplicate BigWig track IDs.")
}

if (!all(file.exists(bw$file))) {
  cat("\nMissing BigWig files:\n")
  print(bw$file[!file.exists(bw$file)])
  stop("One or more required BigWig files are missing.")
}

write.csv(
  bw,
  file.path(outdir, "07_BigWig_manifest.csv"),
  row.names = FALSE
)

cat("\n========================================\n")
cat("07 — EPIGENETIC PERSISTENCE\n")
cat("========================================\n\n")
cat("BigWig tracks:", nrow(bw), "\n")

# ============================================================
# 4. LOAD RNA MEMORY KINASES
# ============================================================

rna <- read.csv(
  rna_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(rna) != 24) {
  warning(paste("Expected 24 RNA memory kinases but found", nrow(rna)))
}

required_rna_cols <- c(
  "ENSEMBL",
  "SYMBOL",
  "RNA_chromatin_priority",
  "Directional_support_score",
  "Day6_memory_log2FC",
  "Day6_memory_FDR_kinase"
)

missing_rna_cols <- setdiff(required_rna_cols, names(rna))

if (length(missing_rna_cols) > 0) {
  stop(
    paste(
      "Step-05 file is missing required columns:",
      paste(missing_rna_cols, collapse = ", ")
    )
  )
}

cat("RNA memory kinases:", nrow(rna), "\n")
cat(
  "RNA HIGH priority:",
  sum(rna$RNA_chromatin_priority == "HIGH"),
  "\n"
)

# ============================================================
# 5. MAP ENSEMBL -> ENTREZ
# ============================================================

entrez_from_ensembl <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = rna$ENSEMBL,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

rna$ENTREZID <- unname(entrez_from_ensembl[rna$ENSEMBL])

missing_entrez <- is.na(rna$ENTREZID) | rna$ENTREZID == ""

if (any(missing_entrez)) {
  entrez_from_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = rna$SYMBOL[missing_entrez],
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )

  rna$ENTREZID[missing_entrez] <- unname(
    entrez_from_symbol[rna$SYMBOL[missing_entrez]]
  )
}

if (any(is.na(rna$ENTREZID) | rna$ENTREZID == "")) {
  cat("\nKinases without ENTREZ mapping:\n")
  print(
    rna[
      is.na(rna$ENTREZID) | rna$ENTREZID == "",
      c("ENSEMBL", "SYMBOL")
    ],
    row.names = FALSE
  )
  stop("Cannot build hg19 regions for all memory kinases.")
}

write.csv(
  rna,
  file.path(outdir, "07_memory24_with_ENTREZ.csv"),
  row.names = FALSE
)

# ============================================================
# 6. hg19 GENE RANGES
# ============================================================

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

gene_gr_all <- GenomicFeatures::genes(
  txdb,
  single.strand.genes.only = TRUE
)

gene_idx <- match(rna$ENTREZID, names(gene_gr_all))

if (anyNA(gene_idx)) {
  cat("\nENTREZ IDs absent from hg19 TxDb:\n")
  print(
    rna[
      is.na(gene_idx),
      c("ENSEMBL", "SYMBOL", "ENTREZID")
    ],
    row.names = FALSE
  )
  stop("One or more candidate genes are absent from hg19 TxDb.")
}

gene_gr <- gene_gr_all[gene_idx]

S4Vectors::mcols(gene_gr)$ENSEMBL <- rna$ENSEMBL
S4Vectors::mcols(gene_gr)$SYMBOL <- rna$SYMBOL
S4Vectors::mcols(gene_gr)$ENTREZID <- rna$ENTREZID
S4Vectors::mcols(gene_gr)$RNA_priority <- rna$RNA_chromatin_priority
S4Vectors::mcols(gene_gr)$RNA_support_score <- rna$Directional_support_score

# ============================================================
# 7. PROMOTERS: +/-2 kb AROUND TSS
# ============================================================

promoter_gr <- GenomicRanges::promoters(
  gene_gr,
  upstream = 2000,
  downstream = 2000
)

promoter_gr <- GenomicRanges::trim(promoter_gr)

S4Vectors::mcols(promoter_gr)$region_type <- "promoter"
S4Vectors::mcols(promoter_gr)$distance_to_TSS <- 0
S4Vectors::mcols(promoter_gr)$region_id <- paste0(
  S4Vectors::mcols(promoter_gr)$SYMBOL,
  "_PROMOTER"
)

# ============================================================
# 8. TSS INFO
# ============================================================

gene_strand <- as.character(GenomicRanges::strand(gene_gr))

tss <- ifelse(
  gene_strand == "+",
  GenomicRanges::start(gene_gr),
  GenomicRanges::end(gene_gr)
)

gene_chr <- as.character(GenomicRanges::seqnames(gene_gr))
chr_lengths <- GenomeInfoDb::seqlengths(txdb)

# ============================================================
# 9. DISTAL 1-kb WINDOWS: 2 kb -> 100 kb FROM TSS
# ============================================================

make_bins <- function(
    chr,
    tss,
    strand,
    symbol,
    ensembl,
    entrez,
    rna_priority,
    rna_support,
    chr_length,
    outer = 100000L,
    promoter_exclusion = 2000L,
    bin_width = 1000L
) {

  pieces <- list()

  left_start <- max(1L, as.integer(tss - outer))
  left_end <- as.integer(tss - promoter_exclusion - 1L)

  if (left_end >= left_start) {
    starts <- seq(left_start, left_end, by = bin_width)
    ends <- pmin(starts + bin_width - 1L, left_end)

    pieces[[length(pieces) + 1L]] <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(start = starts, end = ends),
      strand = "*"
    )
  }

  right_start <- as.integer(tss + promoter_exclusion + 1L)
  right_end <- min(
    as.integer(chr_length),
    as.integer(tss + outer)
  )

  if (right_end >= right_start) {
    starts <- seq(right_start, right_end, by = bin_width)
    ends <- pmin(starts + bin_width - 1L, right_end)

    pieces[[length(pieces) + 1L]] <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(start = starts, end = ends),
      strand = "*"
    )
  }

  if (length(pieces) == 0) {
    return(NULL)
  }

  out <- do.call(c, pieces)

  mid <- floor(
    (
      GenomicRanges::start(out) +
      GenomicRanges::end(out)
    ) / 2
  )

  distance_to_tss <- if (strand == "+") {
    mid - tss
  } else {
    tss - mid
  }

  S4Vectors::mcols(out)$ENSEMBL <- ensembl
  S4Vectors::mcols(out)$SYMBOL <- symbol
  S4Vectors::mcols(out)$ENTREZID <- entrez
  S4Vectors::mcols(out)$RNA_priority <- rna_priority
  S4Vectors::mcols(out)$RNA_support_score <- rna_support
  S4Vectors::mcols(out)$region_type <- "distal"
  S4Vectors::mcols(out)$distance_to_TSS <- distance_to_tss
  S4Vectors::mcols(out)$region_id <- paste0(
    symbol,
    "_DISTAL_",
    as.character(GenomicRanges::seqnames(out)),
    "_",
    GenomicRanges::start(out),
    "_",
    GenomicRanges::end(out)
  )

  out
}

distal_list <- vector("list", length(gene_gr))

for (i in seq_along(gene_gr)) {

  this_chr <- gene_chr[i]
  this_chr_len <- chr_lengths[this_chr]

  if (is.na(this_chr_len)) {
    stop(paste("Missing hg19 chromosome length for", this_chr))
  }

  distal_list[[i]] <- make_bins(
    chr = this_chr,
    tss = tss[i],
    strand = gene_strand[i],
    symbol = rna$SYMBOL[i],
    ensembl = rna$ENSEMBL[i],
    entrez = rna$ENTREZID[i],
    rna_priority = rna$RNA_chromatin_priority[i],
    rna_support = rna$Directional_support_score[i],
    chr_length = this_chr_len
  )
}

distal_list <- distal_list[
  !vapply(distal_list, is.null, logical(1))
]

distal_gr <- do.call(c, distal_list)

cat("Promoter regions:", length(promoter_gr), "\n")
cat("Distal 1-kb regions:", length(distal_gr), "\n")

# ============================================================
# 10. COMBINE REGIONS
# ============================================================

all_regions <- c(promoter_gr, distal_gr)
GenomicRanges::strand(all_regions) <- "*"

region_table <- data.frame(
  region_id = S4Vectors::mcols(all_regions)$region_id,
  ENSEMBL = S4Vectors::mcols(all_regions)$ENSEMBL,
  SYMBOL = S4Vectors::mcols(all_regions)$SYMBOL,
  ENTREZID = S4Vectors::mcols(all_regions)$ENTREZID,
  RNA_priority = S4Vectors::mcols(all_regions)$RNA_priority,
  RNA_support_score = S4Vectors::mcols(all_regions)$RNA_support_score,
  region_type = S4Vectors::mcols(all_regions)$region_type,
  chr = as.character(GenomicRanges::seqnames(all_regions)),
  start = GenomicRanges::start(all_regions),
  end = GenomicRanges::end(all_regions),
  distance_to_TSS = S4Vectors::mcols(all_regions)$distance_to_TSS,
  stringsAsFactors = FALSE
)

if (anyDuplicated(region_table$region_id)) {
  stop("Duplicate regulatory region IDs were generated.")
}

write.csv(
  region_table,
  file.path(outdir, "07_candidate_regulatory_regions.csv"),
  row.names = FALSE
)

# ============================================================
# 11. BIGWIG REGION-MEAN FUNCTION
# ============================================================

bigwig_region_mean <- function(file, regions) {

  out <- numeric(length(regions))
  bf <- rtracklayer::BigWigFile(file)

  bw_seqlevels <- GenomeInfoDb::seqlevels(
    GenomeInfoDb::seqinfo(bf)
  )

  keep <- as.character(
    GenomicRanges::seqnames(regions)
  ) %in% bw_seqlevels

  if (!any(keep)) {
    return(out)
  }

  sub_regions <- regions[keep]

  query_ranges <- GenomicRanges::reduce(
    sub_regions,
    ignore.strand = TRUE
  )

  GenomicRanges::strand(query_ranges) <- "*"

  imported <- rtracklayer::import(
    bf,
    which = query_ranges
  )

  if (length(imported) == 0) {
    return(out)
  }

  score <- S4Vectors::mcols(imported)$score

  if (is.null(score)) {
    stop(paste("No score column found in BigWig:", file))
  }

  score[!is.finite(score)] <- 0

  hits <- GenomicRanges::findOverlaps(
    sub_regions,
    imported,
    ignore.strand = TRUE
  )

  if (length(hits) == 0) {
    return(out)
  }

  q <- S4Vectors::queryHits(hits)
  s <- S4Vectors::subjectHits(hits)

  overlap_start <- pmax(
    GenomicRanges::start(sub_regions)[q],
    GenomicRanges::start(imported)[s]
  )

  overlap_end <- pmin(
    GenomicRanges::end(sub_regions)[q],
    GenomicRanges::end(imported)[s]
  )

  overlap_width <- pmax(
    0,
    overlap_end - overlap_start + 1
  )

  weighted_signal <- score[s] * overlap_width

  summed <- rowsum(
    weighted_signal,
    group = q,
    reorder = FALSE
  )

  local_index <- as.integer(rownames(summed))

  local_mean <- numeric(length(sub_regions))

  local_mean[local_index] <-
    as.numeric(summed[, 1]) /
    GenomicRanges::width(sub_regions)[local_index]

  out[which(keep)] <- local_mean
  out
}

# ============================================================
# 12. QUANTIFY ALL 20 BIGWIG TRACKS
# ============================================================

signal_matrix <- matrix(
  NA_real_,
  nrow = length(all_regions),
  ncol = nrow(bw)
)

rownames(signal_matrix) <- region_table$region_id
colnames(signal_matrix) <- bw$track_id

cat("\n========================================\n")
cat("QUANTIFYING BIGWIG SIGNAL\n")
cat("========================================\n\n")

for (j in seq_len(nrow(bw))) {

  cat(
    sprintf(
      "[%02d/%02d] %s\n",
      j,
      nrow(bw),
      bw$track_id[j]
    )
  )

  signal_matrix[, j] <- bigwig_region_mean(
    bw$file[j],
    all_regions
  )

  gc()
}

if (anyNA(signal_matrix)) {
  stop("NA values remain in quantified BigWig signal matrix.")
}

signal_df <- cbind(
  region_table,
  as.data.frame(signal_matrix, check.names = FALSE)
)

write.csv(
  signal_df,
  file.path(outdir, "07_all_region_BigWig_signal.csv"),
  row.names = FALSE
)

# ============================================================
# 13. MARK-SPECIFIC EFFECTS
# ============================================================

safe_ratio <- function(numerator, denominator, epsilon = 1e-8) {

  out <- rep(NA_real_, length(numerator))

  ok <-
    is.finite(numerator) &
    is.finite(denominator) &
    abs(denominator) > epsilon

  out[ok] <- numerator[ok] / denominator[ok]
  out
}

calculate_mark_table <- function(
    mark_name,
    region_table,
    signal_matrix,
    bw
) {

  select_tracks <- function(day, condition) {
    bw$track_id[
      bw$mark == mark_name &
      bw$day == day &
      bw$condition == condition
    ]
  }

  d1_rpmi_tracks <- select_tracks("Day1", "RPMI")
  d1_lps_tracks <- select_tracks("Day1", "LPS")
  d6_rpmi_tracks <- select_tracks("Day6", "RPMI")
  d6_lps_tracks <- select_tracks("Day6", "LPS")
  d6_rescue_tracks <- select_tracks("Day6", "RESCUE")

  expected_n <- c(
    length(d1_rpmi_tracks),
    length(d1_lps_tracks),
    length(d6_rpmi_tracks),
    length(d6_lps_tracks),
    length(d6_rescue_tracks)
  )

  if (!all(expected_n == 2)) {
    stop(
      paste(
        "Expected exactly two tracks per condition for",
        mark_name
      )
    )
  }

  d1_rpmi <- rowMeans(
    signal_matrix[, d1_rpmi_tracks, drop = FALSE]
  )
  d1_lps <- rowMeans(
    signal_matrix[, d1_lps_tracks, drop = FALSE]
  )
  d6_rpmi <- rowMeans(
    signal_matrix[, d6_rpmi_tracks, drop = FALSE]
  )
  d6_lps <- rowMeans(
    signal_matrix[, d6_lps_tracks, drop = FALSE]
  )
  d6_rescue <- rowMeans(
    signal_matrix[, d6_rescue_tracks, drop = FALSE]
  )

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

  rescue_closer_to_rpmi <-
    abs(d6_rescue - d6_rpmi) <
    abs(d6_lps - d6_rpmi)

  persistent_bg_reversible <-
    same_direction &
    bg_opposes_day6 &
    rescue_closer_to_rpmi

  rescue_fraction <- safe_ratio(
    d6_lps - d6_rescue,
    d6_lps - d6_rpmi
  )

  effect_strength <- abs(d1_delta) + abs(d6_delta)

  out <- region_table
  out$mark <- mark_name

  out$Day1_RPMI_mean <- d1_rpmi
  out$Day1_LPS_mean <- d1_lps
  out$Day1_LPS_minus_RPMI <- d1_delta

  out$Day6_RPMI_mean <- d6_rpmi
  out$Day6_LPS_mean <- d6_lps
  out$Day6_RESCUE_mean <- d6_rescue
  out$Day6_LPS_minus_RPMI <- d6_delta
  out$Day6_RESCUE_minus_LPS <- bg_shift

  out$cross_time_same_direction <- same_direction
  out$BG_shift_opposes_Day6 <- bg_opposes_day6
  out$BG_closer_to_RPMI <- rescue_closer_to_rpmi
  out$persistent_BG_reversible <- persistent_bg_reversible
  out$BG_rescue_fraction <- rescue_fraction
  out$effect_strength <- effect_strength

  out$strength_percentile <- ave(
    out$effect_strength,
    out$region_type,
    FUN = function(x) {
      rank(
        x,
        ties.method = "average",
        na.last = "keep"
      ) / sum(!is.na(x))
    }
  )

  out
}

ac <- calculate_mark_table(
  "H3K27ac",
  region_table,
  signal_matrix,
  bw
)

me1 <- calculate_mark_table(
  "H3K4me1",
  region_table,
  signal_matrix,
  bw
)

mark_level <- rbind(ac, me1)

write.csv(
  mark_level,
  file.path(outdir, "07_region_mark_level_chromatin.csv"),
  row.names = FALSE
)

write.csv(
  mark_level[
    mark_level$region_type == "promoter",
    ,
    drop = FALSE
  ],
  file.path(outdir, "07_promoter_mark_level_chromatin.csv"),
  row.names = FALSE
)

write.csv(
  mark_level[
    mark_level$region_type == "distal",
    ,
    drop = FALSE
  ],
  file.path(outdir, "07_distal_mark_level_chromatin.csv"),
  row.names = FALSE
)

# ============================================================
# 14. INTEGRATE H3K27ac + H3K4me1 AT SAME REGION
# ============================================================

integrate_marks <- function(ac, me1) {

  idx <- match(ac$region_id, me1$region_id)

  if (anyNA(idx)) {
    stop("H3K27ac/H3K4me1 region sets do not match.")
  }

  out <- ac[
    ,
    c(
      "region_id",
      "ENSEMBL",
      "SYMBOL",
      "ENTREZID",
      "RNA_priority",
      "RNA_support_score",
      "region_type",
      "chr",
      "start",
      "end",
      "distance_to_TSS"
    )
  ]

  out$H3K27ac_Day1_delta <- ac$Day1_LPS_minus_RPMI
  out$H3K27ac_Day6_delta <- ac$Day6_LPS_minus_RPMI
  out$H3K27ac_BG_shift <- ac$Day6_RESCUE_minus_LPS
  out$H3K27ac_BG_rescue_fraction <- ac$BG_rescue_fraction
  out$H3K27ac_persistent_BG_reversible <-
    ac$persistent_BG_reversible
  out$H3K27ac_strength_percentile <- ac$strength_percentile

  out$H3K4me1_Day1_delta <- me1$Day1_LPS_minus_RPMI[idx]
  out$H3K4me1_Day6_delta <- me1$Day6_LPS_minus_RPMI[idx]
  out$H3K4me1_BG_shift <- me1$Day6_RESCUE_minus_LPS[idx]
  out$H3K4me1_BG_rescue_fraction <- me1$BG_rescue_fraction[idx]
  out$H3K4me1_persistent_BG_reversible <-
    me1$persistent_BG_reversible[idx]
  out$H3K4me1_strength_percentile <- me1$strength_percentile[idx]

  out$mark_support_count <-
    as.integer(out$H3K27ac_persistent_BG_reversible) +
    as.integer(out$H3K4me1_persistent_BG_reversible)

  out$dual_mark_support <- out$mark_support_count == 2
  out$any_mark_support <- out$mark_support_count >= 1

  out$combined_strength_percentile <- rowMeans(
    cbind(
      out$H3K27ac_strength_percentile,
      out$H3K4me1_strength_percentile
    ),
    na.rm = TRUE
  )

  out
}

integrated_regions <- integrate_marks(ac, me1)

promoter_integrated <- integrated_regions[
  integrated_regions$region_type == "promoter",
  ,
  drop = FALSE
]

distal_integrated <- integrated_regions[
  integrated_regions$region_type == "distal",
  ,
  drop = FALSE
]

write.csv(
  promoter_integrated,
  file.path(outdir, "07_promoter_integrated_marks.csv"),
  row.names = FALSE
)

write.csv(
  distal_integrated,
  file.path(outdir, "07_distal_integrated_marks.csv"),
  row.names = FALSE
)

# ============================================================
# 15. BEST DISTAL REGION PER GENE
# ============================================================

distal_ordered <- distal_integrated[
  order(
    distal_integrated$SYMBOL,
    -distal_integrated$mark_support_count,
    -distal_integrated$combined_strength_percentile
  ),
  ,
  drop = FALSE
]

best_distal <- distal_ordered[
  !duplicated(distal_ordered$SYMBOL),
  ,
  drop = FALSE
]

write.csv(
  best_distal,
  file.path(outdir, "07_best_distal_region_per_gene.csv"),
  row.names = FALSE
)

# ============================================================
# 16. GENE-LEVEL CHROMATIN SUMMARY
# ============================================================

gene_summary <- rna[
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

pidx <- match(gene_summary$SYMBOL, promoter_integrated$SYMBOL)
didx <- match(gene_summary$SYMBOL, best_distal$SYMBOL)

if (anyNA(pidx)) {
  stop("Missing promoter result for one or more kinases.")
}

if (anyNA(didx)) {
  stop("Missing distal result for one or more kinases.")
}

gene_summary$Promoter_mark_support_count <-
  promoter_integrated$mark_support_count[pidx]

gene_summary$Promoter_dual_mark_support <-
  promoter_integrated$dual_mark_support[pidx]

gene_summary$Promoter_any_mark_support <-
  promoter_integrated$any_mark_support[pidx]

gene_summary$Promoter_combined_strength_percentile <-
  promoter_integrated$combined_strength_percentile[pidx]

gene_summary$Best_distal_region_id <- best_distal$region_id[didx]
gene_summary$Best_distal_chr <- best_distal$chr[didx]
gene_summary$Best_distal_start <- best_distal$start[didx]
gene_summary$Best_distal_end <- best_distal$end[didx]
gene_summary$Best_distal_distance_to_TSS <-
  best_distal$distance_to_TSS[didx]

gene_summary$Distal_mark_support_count <-
  best_distal$mark_support_count[didx]

gene_summary$Distal_dual_mark_support <-
  best_distal$dual_mark_support[didx]

gene_summary$Distal_any_mark_support <-
  best_distal$any_mark_support[didx]

gene_summary$Distal_combined_strength_percentile <-
  best_distal$combined_strength_percentile[didx]

gene_summary$Any_chromatin_support <-
  gene_summary$Promoter_any_mark_support |
  gene_summary$Distal_any_mark_support

gene_summary$Any_dual_mark_support <-
  gene_summary$Promoter_dual_mark_support |
  gene_summary$Distal_dual_mark_support

gene_summary$Chromatin_support_class <- ifelse(
  gene_summary$Any_dual_mark_support,
  "DUAL_MARK_PERSISTENT_BG_REVERSIBLE",
  ifelse(
    gene_summary$Any_chromatin_support,
    "SINGLE_MARK_PERSISTENT_BG_REVERSIBLE",
    "NO_DIRECTIONAL_CHROMATIN_SUPPORT"
  )
)

gene_summary$RNA_plus_chromatin_priority <- ifelse(
  gene_summary$RNA_chromatin_priority == "HIGH" &
    gene_summary$Any_dual_mark_support,
  "TOP_RNA_HIGH_PLUS_DUAL_CHROMATIN",
  ifelse(
    gene_summary$RNA_chromatin_priority == "HIGH" &
      gene_summary$Any_chromatin_support,
    "RNA_HIGH_PLUS_SINGLE_MARK_CHROMATIN",
    ifelse(
      gene_summary$RNA_chromatin_priority == "HIGH",
      "RNA_HIGH_NO_DIRECTIONAL_CHROMATIN",
      ifelse(
        gene_summary$Any_dual_mark_support,
        "NONHIGH_RNA_PLUS_DUAL_CHROMATIN",
        ifelse(
          gene_summary$Any_chromatin_support,
          "NONHIGH_RNA_PLUS_SINGLE_MARK_CHROMATIN",
          "EXPLORATORY"
        )
      )
    )
  )
)

priority_order <- c(
  "TOP_RNA_HIGH_PLUS_DUAL_CHROMATIN" = 1,
  "RNA_HIGH_PLUS_SINGLE_MARK_CHROMATIN" = 2,
  "RNA_HIGH_NO_DIRECTIONAL_CHROMATIN" = 3,
  "NONHIGH_RNA_PLUS_DUAL_CHROMATIN" = 4,
  "NONHIGH_RNA_PLUS_SINGLE_MARK_CHROMATIN" = 5,
  "EXPLORATORY" = 6
)

gene_summary$priority_rank <- unname(
  priority_order[
    gene_summary$RNA_plus_chromatin_priority
  ]
)

gene_summary <- gene_summary[
  order(
    gene_summary$priority_rank,
    gene_summary$Day6_memory_FDR_kinase
  ),
  ,
  drop = FALSE
]

gene_summary$priority_rank <- NULL

write.csv(
  gene_summary,
  file.path(outdir, "07_gene_level_chromatin_summary.csv"),
  row.names = FALSE
)

high_summary <- gene_summary[
  gene_summary$RNA_chromatin_priority == "HIGH",
  ,
  drop = FALSE
]

write.csv(
  high_summary,
  file.path(outdir, "07_HIGH_RNA_priority_chromatin_summary.csv"),
  row.names = FALSE
)

top_integrated <- gene_summary[
  gene_summary$RNA_plus_chromatin_priority ==
    "TOP_RNA_HIGH_PLUS_DUAL_CHROMATIN",
  ,
  drop = FALSE
]

write.csv(
  top_integrated,
  file.path(outdir, "07_TOP_RNA_plus_dual_chromatin.csv"),
  row.names = FALSE
)

# ============================================================
# 17. SIMPLE QC FIGURE
# ============================================================

support_matrix <- cbind(
  Promoter_H3K27ac =
    promoter_integrated$H3K27ac_persistent_BG_reversible[
      match(gene_summary$SYMBOL, promoter_integrated$SYMBOL)
    ],
  Promoter_H3K4me1 =
    promoter_integrated$H3K4me1_persistent_BG_reversible[
      match(gene_summary$SYMBOL, promoter_integrated$SYMBOL)
    ],
  Distal_H3K27ac =
    best_distal$H3K27ac_persistent_BG_reversible[
      match(gene_summary$SYMBOL, best_distal$SYMBOL)
    ],
  Distal_H3K4me1 =
    best_distal$H3K4me1_persistent_BG_reversible[
      match(gene_summary$SYMBOL, best_distal$SYMBOL)
    ]
)

support_matrix[is.na(support_matrix)] <- FALSE
storage.mode(support_matrix) <- "numeric"
rownames(support_matrix) <- gene_summary$SYMBOL

png(
  file.path(
    figdir,
    "07A_chromatin_support_QC_matrix.png"
  ),
  width = 2200,
  height = 2400,
  res = 220
)

par(mar = c(10, 10, 5, 3))

image(
  x = seq_len(ncol(support_matrix)),
  y = seq_len(nrow(support_matrix)),
  z = t(
    support_matrix[
      nrow(support_matrix):1,
      ,
      drop = FALSE
    ]
  ),
  col = c("white", "black"),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "Directional chromatin support at memory-kinase loci"
)

axis(
  1,
  at = seq_len(ncol(support_matrix)),
  labels = colnames(support_matrix),
  las = 2,
  cex.axis = 0.75
)

axis(
  2,
  at = seq_len(nrow(support_matrix)),
  labels = rev(rownames(support_matrix)),
  las = 2,
  cex.axis = 0.75
)

abline(
  v = seq(
    1.5,
    ncol(support_matrix) - 0.5,
    by = 1
  )
)

abline(
  h = seq(
    1.5,
    nrow(support_matrix) - 0.5,
    by = 1
  )
)

dev.off()

# ============================================================
# 18. SAVE RDS
# ============================================================

saveRDS(
  list(
    BigWig_manifest = bw,
    RNA_memory24 = rna,
    promoter_regions = promoter_gr,
    distal_regions = distal_gr,
    region_signal = signal_matrix,
    mark_level = mark_level,
    promoter_integrated = promoter_integrated,
    distal_integrated = distal_integrated,
    best_distal = best_distal,
    gene_summary = gene_summary
  ),
  file.path(
    outdir,
    "07_epigenetic_persistence.rds"
  )
)

# ============================================================
# 19. FINAL SUMMARY
# ============================================================

n_any <- sum(
  gene_summary$Any_chromatin_support,
  na.rm = TRUE
)

n_dual <- sum(
  gene_summary$Any_dual_mark_support,
  na.rm = TRUE
)

n_high_any <- sum(
  high_summary$Any_chromatin_support,
  na.rm = TRUE
)

n_high_dual <- sum(
  high_summary$Any_dual_mark_support,
  na.rm = TRUE
)

n_top <- nrow(top_integrated)

summary_lines <- c(
  "GSE85246 — Step 07 Epigenetic Persistence",
  "",
  paste("RNA memory kinases analysed:", nrow(gene_summary)),
  paste("RNA HIGH-priority kinases:", nrow(high_summary)),
  "",
  paste("Any directional chromatin support:", n_any),
  paste("Any dual-mark chromatin support:", n_dual),
  "",
  paste("RNA HIGH with any chromatin support:", n_high_any),
  paste("RNA HIGH with dual-mark chromatin support:", n_high_dual),
  "",
  paste("Top RNA-HIGH + dual-chromatin candidates:", n_top),
  "",
  "Top integrated symbols:",
  if (n_top > 0) {
    paste(top_integrated$SYMBOL, collapse = ", ")
  } else {
    "None"
  },
  "",
  paste(
    "IMPORTANT:",
    "Chromatin support is directional/descriptive and reflects",
    "same-locus cross-time patterns, not longitudinal donor pairing",
    "or causal inheritance."
  )
)

writeLines(
  summary_lines,
  file.path(
    outdir,
    "07_epigenetic_persistence_summary.txt"
  )
)

cat("\n========================================\n")
cat("07 EPIGENETIC PERSISTENCE COMPLETE\n")
cat("========================================\n\n")

cat(
  "RNA memory kinases analysed:",
  nrow(gene_summary),
  "\n"
)

cat(
  "RNA HIGH-priority kinases:",
  nrow(high_summary),
  "\n"
)

cat(
  "\nAny directional chromatin support:",
  n_any,
  "/",
  nrow(gene_summary),
  "\n"
)

cat(
  "Any dual-mark chromatin support:",
  n_dual,
  "/",
  nrow(gene_summary),
  "\n"
)

cat(
  "\nRNA HIGH with any chromatin support:",
  n_high_any,
  "/",
  nrow(high_summary),
  "\n"
)

cat(
  "RNA HIGH with dual-mark chromatin support:",
  n_high_dual,
  "/",
  nrow(high_summary),
  "\n"
)

cat(
  "\nTOP RNA-HIGH + dual-chromatin candidates:",
  n_top,
  "\n"
)

if (n_top > 0) {
  cat(
    paste(
      top_integrated$SYMBOL,
      collapse = ", "
    ),
    "\n"
  )
}

cat(
  "\nIMPORTANT:\n",
  "Chromatin evidence is directional/descriptive; ",
  "Day1 and Day6 ChIP donor sets are not longitudinally paired.\n",
  sep = ""
)

cat(
  "\nNext canonical step:\n",
  "scripts/08_make_final_chromatin_figures.R\n",
  sep = ""
)

cat("\n========================================\n")
