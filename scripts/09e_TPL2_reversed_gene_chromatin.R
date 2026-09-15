# ============================================================
# 09e_TPL2_reversed_gene_chromatin.R
#
# GSE85246 / GSE85245
#
# Question:
# Do the persistent genes significantly reversed by independent
# TPL2/MAP3K8 inhibition also show persistent/reversible
# H3K27ac and/or H3K4me1 changes in GSE85246?
#
# Targets originate from 09d:
#   MAP3K7CL
#   TNIP3
#   TTC39B
#
# Chromatin logic is intentionally matched to corrected Step 07b.
# ============================================================

rm(list = ls())

# ============================================================
# 1. PACKAGES
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
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste(
      "Missing packages:",
      paste(missing_packages, collapse = ", ")
    )
  )
}

# ============================================================
# 2. INPUT / OUTPUT
# ============================================================

reversal_file <-
  "results/09_MAP3K8_validation/09d_significant_TPL2_reversal_genes.csv"

existing_signal_file <-
  "results/07_epigenetic_persistence/07_all_region_BigWig_signal.csv"

outdir <-
  "results/09_MAP3K8_validation"

figdir <-
  "figures/09_MAP3K8_validation"

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

if (!file.exists(reversal_file)) {
  stop("Missing 09d TPL2 reversal table.")
}

if (!file.exists(existing_signal_file)) {
  stop("Missing Step-07 BigWig signal matrix.")
}

# ============================================================
# 3. LOAD THE THREE FUNCTIONAL TARGETS
# ============================================================

functional <- read.csv(
  reversal_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

target_symbols <- unique(functional$SYMBOL)

expected_targets <- c(
  "MAP3K7CL",
  "TNIP3",
  "TTC39B"
)

if (!setequal(target_symbols, expected_targets)) {
  warning(
    paste(
      "09d target set differs from expected:",
      paste(target_symbols, collapse = ", ")
    )
  )
}

cat("\n========================================\n")
cat("09e — TPL2-REVERSED GENE CHROMATIN\n")
cat("========================================\n\n")

cat(
  "Targets:",
  paste(target_symbols, collapse = ", "),
  "\n"
)

# ============================================================
# 4. BIGWIG MANIFEST
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

  day = c(
    rep("Day1", 8),
    rep("Day6", 12)
  ),

  condition = c(
    "RPMI", "RPMI", "LPS", "LPS",
    "RPMI", "RPMI", "LPS", "LPS",

    "RPMI", "RPMI",
    "LPS", "LPS",
    "RESCUE", "RESCUE",

    "RPMI", "RPMI",
    "LPS", "LPS",
    "RESCUE", "RESCUE"
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

if (!all(file.exists(bw$file))) {
  cat("\nMissing BigWigs:\n")
  print(bw$file[!file.exists(bw$file)])
  stop("Required BigWig files missing.")
}

cat(
  "BigWig tracks:",
  nrow(bw),
  "\n"
)

# ============================================================
# 5. SYMBOL -> ENTREZ / ENSEMBL
# ============================================================

entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = target_symbols,
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)

ensembl <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = target_symbols,
  keytype = "SYMBOL",
  column = "ENSEMBL",
  multiVals = "first"
)

targets <- data.frame(
  SYMBOL = target_symbols,
  ENTREZID = unname(entrez[target_symbols]),
  ENSEMBL = unname(ensembl[target_symbols]),
  stringsAsFactors = FALSE
)

print(
  targets,
  row.names = FALSE
)

if (anyNA(targets$ENTREZID)) {
  stop("One or more targets failed SYMBOL -> ENTREZ mapping.")
}

# ============================================================
# 6. hg19 GENE RANGES
# ============================================================

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

gene_gr_all <- GenomicFeatures::genes(
  txdb,
  single.strand.genes.only = TRUE
)

idx <- match(
  targets$ENTREZID,
  names(gene_gr_all)
)

if (anyNA(idx)) {
  print(
    targets[is.na(idx), ],
    row.names = FALSE
  )
  stop("Target absent from hg19 TxDb.")
}

gene_gr <- gene_gr_all[idx]

S4Vectors::mcols(gene_gr)$SYMBOL <-
  targets$SYMBOL

S4Vectors::mcols(gene_gr)$ENTREZID <-
  targets$ENTREZID

S4Vectors::mcols(gene_gr)$ENSEMBL <-
  targets$ENSEMBL

S4Vectors::mcols(gene_gr)$RNA_priority <-
  "TPL2_reversed_memory"

S4Vectors::mcols(gene_gr)$RNA_support_score <-
  NA_real_

# ============================================================
# 7. PROMOTERS ±2 kb
# ============================================================

promoter_gr <- GenomicRanges::promoters(
  gene_gr,
  upstream = 2000,
  downstream = 2000
)

promoter_gr <-
  GenomicRanges::trim(
    promoter_gr
  )

S4Vectors::mcols(promoter_gr)$region_type <-
  "promoter"

S4Vectors::mcols(promoter_gr)$distance_to_TSS <-
  0

S4Vectors::mcols(promoter_gr)$region_id <-
  paste0(
    S4Vectors::mcols(promoter_gr)$SYMBOL,
    "_PROMOTER"
  )

# ============================================================
# 8. TSS
# ============================================================

gene_strand <-
  as.character(
    GenomicRanges::strand(gene_gr)
  )

tss <- ifelse(
  gene_strand == "+",
  GenomicRanges::start(gene_gr),
  GenomicRanges::end(gene_gr)
)

gene_chr <-
  as.character(
    GenomicRanges::seqnames(gene_gr)
  )

chr_lengths <-
  GenomeInfoDb::seqlengths(txdb)

# ============================================================
# 9. DISTAL WINDOWS
# Same as Step 07:
# 1 kb windows, 2 kb -> 100 kb from TSS
# ============================================================

make_bins <- function(
    chr,
    tss,
    strand,
    symbol,
    ensembl,
    entrez,
    chr_length,
    outer = 100000L,
    promoter_exclusion = 2000L,
    bin_width = 1000L
) {

  pieces <- list()

  left_start <-
    max(
      1L,
      as.integer(tss - outer)
    )

  left_end <-
    as.integer(
      tss -
      promoter_exclusion -
      1L
    )

  if (left_end >= left_start) {

    starts <-
      seq(
        left_start,
        left_end,
        by = bin_width
      )

    ends <-
      pmin(
        starts + bin_width - 1L,
        left_end
      )

    pieces[[length(pieces) + 1L]] <-
      GenomicRanges::GRanges(
        seqnames = chr,
        ranges =
          IRanges::IRanges(
            start = starts,
            end = ends
          ),
        strand = "*"
      )
  }

  right_start <-
    as.integer(
      tss +
      promoter_exclusion +
      1L
    )

  right_end <-
    min(
      as.integer(chr_length),
      as.integer(tss + outer)
    )

  if (right_end >= right_start) {

    starts <-
      seq(
        right_start,
        right_end,
        by = bin_width
      )

    ends <-
      pmin(
        starts + bin_width - 1L,
        right_end
      )

    pieces[[length(pieces) + 1L]] <-
      GenomicRanges::GRanges(
        seqnames = chr,
        ranges =
          IRanges::IRanges(
            start = starts,
            end = ends
          ),
        strand = "*"
      )
  }

  if (length(pieces) == 0) {
    return(NULL)
  }

  out <- do.call(
    c,
    pieces
  )

  mid <-
    floor(
      (
        GenomicRanges::start(out) +
        GenomicRanges::end(out)
      ) / 2
    )

  distance_to_tss <-
    if (strand == "+") {
      mid - tss
    } else {
      tss - mid
    }

  S4Vectors::mcols(out)$SYMBOL <- symbol
  S4Vectors::mcols(out)$ENTREZID <- entrez
  S4Vectors::mcols(out)$ENSEMBL <- ensembl

  S4Vectors::mcols(out)$RNA_priority <-
    "TPL2_reversed_memory"

  S4Vectors::mcols(out)$RNA_support_score <-
    NA_real_

  S4Vectors::mcols(out)$region_type <-
    "distal"

  S4Vectors::mcols(out)$distance_to_TSS <-
    distance_to_tss

  S4Vectors::mcols(out)$region_id <-
    paste0(
      symbol,
      "_DISTAL_",
      as.character(
        GenomicRanges::seqnames(out)
      ),
      "_",
      GenomicRanges::start(out),
      "_",
      GenomicRanges::end(out)
    )

  out
}

distal_list <-
  vector(
    "list",
    length(gene_gr)
  )

for (i in seq_along(gene_gr)) {

  this_chr <-
    gene_chr[i]

  distal_list[[i]] <-
    make_bins(
      chr = this_chr,
      tss = tss[i],
      strand = gene_strand[i],
      symbol = targets$SYMBOL[i],
      ensembl = targets$ENSEMBL[i],
      entrez = targets$ENTREZID[i],
      chr_length =
        chr_lengths[this_chr]
    )
}

distal_gr <-
  do.call(
    c,
    distal_list
  )

cat(
  "Promoters:",
  length(promoter_gr),
  "\n"
)

cat(
  "New distal windows:",
  length(distal_gr),
  "\n"
)

# ============================================================
# 10. COMBINE NEW TARGET REGIONS
# ============================================================

all_regions <-
  c(
    promoter_gr,
    distal_gr
  )

GenomicRanges::strand(all_regions) <-
  "*"

region_table <- data.frame(
  region_id =
    S4Vectors::mcols(all_regions)$region_id,

  ENSEMBL =
    S4Vectors::mcols(all_regions)$ENSEMBL,

  SYMBOL =
    S4Vectors::mcols(all_regions)$SYMBOL,

  ENTREZID =
    S4Vectors::mcols(all_regions)$ENTREZID,

  RNA_priority =
    S4Vectors::mcols(all_regions)$RNA_priority,

  RNA_support_score =
    S4Vectors::mcols(all_regions)$RNA_support_score,

  region_type =
    S4Vectors::mcols(all_regions)$region_type,

  chr =
    as.character(
      GenomicRanges::seqnames(all_regions)
    ),

  start =
    GenomicRanges::start(all_regions),

  end =
    GenomicRanges::end(all_regions),

  distance_to_TSS =
    S4Vectors::mcols(all_regions)$distance_to_TSS,

  stringsAsFactors = FALSE
)

# ============================================================
# 11. BIGWIG REGION-MEAN FUNCTION
# Exact Step-07 logic
# ============================================================

bigwig_region_mean <- function(
    file,
    regions
) {

  out <-
    numeric(
      length(regions)
    )

  bf <-
    rtracklayer::BigWigFile(
      file
    )

  bw_seqlevels <-
    GenomeInfoDb::seqlevels(
      GenomeInfoDb::seqinfo(bf)
    )

  keep <-
    as.character(
      GenomicRanges::seqnames(regions)
    ) %in%
    bw_seqlevels

  if (!any(keep)) {
    return(out)
  }

  sub_regions <-
    regions[keep]

  query_ranges <-
    GenomicRanges::reduce(
      sub_regions,
      ignore.strand = TRUE
    )

  GenomicRanges::strand(query_ranges) <-
    "*"

  imported <-
    rtracklayer::import(
      bf,
      which = query_ranges
    )

  if (length(imported) == 0) {
    return(out)
  }

  score <-
    S4Vectors::mcols(imported)$score

  score[!is.finite(score)] <-
    0

  hits <-
    GenomicRanges::findOverlaps(
      sub_regions,
      imported,
      ignore.strand = TRUE
    )

  if (length(hits) == 0) {
    return(out)
  }

  q <-
    S4Vectors::queryHits(hits)

  s <-
    S4Vectors::subjectHits(hits)

  overlap_start <-
    pmax(
      GenomicRanges::start(sub_regions)[q],
      GenomicRanges::start(imported)[s]
    )

  overlap_end <-
    pmin(
      GenomicRanges::end(sub_regions)[q],
      GenomicRanges::end(imported)[s]
    )

  overlap_width <-
    pmax(
      0,
      overlap_end -
      overlap_start +
      1
    )

  weighted_signal <-
    score[s] *
    overlap_width

  summed <-
    rowsum(
      weighted_signal,
      group = q,
      reorder = FALSE
    )

  local_index <-
    as.integer(
      rownames(summed)
    )

  local_mean <-
    numeric(
      length(sub_regions)
    )

  local_mean[local_index] <-
    as.numeric(
      summed[, 1]
    ) /
    GenomicRanges::width(
      sub_regions
    )[local_index]

  out[which(keep)] <-
    local_mean

  out
}

# ============================================================
# 12. QUANTIFY ALL 20 TRACKS
# ============================================================

signal_matrix <- matrix(
  NA_real_,
  nrow = length(all_regions),
  ncol = nrow(bw)
)

rownames(signal_matrix) <-
  region_table$region_id

colnames(signal_matrix) <-
  bw$track_id

cat("\n========================================\n")
cat("QUANTIFYING NEW TARGET REGIONS\n")
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

  signal_matrix[, j] <-
    bigwig_region_mean(
      bw$file[j],
      all_regions
    )

  gc()
}

if (anyNA(signal_matrix)) {
  stop("NA remains in target signal matrix.")
}

target_signal <- cbind(
  region_table,
  as.data.frame(
    signal_matrix,
    check.names = FALSE
  )
)

write.csv(
  target_signal,
  file.path(
    outdir,
    "09e_target_all_region_BigWig_signal.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 13. CONDITION-BLIND FIXED DISTAL SELECTION
#
# Use existing 24-gene distal windows as background
# plus these new target windows.
#
# No Day1/Day6/BG effect is used in selection.
# ============================================================

existing_sig <- read.csv(
  existing_signal_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta_cols <- c(
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

track_cols <-
  bw$track_id

existing_distal <-
  existing_sig[
    existing_sig$region_type ==
      "distal" &
    !existing_sig$SYMBOL %in%
      target_symbols,
    c(
      meta_cols,
      track_cols
    ),
    drop = FALSE
  ]

new_distal <-
  target_signal[
    target_signal$region_type ==
      "distal",
    c(
      meta_cols,
      track_cols
    ),
    drop = FALSE
  ]

selection_pool <-
  rbind(
    existing_distal,
    new_distal
  )

h3k27ac_tracks <-
  bw$track_id[
    bw$mark ==
      "H3K27ac"
  ]

h3k4me1_tracks <-
  bw$track_id[
    bw$mark ==
      "H3K4me1"
  ]

selection_pool$H3K27ac_condition_blind_mean <-
  rowMeans(
    selection_pool[
      ,
      h3k27ac_tracks,
      drop = FALSE
    ]
  )

selection_pool$H3K4me1_condition_blind_mean <-
  rowMeans(
    selection_pool[
      ,
      h3k4me1_tracks,
      drop = FALSE
    ]
  )

percentile_rank <- function(x) {

  ok <- is.finite(x)

  out <-
    rep(
      NA_real_,
      length(x)
    )

  out[ok] <-
    rank(
      x[ok],
      ties.method = "average"
    ) /
    sum(ok)

  out
}

selection_pool$H3K27ac_condition_blind_percentile <-
  percentile_rank(
    selection_pool$H3K27ac_condition_blind_mean
  )

selection_pool$H3K4me1_condition_blind_percentile <-
  percentile_rank(
    selection_pool$H3K4me1_condition_blind_mean
  )

selection_pool$condition_blind_combined_score <-
  rowMeans(
    cbind(
      selection_pool$H3K27ac_condition_blind_percentile,
      selection_pool$H3K4me1_condition_blind_percentile
    )
  )

target_distal_scored <-
  selection_pool[
    selection_pool$SYMBOL %in%
      target_symbols,
    ,
    drop = FALSE
  ]

target_distal_scored$abs_distance_to_TSS <-
  abs(
    target_distal_scored$distance_to_TSS
  )

target_distal_scored <-
  target_distal_scored[
    order(
      target_distal_scored$SYMBOL,
      -target_distal_scored$condition_blind_combined_score,
      target_distal_scored$abs_distance_to_TSS,
      target_distal_scored$chr,
      target_distal_scored$start
    ),
    ,
    drop = FALSE
  ]

selected <-
  target_distal_scored[
    !duplicated(
      target_distal_scored$SYMBOL
    ),
    ,
    drop = FALSE
  ]

if (nrow(selected) != length(target_symbols)) {
  stop(
    "Did not obtain exactly one fixed distal window per target."
  )
}

write.csv(
  selected,
  file.path(
    outdir,
    "09e_condition_blind_selected_distal_windows.csv"
  ),
  row.names = FALSE
)

cat("\nSelected fixed distal windows:\n\n")

print(
  selected[
    ,
    c(
      "SYMBOL",
      "chr",
      "start",
      "end",
      "distance_to_TSS",
      "condition_blind_combined_score"
    )
  ],
  row.names = FALSE
)

# ============================================================
# 14. MARK EFFECT CALCULATION
# Exact corrected 07b criteria
# ============================================================

get_tracks <- function(
    mark,
    day,
    condition
) {

  out <-
    bw$track_id[
      bw$mark == mark &
      bw$day == day &
      bw$condition == condition
    ]

  if (length(out) != 2) {
    stop(
      paste(
        "Expected two tracks:",
        mark,
        day,
        condition
      )
    )
  }

  out
}

calc_mark <- function(
    d,
    mark
) {

  mean_group <- function(
      day,
      condition
  ) {

    tracks <-
      get_tracks(
        mark,
        day,
        condition
      )

    rowMeans(
      d[
        ,
        tracks,
        drop = FALSE
      ]
    )
  }

  d1_rpmi <-
    mean_group(
      "Day1",
      "RPMI"
    )

  d1_lps <-
    mean_group(
      "Day1",
      "LPS"
    )

  d6_rpmi <-
    mean_group(
      "Day6",
      "RPMI"
    )

  d6_lps <-
    mean_group(
      "Day6",
      "LPS"
    )

  d6_rescue <-
    mean_group(
      "Day6",
      "RESCUE"
    )

  d1_delta <-
    d1_lps -
    d1_rpmi

  d6_delta <-
    d6_lps -
    d6_rpmi

  bg_shift <-
    d6_rescue -
    d6_lps

  same_direction <-
    is.finite(d1_delta) &
    is.finite(d6_delta) &
    (
      d1_delta *
      d6_delta >
      0
    )

  bg_opposes_day6 <-
    is.finite(bg_shift) &
    is.finite(d6_delta) &
    (
      bg_shift *
      d6_delta <
      0
    )

  rescue_closer <-
    abs(
      d6_rescue -
      d6_rpmi
    ) <
    abs(
      d6_lps -
      d6_rpmi
    )

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

    cross_time_same_direction =
      same_direction,

    BG_shift_opposes_Day6 =
      bg_opposes_day6,

    BG_closer_to_RPMI =
      rescue_closer,

    persistent_BG_reversible =
      persistent_bg_reversible,

    stringsAsFactors = FALSE
  )
}

# ============================================================
# 15. REGION-LEVEL INTEGRATION
# ============================================================

make_region_result <- function(
    d
) {

  ac <-
    calc_mark(
      d,
      "H3K27ac"
    )

  me1 <-
    calc_mark(
      d,
      "H3K4me1"
    )

  out <-
    d[
      ,
      c(
        "region_id",
        "ENSEMBL",
        "SYMBOL",
        "chr",
        "start",
        "end",
        "distance_to_TSS"
      ),
      drop = FALSE
    ]

  out$H3K27ac_Day1_delta <-
    ac$Day1_delta

  out$H3K27ac_Day6_delta <-
    ac$Day6_delta

  out$H3K27ac_BG_shift <-
    ac$BG_shift

  out$H3K27ac_BG_closer_to_RPMI <-
    ac$BG_closer_to_RPMI

  out$H3K27ac_support <-
    ac$persistent_BG_reversible

  out$H3K4me1_Day1_delta <-
    me1$Day1_delta

  out$H3K4me1_Day6_delta <-
    me1$Day6_delta

  out$H3K4me1_BG_shift <-
    me1$BG_shift

  out$H3K4me1_BG_closer_to_RPMI <-
    me1$BG_closer_to_RPMI

  out$H3K4me1_support <-
    me1$persistent_BG_reversible

  out$mark_support_count <-
    as.integer(
      out$H3K27ac_support
    ) +
    as.integer(
      out$H3K4me1_support
    )

  out$any_mark_support <-
    out$mark_support_count >= 1

  out$dual_mark_support <-
    out$mark_support_count == 2

  out
}

promoter_rows <-
  target_signal[
    target_signal$region_type ==
      "promoter",
    ,
    drop = FALSE
  ]

promoter_result <-
  make_region_result(
    promoter_rows
  )

fixed_result <-
  make_region_result(
    selected
  )

write.csv(
  promoter_result,
  file.path(
    outdir,
    "09e_promoter_chromatin_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  fixed_result,
  file.path(
    outdir,
    "09e_fixed_distal_chromatin_results.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 16. FINAL FUNCTIONAL + CHROMATIN TABLE
# ============================================================

evidence_cols <- intersect(
  c(
    "SYMBOL",
    "Memory_log2FC",
    "Memory_FDR",
    "AcuteLPS_log2FC",
    "AcuteLPS_FDR",
    "TPL2i_log2FC",
    "TPL2i_FDR",
    "MEKi_log2FC",
    "MEKi_FDR"
  ),
  names(functional)
)

final <- functional[
  ,
  evidence_cols,
  drop = FALSE
]

pidx <-
  match(
    final$SYMBOL,
    promoter_result$SYMBOL
  )

didx <-
  match(
    final$SYMBOL,
    fixed_result$SYMBOL
  )

final$Promoter_H3K27ac_support <-
  promoter_result$H3K27ac_support[pidx]

final$Promoter_H3K4me1_support <-
  promoter_result$H3K4me1_support[pidx]

final$Promoter_any_mark_support <-
  promoter_result$any_mark_support[pidx]

final$Promoter_dual_mark_support <-
  promoter_result$dual_mark_support[pidx]

final$Fixed_distal_region <-
  fixed_result$region_id[didx]

final$Fixed_distal_chr <-
  fixed_result$chr[didx]

final$Fixed_distal_start <-
  fixed_result$start[didx]

final$Fixed_distal_end <-
  fixed_result$end[didx]

final$Fixed_distal_H3K27ac_support <-
  fixed_result$H3K27ac_support[didx]

final$Fixed_distal_H3K4me1_support <-
  fixed_result$H3K4me1_support[didx]

final$Fixed_distal_any_mark_support <-
  fixed_result$any_mark_support[didx]

final$Fixed_distal_dual_mark_support <-
  fixed_result$dual_mark_support[didx]

final$Any_corrected_chromatin_support <-
  final$Promoter_any_mark_support |
  final$Fixed_distal_any_mark_support

final$Any_corrected_dual_mark_support <-
  final$Promoter_dual_mark_support |
  final$Fixed_distal_dual_mark_support

final$Integrated_class <- ifelse(
  final$Any_corrected_dual_mark_support,
  "TPL2 reversal + dual-mark chromatin",
  ifelse(
    final$Any_corrected_chromatin_support,
    "TPL2 reversal + single-mark chromatin",
    "TPL2 reversal only"
  )
)

write.csv(
  final,
  file.path(
    outdir,
    "09e_FINAL_TPL2_reversal_chromatin_integration.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 17. CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("FINAL 09e SUMMARY\n")
cat("========================================\n\n")

print(
  final[
    ,
    c(
      "SYMBOL",
      "Memory_log2FC",
      "TPL2i_log2FC",
      "TPL2i_FDR",
      "Promoter_H3K27ac_support",
      "Promoter_H3K4me1_support",
      "Fixed_distal_H3K27ac_support",
      "Fixed_distal_H3K4me1_support",
      "Any_corrected_chromatin_support",
      "Any_corrected_dual_mark_support",
      "Integrated_class"
    )
  ],
  row.names = FALSE
)

cat(
  "\nGenes with any corrected chromatin support:",
  sum(
    final$Any_corrected_chromatin_support
  ),
  "/",
  nrow(final),
  "\n"
)

cat(
  "Genes with corrected dual-mark support:",
  sum(
    final$Any_corrected_dual_mark_support
  ),
  "/",
  nrow(final),
  "\n"
)

cat("\n========================================\n")
cat("09e COMPLETE\n")
cat("========================================\n\n")

cat(
  "Interpretation remains directional/descriptive.\n",
  "This does not prove that MAP3K8 causes chromatin persistence.\n"
)
