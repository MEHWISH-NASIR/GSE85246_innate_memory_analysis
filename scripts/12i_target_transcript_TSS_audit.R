# ============================================================
# 12i — TARGET TRANSCRIPT / TSS PROMOTER AUDIT
#
# Purpose:
#   Determine whether the remaining discrepancy with Tobias
#   is caused by promoter/TSS definition.
#
# Genes:
#   JAK3
#   EPHB2
#   MET
#   MAP3K8
#   BMPR1A
#
# For every UCSC hg19 knownGene transcript:
#   TSS +/- 2 kb
#
# Existing Day-6 BigWigs are reused.
# NO new downloads.
#
# This is a coordinate/direction audit, NOT the final
# formal significance test.
# ============================================================

rm(list = ls())

required <- c(
  "rtracklayer",
  "GenomicRanges",
  "GenomicFeatures",
  "GenomeInfoDb",
  "TxDb.Hsapiens.UCSC.hg19.knownGene"
)

missing <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0) {
  stop(
    paste(
      "Missing packages:",
      paste(missing, collapse = ", ")
    )
  )
}

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)

outdir <- "results/12_kinase_reevaluation"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

rds_file <- file.path(
  outdir,
  "12e_day6_promoter_signal_rescaled.rds"
)

old_promoter_file <- file.path(
  outdir,
  "12d_hg19_TSS_2kb_promoters.csv"
)

if (!file.exists(rds_file)) {
  stop("Missing Step 12e RDS.")
}

if (!file.exists(old_promoter_file)) {
  stop("Missing Step 12d promoter manifest.")
}

x <- readRDS(rds_file)

sample_manifest <- x$sample_manifest
scaling_factors <- x$scaling_factors

old_promoters <- read.csv(
  old_promoter_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("12i — TARGET TRANSCRIPT/TSS AUDIT\n")
cat("========================================\n\n")


# ============================================================
# 1. TARGET GENES
# ============================================================

targets <- data.frame(
  SYMBOL = c(
    "JAK3",
    "EPHB2",
    "MET",
    "MAP3K8",
    "BMPR1A"
  ),
  ENTREZID = c(
    "3718",
    "2048",
    "4233",
    "1326",
    "657"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 2. LOAD hg19 UCSC KNOWNGENE TRANSCRIPTS
# ============================================================

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

tx_by_gene <- GenomicFeatures::transcriptsBy(
  txdb,
  by = "gene"
)

cat(
  "Genes represented in TxDb:",
  length(tx_by_gene),
  "\n"
)


# ============================================================
# 3. EXTRACT ALL TRANSCRIPTS FOR TARGET GENES
# ============================================================

transcript_tables <- list()
promoter_ranges <- list()

for (i in seq_len(nrow(targets))) {

  symbol <- targets$SYMBOL[i]
  entrez <- targets$ENTREZID[i]

  if (!(entrez %in% names(tx_by_gene))) {
    stop(
      paste(
        "Gene missing from TxDb:",
        symbol,
        entrez
      )
    )
  }

  tx <- tx_by_gene[[entrez]]

  tx_id <- as.character(
    mcols(tx)$tx_id
  )

  tx_name <- as.character(
    mcols(tx)$tx_name
  )

  if (is.null(tx_name) ||
      length(tx_name) != length(tx)) {
    tx_name <- rep(
      NA_character_,
      length(tx)
    )
  }

  # Transcript-specific TSS
  tss <- ifelse(
    as.character(strand(tx)) == "+",
    start(tx),
    end(tx)
  )

  prom <- GenomicRanges::promoters(
    tx,
    upstream = 2000,
    downstream = 2000
  )

  prom <- GenomicRanges::trim(
    prom
  )

  df <- data.frame(
    SYMBOL = symbol,
    ENTREZID = entrez,
    tx_id = tx_id,
    tx_name = tx_name,
    chr =
      as.character(seqnames(prom)),
    promoter_start =
      start(prom),
    promoter_end =
      end(prom),
    strand =
      as.character(strand(prom)),
    TSS = tss,
    stringsAsFactors = FALSE
  )

  transcript_tables[[symbol]] <- df

  mcols(prom)$SYMBOL <- symbol
  mcols(prom)$ENTREZID <- entrez
  mcols(prom)$tx_id <- tx_id
  mcols(prom)$tx_name <- tx_name
  mcols(prom)$TSS <- tss

  promoter_ranges[[symbol]] <- prom
}

transcript_df <- do.call(
  rbind,
  transcript_tables
)

promoter_gr <- do.call(
  c,
  unname(promoter_ranges)
)

rownames(transcript_df) <- NULL


# ============================================================
# 4. COLLAPSE TRANSCRIPTS SHARING IDENTICAL TSS/PROMOTER
# ============================================================

key <- paste(
  transcript_df$SYMBOL,
  transcript_df$chr,
  transcript_df$promoter_start,
  transcript_df$promoter_end,
  transcript_df$strand,
  sep = "|"
)

split_idx <- split(
  seq_len(nrow(transcript_df)),
  key
)

collapsed <- lapply(
  split_idx,
  function(idx) {

    z <- transcript_df[idx, ]

    data.frame(
      SYMBOL = z$SYMBOL[1],
      ENTREZID = z$ENTREZID[1],
      chr = z$chr[1],
      promoter_start =
        z$promoter_start[1],
      promoter_end =
        z$promoter_end[1],
      strand =
        z$strand[1],
      TSS =
        z$TSS[1],
      N_transcripts =
        length(idx),
      tx_id =
        paste(
          unique(z$tx_id),
          collapse = ";"
        ),
      tx_name =
        paste(
          unique(
            z$tx_name[
              !is.na(z$tx_name)
            ]
          ),
          collapse = ";"
        ),
      stringsAsFactors = FALSE
    )
  }
)

promoter_manifest <- do.call(
  rbind,
  collapsed
)

rownames(promoter_manifest) <- NULL

promoter_manifest <- promoter_manifest[
  order(
    promoter_manifest$SYMBOL,
    promoter_manifest$chr,
    promoter_manifest$TSS
  ),
]

promoter_manifest$Promoter_ID <- paste0(
  promoter_manifest$SYMBOL,
  "_TSS",
  ave(
    promoter_manifest$TSS,
    promoter_manifest$SYMBOL,
    FUN = seq_along
  )
)


# ============================================================
# 5. COMPARE WITH PREVIOUS 12d GENE-LEVEL PROMOTER
# ============================================================

promoter_manifest$Matches_12d_gene_promoter <- FALSE
promoter_manifest$Old_12d_start <- NA_integer_
promoter_manifest$Old_12d_end <- NA_integer_

for (i in seq_len(nrow(promoter_manifest))) {

  z <- old_promoters[
    old_promoters$SYMBOL ==
      promoter_manifest$SYMBOL[i],
  ]

  if (nrow(z) == 1) {

    promoter_manifest$Old_12d_start[i] <-
      z$start

    promoter_manifest$Old_12d_end[i] <-
      z$end

    promoter_manifest$
      Matches_12d_gene_promoter[i] <-
      (
        promoter_manifest$chr[i] ==
          z$chr &&
        promoter_manifest$
          promoter_start[i] ==
          z$start &&
        promoter_manifest$
          promoter_end[i] ==
          z$end
      )
  }
}

write.csv(
  promoter_manifest,
  file.path(
    outdir,
    "12i_target_transcript_promoters.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. REBUILD GRanges FROM UNIQUE PROMOTERS
# ============================================================

audit_gr <- GRanges(
  seqnames =
    promoter_manifest$chr,
  ranges = IRanges(
    start =
      promoter_manifest$
      promoter_start,
    end =
      promoter_manifest$
      promoter_end
  ),
  strand =
    promoter_manifest$strand
)


# ============================================================
# 7. QUANTIFY EXISTING 18 BIGWIGS
# ============================================================

quantify_bigwig <- function(
    bw_file,
    ranges
) {

  bw <- BigWigFile(
    bw_file
  )

  output <- rep(
    NA_real_,
    length(ranges)
  )

  chromosomes <- unique(
    as.character(
      seqnames(ranges)
    )
  )

  for (chr in chromosomes) {

    idx <- which(
      as.character(
        seqnames(ranges)
      ) == chr
    )

    gr_chr <- ranges[idx]

    ans <- rtracklayer::summary(
      bw,
      gr_chr,
      size = 1L,
      type = "mean",
      defaultValue = 0,
      as = "matrix"
    )

    ans <- as.numeric(ans)

    if (length(ans) != length(idx)) {
      stop(
        paste(
          "BigWig summary mismatch:",
          basename(bw_file),
          chr,
          "expected",
          length(idx),
          "obtained",
          length(ans)
        )
      )
    }

    output[idx] <- ans
  }

  output
}


raw_signal <- matrix(
  NA_real_,
  nrow = nrow(promoter_manifest),
  ncol = nrow(sample_manifest)
)

rownames(raw_signal) <-
  promoter_manifest$Promoter_ID

colnames(raw_signal) <-
  sample_manifest$sample_id

cat("\n========================================\n")
cat("QUANTIFYING TRANSCRIPT PROMOTERS\n")
cat("========================================\n\n")

for (i in seq_len(nrow(sample_manifest))) {

  cat(
    sprintf(
      "[%02d/%02d] %s\n",
      i,
      nrow(sample_manifest),
      sample_manifest$sample_id[i]
    )
  )

  raw_signal[, i] <-
    quantify_bigwig(
      sample_manifest$file[i],
      audit_gr
    )
}


# ============================================================
# 8. APPLY BOTH STEP-12e SCALING SCHEMES
# ============================================================

peer_signal <- raw_signal
rpm1_signal <- raw_signal

marks <- c(
  "H3K27ac",
  "H3K4me1",
  "H3K4me3"
)

for (m in marks) {

  sf <- scaling_factors[
    scaling_factors$mark == m,
  ]

  if (nrow(sf) != 1) {
    stop(
      paste(
        "Scaling factor missing:",
        m
      )
    )
  }

  col <- paste(
    m,
    "RPMI",
    "rep2",
    sep = "_"
  )

  peer_signal[, col] <-
    raw_signal[, col] *
    sf$primary_scaling_factor

  rpm1_signal[, col] <-
    raw_signal[, col] *
    sf$sensitivity_scaling_factor
}

matrices <- list(
  PeerMedian = peer_signal,
  RPMIrep1Matched = rpm1_signal
)


# ============================================================
# 9. EFFECT AUDIT
#
# Use log2(signal + 0.01), because 12h showed that this
# transformation reproduced JAK3 H3K4me3 almost exactly.
#
# This step focuses on DIRECTION and TSS dependence.
# ============================================================

effects <- list()

counter <- 0L

for (scaling_name in names(matrices)) {

  mat <- matrices[[scaling_name]]

  logmat <- log2(
    mat + 0.01
  )

  for (mark_name in marks) {

    idx <- which(
      sample_manifest$mark ==
        mark_name
    )

    condition <-
      sample_manifest$condition[idx]

    replicate <-
      sample_manifest$replicate[idx]

    m <- logmat[
      ,
      idx,
      drop = FALSE
    ]

    get_col <- function(
        cond,
        rep
    ) {

      which(
        condition == cond &
        replicate == rep
      )
    }

    rpm1 <- m[
      ,
      get_col(
        "RPMI",
        "rep1"
      )
    ]

    rpm2 <- m[
      ,
      get_col(
        "RPMI",
        "rep2"
      )
    ]

    lps1 <- m[
      ,
      get_col(
        "LPS",
        "rep1"
      )
    ]

    lps2 <- m[
      ,
      get_col(
        "LPS",
        "rep2"
      )
    ]

    bg1 <- m[
      ,
      get_col(
        "BG",
        "rep1"
      )
    ]

    bg2 <- m[
      ,
      get_col(
        "BG",
        "rep2"
      )
    ]

    rpm_mean <- rowMeans(
      cbind(rpm1, rpm2)
    )

    lps_mean <- rowMeans(
      cbind(lps1, lps2)
    )

    bg_mean <- rowMeans(
      cbind(bg1, bg2)
    )

    z <- data.frame(
      promoter_manifest,
      Scaling =
        scaling_name,
      Mark =
        mark_name,
      RPMI_rep1_log =
        rpm1,
      RPMI_rep2_log =
        rpm2,
      LPS_rep1_log =
        lps1,
      LPS_rep2_log =
        lps2,
      BG_rep1_log =
        bg1,
      BG_rep2_log =
        bg2,
      RPMI_mean_log =
        rpm_mean,
      LPS_mean_log =
        lps_mean,
      BG_mean_log =
        bg_mean,
      LPS_minus_RPMI =
        lps_mean - rpm_mean,
      BG_minus_RPMI =
        bg_mean - rpm_mean,
      BG_minus_LPS =
        bg_mean - lps_mean,
      LPS_vs_RPMI_direction =
        ifelse(
          lps_mean - rpm_mean > 0,
          "UP",
          ifelse(
            lps_mean - rpm_mean < 0,
            "DOWN",
            "NO_CHANGE"
          )
        ),
      stringsAsFactors = FALSE
    )

    counter <- counter + 1L

    effects[[counter]] <- z
  }
}

effects <- do.call(
  rbind,
  effects
)

rownames(effects) <- NULL

write.csv(
  effects,
  file.path(
    outdir,
    "12i_target_transcript_TSS_effects.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. SAVE RAW SIGNAL MATRIX
# ============================================================

raw_table <- cbind(
  promoter_manifest,
  as.data.frame(
    raw_signal,
    check.names = FALSE
  )
)

write.csv(
  raw_table,
  file.path(
    outdir,
    "12i_target_transcript_raw_signals.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. CONSOLE REPORT
# ============================================================

cat("\n========================================\n")
cat("UNIQUE TRANSCRIPT-TSS COUNTS\n")
cat("========================================\n\n")

print(
  table(
    promoter_manifest$SYMBOL
  )
)

cat("\n========================================\n")
cat("12d PROMOTER MATCH CHECK\n")
cat("========================================\n\n")

print(
  promoter_manifest[
    ,
    c(
      "Promoter_ID",
      "SYMBOL",
      "chr",
      "TSS",
      "promoter_start",
      "promoter_end",
      "N_transcripts",
      "Matches_12d_gene_promoter",
      "Old_12d_start",
      "Old_12d_end"
    )
  ],
  row.names = FALSE
)


show_gene_mark <- function(
    gene,
    mark
) {

  z <- effects[
    effects$SYMBOL == gene &
    effects$Mark == mark &
    effects$Scaling == "PeerMedian",
  ]

  z <- z[
    order(
      -z$LPS_minus_RPMI
    ),
  ]

  print(
    z[
      ,
      c(
        "Promoter_ID",
        "tx_name",
        "TSS",
        "N_transcripts",
        "Matches_12d_gene_promoter",
        "RPMI_mean_log",
        "LPS_mean_log",
        "BG_mean_log",
        "LPS_minus_RPMI",
        "BG_minus_RPMI",
        "BG_minus_LPS",
        "LPS_vs_RPMI_direction"
      )
    ],
    row.names = FALSE
  )
}


cat("\n========================================\n")
cat("EPHB2 H3K27ac — ALL TSSs\n")
cat("========================================\n\n")

show_gene_mark(
  "EPHB2",
  "H3K27ac"
)

cat("\n========================================\n")
cat("JAK3 H3K27ac — ALL TSSs\n")
cat("========================================\n\n")

show_gene_mark(
  "JAK3",
  "H3K27ac"
)

cat("\n========================================\n")
cat("JAK3 H3K4me3 — ALL TSSs\n")
cat("========================================\n\n")

show_gene_mark(
  "JAK3",
  "H3K4me3"
)

cat("\n========================================\n")
cat("MET H3K27ac — ALL TSSs\n")
cat("========================================\n\n")

show_gene_mark(
  "MET",
  "H3K27ac"
)


# ============================================================
# 12. SIMPLE DIAGNOSTIC SUMMARY
# ============================================================

diagnostic <- do.call(
  rbind,
  lapply(
    c(
      "JAK3",
      "EPHB2",
      "MET",
      "MAP3K8",
      "BMPR1A"
    ),
    function(g) {

      z <- effects[
        effects$SYMBOL == g &
        effects$Scaling == "PeerMedian" &
        effects$Mark == "H3K27ac",
      ]

      data.frame(
        SYMBOL = g,
        N_unique_TSS =
          nrow(z),
        N_TSS_LPS_UP =
          sum(
            z$LPS_minus_RPMI > 0
          ),
        N_TSS_LPS_DOWN =
          sum(
            z$LPS_minus_RPMI < 0
          ),
        Max_LPS_minus_RPMI =
          max(
            z$LPS_minus_RPMI
          ),
        Min_LPS_minus_RPMI =
          min(
            z$LPS_minus_RPMI
          ),
        Any_alternative_TSS_UP =
          any(
            z$LPS_minus_RPMI > 0
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  diagnostic,
  file.path(
    outdir,
    "12i_TSS_direction_diagnostic.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("H3K27ac TSS DIRECTION SUMMARY\n")
cat("========================================\n\n")

print(
  diagnostic,
  row.names = FALSE
)

cat("\nSaved:\n")

cat(
  " ",
  file.path(
    outdir,
    "12i_target_transcript_promoters.csv"
  ),
  "\n"
)

cat(
  " ",
  file.path(
    outdir,
    "12i_target_transcript_TSS_effects.csv"
  ),
  "\n"
)

cat(
  " ",
  file.path(
    outdir,
    "12i_TSS_direction_diagnostic.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("12i COMPLETE\n")
cat("========================================\n")
