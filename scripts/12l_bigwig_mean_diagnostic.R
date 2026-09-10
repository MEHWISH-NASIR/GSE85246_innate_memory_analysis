rm(list = ls())

library(rtracklayer)
library(GenomicRanges)

outdir <- "results/12_kinase_reevaluation"

prom <- read.csv(
  file.path(outdir, "12d_hg19_TSS_2kb_promoters.csv"),
  stringsAsFactors = FALSE
)

targets <- c("JAK3", "EPHB2")

prom <- prom[prom$SYMBOL %in% targets, ]

files <- list.files(
  "data/chip/day6_formal/H3K27ac",
  pattern = "\\.bw$",
  full.names = TRUE
)

files <- sort(files)

peer_factor <- 1.16507065275444

res <- list()
k <- 1

for (g in seq_len(nrow(prom))) {

  p <- GRanges(
    seqnames = prom$chr[g],
    ranges = IRanges(
      start = prom$start[g],
      end   = prom$end[g]
    ),
    strand = "*"
  )

  for (f in files) {

    bw <- BigWigFile(f)

    # Current 12d-style summary
    s <- rtracklayer::summary(
      bw,
      p,
      size = 1L,
      type = "mean",
      defaultValue = 0,
      as = "matrix"
    )

    summary_mean <- as.numeric(s)

    # Exact BigWig intervals
    x <- rtracklayer::import(
      bw,
      which = p,
      as = "GRanges"
    )

    if (length(x) == 0) {

      covered <- 0
      signal_sum <- 0

    } else {

      ov_start <- pmax(
        start(x),
        start(p)
      )

      ov_end <- pmin(
        end(x),
        end(p)
      )

      ov_width <- pmax(
        0,
        ov_end - ov_start + 1
      )

      covered <- sum(ov_width)

      signal_sum <- sum(
        ov_width * x$score
      )
    }

    mean0 <- signal_sum / width(p)

    mean_covered <- if (covered > 0) {
      signal_sum / covered
    } else {
      0
    }

    sample <- basename(f)

    condition <- if (
      grepl("^LPS_", sample)
    ) {
      "LPS"
    } else if (
      grepl("^RPMI_", sample)
    ) {
      "RPMI"
    } else {
      "BG"
    }

    replicate <- if (
      grepl("_rep1_", sample)
    ) {
      "rep1"
    } else {
      "rep2"
    }

    scale_factor <- if (
      condition == "RPMI" &&
      replicate == "rep2"
    ) {
      peer_factor
    } else {
      1
    }

    res[[k]] <- data.frame(
      SYMBOL = prom$SYMBOL[g],
      chr = prom$chr[g],
      start = prom$start[g],
      end = prom$end[g],
      sample = sample,
      condition = condition,
      replicate = replicate,
      coverage_fraction =
        covered / width(p),
      rtrack_summary =
        summary_mean,
      exact_mean0 =
        mean0,
      exact_mean_covered =
        mean_covered,
      corrected_summary =
        summary_mean * scale_factor,
      corrected_mean0 =
        mean0 * scale_factor,
      corrected_mean_covered =
        mean_covered * scale_factor,
      stringsAsFactors = FALSE
    )

    k <- k + 1
  }
}

res <- do.call(rbind, res)

write.csv(
  res,
  file.path(
    outdir,
    "12l_H3K27ac_mean_definition_diagnostic.csv"
  ),
  row.names = FALSE
)

cat("\n===== RAW / EXACT SIGNAL COMPARISON =====\n\n")

print(
  res[
    ,
    c(
      "SYMBOL",
      "condition",
      "replicate",
      "coverage_fraction",
      "rtrack_summary",
      "exact_mean0",
      "exact_mean_covered"
    )
  ],
  row.names = FALSE
)

cat("\n===== CORRECTED LPS-vs-RPMI EFFECTS =====\n\n")

metrics <- c(
  "corrected_summary",
  "corrected_mean0",
  "corrected_mean_covered"
)

effects <- list()
z <- 1

for (gene in targets) {

  d <- res[res$SYMBOL == gene, ]

  for (m in metrics) {

    lps <- d[
      d$condition == "LPS",
      m
    ]

    rpmi <- d[
      d$condition == "RPMI",
      m
    ]

    effect <-
      mean(log2(lps + 1)) -
      mean(log2(rpmi + 1))

    effects[[z]] <- data.frame(
      SYMBOL = gene,
      metric = m,
      log2FC_LPS_vs_RPMI = effect
    )

    z <- z + 1
  }
}

effects <- do.call(rbind, effects)

write.csv(
  effects,
  file.path(
    outdir,
    "12l_H3K27ac_effect_comparison.csv"
  ),
  row.names = FALSE
)

print(
  effects,
  row.names = FALSE
)
