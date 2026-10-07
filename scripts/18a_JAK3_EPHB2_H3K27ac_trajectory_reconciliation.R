
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(limma)
})

cat("\n========================================\n")
cat("JAK3 / EPHB2 H3K27ac TRAJECTORY RECONCILIATION\n")
cat("DAY1 VS DAY6\n")
cat("========================================\n\n")

outdir <- "results/18_JAK3_EPHB2_reconciliation/H3K27ac"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Promoter definitions
# Same hg19 TSS +/- 2 kb windows used in prior chromatin audit
# ------------------------------------------------------------

prom <- data.frame(
  gene = c("JAK3", "EPHB2"),
  chr = c("chr19", "chr1"),
  start = c(17956881, 23035332),
  end = c(17960880, 23039331),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Sample manifest
# ------------------------------------------------------------

samples <- data.frame(
  sample = c(
    "RPMI_d1_rep1",
    "RPMI_d1_rep2",
    "LPS_d1_rep1",
    "LPS_d1_rep2",
    "RPMI_d6_rep1",
    "RPMI_d6_rep2",
    "LPS_d6_rep1",
    "LPS_d6_rep2"
  ),
  condition = c(
    "RPMI","RPMI","LPS","LPS",
    "RPMI","RPMI","LPS","LPS"
  ),
  day = c(
    "d1","d1","d1","d1",
    "d6","d6","d6","d6"
  ),
  replicate = c(
    "rep1","rep2","rep1","rep2",
    "rep1","rep2","rep1","rep2"
  ),
  file = c(
    "data/chip/day1/H3K27ac/RPMI_rep1_H3K27ac.bw",
    "data/chip/day1/H3K27ac/RPMI_rep2_H3K27ac.bw",
    "data/chip/day1/H3K27ac/LPS_rep1_H3K27ac.bw",
    "data/chip/day1/H3K27ac/LPS_rep2_H3K27ac.bw",
    "data/chip/day6_formal/H3K27ac/RPMI_rep1_H3K27ac.Normalized.bw",
    "data/chip/day6_formal/H3K27ac/RPMI_rep2_H3K27ac.NotNormalized.bw",
    "data/chip/day6_formal/H3K27ac/LPS_rep1_H3K27ac.Normalized.bw",
    "data/chip/day6_formal/H3K27ac/LPS_rep2_H3K27ac.Normalized.bw"
  ),
  stringsAsFactors = FALSE
)

if (!all(file.exists(samples$file))) {
  stop(
    paste(
      "Missing files:",
      paste(samples$file[!file.exists(samples$file)], collapse = "\n")
    )
  )
}

# ------------------------------------------------------------
# 3. Primary and sensitivity scaling factors
# Day6 RPMI rep2 only
# ------------------------------------------------------------

primary_factor <- 1.16507065275444
sensitivity_factor <- 1.52081534357156

# ------------------------------------------------------------
# 4. Exact promoter mean0 extractor
# uncovered bases contribute zero
# ------------------------------------------------------------

extract_mean0 <- function(bw_file, chr, start, end) {

  region <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = start, end = end)
  )

  x <- import(
    BigWigFile(bw_file),
    which = region
  )

  region_width <- width(region)

  if (length(x) == 0) {
    return(0)
  }

  ov_start <- pmax(start(x), start(region))
  ov_end <- pmin(end(x), end(region))
  ov_width <- pmax(0, ov_end - ov_start + 1)

  signal_sum <- sum(
    ov_width * x$score,
    na.rm = TRUE
  )

  signal_sum / region_width
}

# ------------------------------------------------------------
# 5. Extract individual replicate values
# ------------------------------------------------------------

raw <- list()
k <- 1

for (g in seq_len(nrow(prom))) {

  for (i in seq_len(nrow(samples))) {

    v <- extract_mean0(
      samples$file[i],
      prom$chr[g],
      prom$start[g],
      prom$end[g]
    )

    raw[[k]] <- data.frame(
      gene = prom$gene[g],
      chr = prom$chr[g],
      start = prom$start[g],
      end = prom$end[g],
      sample = samples$sample[i],
      condition = samples$condition[i],
      day = samples$day[i],
      replicate = samples$replicate[i],
      raw_mean0 = v,
      stringsAsFactors = FALSE
    )

    k <- k + 1
  }
}

raw <- do.call(rbind, raw)

# ------------------------------------------------------------
# 6. Generate primary + sensitivity corrected datasets
# ------------------------------------------------------------

make_scaled <- function(dat, factor_value, method_name) {

  x <- dat

  x$scale_factor <- 1

  hit <- (
    x$day == "d6" &
    x$condition == "RPMI" &
    x$replicate == "rep2"
  )

  x$scale_factor[hit] <- factor_value

  x$corrected_mean0 <-
    x$raw_mean0 * x$scale_factor

  x$log2_signal <-
    log2(x$corrected_mean0 + 1)

  x$scaling_method <- method_name

  x
}

primary <- make_scaled(
  raw,
  primary_factor,
  "PeerMedian_primary"
)

sensitivity <- make_scaled(
  raw,
  sensitivity_factor,
  "RPMI_rep1_matching_sensitivity"
)

all_values <- rbind(
  primary,
  sensitivity
)

write.csv(
  all_values,
  file.path(
    outdir,
    "18a_JAK3_EPHB2_H3K27ac_individual_replicates.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Descriptive effect-size summaries
# ------------------------------------------------------------

summarise_effects <- function(dat) {

  genes <- unique(dat$gene)

  res <- list()
  k <- 1

  for (g in genes) {

    z <- dat[dat$gene == g, ]

    for (d in c("d1", "d6")) {

      zz <- z[z$day == d, ]

      lps <- zz$log2_signal[
        zz$condition == "LPS"
      ]

      rpmi <- zz$log2_signal[
        zz$condition == "RPMI"
      ]

      effect <-
        mean(lps) - mean(rpmi)

      # Welch-style uncertainty around difference in means
      se <- sqrt(
        var(lps) / length(lps) +
        var(rpmi) / length(rpmi)
      )

      df_num <-
        (
          var(lps) / length(lps) +
          var(rpmi) / length(rpmi)
        )^2

      df_den <-
        (
          (var(lps) / length(lps))^2 /
          (length(lps) - 1)
        ) +
        (
          (var(rpmi) / length(rpmi))^2 /
          (length(rpmi) - 1)
        )

      df <- df_num / df_den

      tcrit <- qt(0.975, df = df)

      res[[k]] <- data.frame(
        gene = g,
        day = d,
        mean_LPS = mean(lps),
        mean_RPMI = mean(rpmi),
        effect_LPS_minus_RPMI = effect,
        SE = se,
        df = df,
        CI_low = effect - tcrit * se,
        CI_high = effect + tcrit * se,
        stringsAsFactors = FALSE
      )

      k <- k + 1
    }
  }

  do.call(rbind, res)
}

primary_effects <- summarise_effects(primary)
primary_effects$scaling_method <- "PeerMedian_primary"

sensitivity_effects <- summarise_effects(sensitivity)
sensitivity_effects$scaling_method <-
  "RPMI_rep1_matching_sensitivity"

effect_table <- rbind(
  primary_effects,
  sensitivity_effects
)

write.csv(
  effect_table,
  file.path(
    outdir,
    "18a_JAK3_EPHB2_H3K27ac_effect_sizes.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Limma model within each gene
#
# Four groups:
# RPMI_d1, LPS_d1, RPMI_d6, LPS_d6
#
# contrasts:
# Day1_effect      = LPS_d1 - RPMI_d1
# Day6_effect      = LPS_d6 - RPMI_d6
# Trajectory       = Day6_effect - Day1_effect
# ------------------------------------------------------------

run_limma <- function(dat, method_name) {

  genes <- c("JAK3", "EPHB2")

  mat <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = 8,
    dimnames = list(
      genes,
      samples$sample
    )
  )

  for (g in genes) {

    z <- dat[dat$gene == g, ]

    mat[g, z$sample] <- z$log2_signal
  }

  group <- factor(
    c(
      "RPMI_d1",
      "RPMI_d1",
      "LPS_d1",
      "LPS_d1",
      "RPMI_d6",
      "RPMI_d6",
      "LPS_d6",
      "LPS_d6"
    ),
    levels = c(
      "RPMI_d1",
      "LPS_d1",
      "RPMI_d6",
      "LPS_d6"
    )
  )

  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)

  fit <- lmFit(
    mat,
    design
  )

  cont <- makeContrasts(
    Day1_effect =
      LPS_d1 - RPMI_d1,

    Day6_effect =
      LPS_d6 - RPMI_d6,

    Trajectory_change =
      (LPS_d6 - RPMI_d6) -
      (LPS_d1 - RPMI_d1),

    levels = design
  )

  fit2 <- contrasts.fit(
    fit,
    cont
  )

  fit2 <- eBayes(
    fit2,
    trend = FALSE
  )

  out <- list()

  for (cc in colnames(cont)) {

    tt <- topTable(
      fit2,
      coef = cc,
      number = Inf,
      adjust.method = "BH",
      sort.by = "none",
      confint = TRUE
    )

    tt$gene <- rownames(tt)
    tt$contrast <- cc
    tt$scaling_method <- method_name

    out[[cc]] <- tt
  }

  do.call(rbind, out)
}

limma_primary <- run_limma(
  primary,
  "PeerMedian_primary"
)

limma_sensitivity <- run_limma(
  sensitivity,
  "RPMI_rep1_matching_sensitivity"
)

limma_results <- rbind(
  limma_primary,
  limma_sensitivity
)

write.csv(
  limma_results,
  file.path(
    outdir,
    "18a_JAK3_EPHB2_H3K27ac_limma_trajectory_results.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Direct JAK3-vs-EPHB2 trajectory comparison
#
# IMPORTANT:
# genes are measured in the same eight sample tracks.
# Here we calculate each gene's interaction estimate first
# and compare the estimates descriptively.
#
# A fully independent gene-by-time interaction test is not
# treated as primary inference because there are only two
# biological replicates per group and genes are not independent
# replicate units.
# ------------------------------------------------------------

primary_traj <- limma_primary[
  limma_primary$contrast == "Trajectory_change",
]

jak <- primary_traj$logFC[
  primary_traj$gene == "JAK3"
]

eph <- primary_traj$logFC[
  primary_traj$gene == "EPHB2"
]

trajectory_compare <- data.frame(
  JAK3_trajectory_change = jak,
  EPHB2_trajectory_change = eph,
  EPHB2_minus_JAK3 = eph - jak,
  Interpretation =
    "Descriptive between-gene difference; primary inference is gene-specific trajectory contrast",
  stringsAsFactors = FALSE
)

write.csv(
  trajectory_compare,
  file.path(
    outdir,
    "18a_JAK3_vs_EPHB2_trajectory_comparison.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Explicit operational trajectory definitions
#
# These definitions intentionally do NOT rely on one time point
# being significant and the other being non-significant.
# ------------------------------------------------------------

definitions <- data.frame(
  Label = c(
    "Persistent",
    "Washout-emergent",
    "Transient",
    "No_supported_trajectory"
  ),
  Operational_definition = c(
    paste(
      "Positive LPS-vs-RPMI effect at Day1 and Day6;",
      "Day6 effect remains directionally positive;",
      "trajectory-change estimate is used to quantify strengthening or weakening."
    ),
    paste(
      "Day6 LPS-vs-RPMI effect is positive and materially larger than Day1;",
      "classification requires effect-size evidence from the direct trajectory contrast,",
      "not merely different significance labels across time points."
    ),
    paste(
      "Positive Day1 LPS-vs-RPMI effect followed by substantial attenuation,",
      "loss, or reversal by Day6."
    ),
    paste(
      "Effect sizes and trajectory contrast do not support one of the above patterns."
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  definitions,
  file.path(
    outdir,
    "18a_trajectory_definitions.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Console report
# ------------------------------------------------------------

cat("\n===== INDIVIDUAL REPLICATE VALUES: PRIMARY SCALING =====\n")
print(
  primary[
    ,
    c(
      "gene",
      "day",
      "condition",
      "replicate",
      "raw_mean0",
      "scale_factor",
      "corrected_mean0",
      "log2_signal"
    )
  ],
  row.names = FALSE
)

cat("\n===== EFFECT SIZE SUMMARY =====\n")
print(
  primary_effects,
  row.names = FALSE
)

cat("\n===== LIMMA TRAJECTORY RESULTS: PRIMARY =====\n")
print(
  limma_primary[
    ,
    c(
      "gene",
      "contrast",
      "logFC",
      "CI.L",
      "CI.R",
      "P.Value",
      "adj.P.Val"
    )
  ],
  row.names = FALSE
)

cat("\n===== JAK3 VS EPHB2 TRAJECTORY =====\n")
print(
  trajectory_compare,
  row.names = FALSE
)

cat("\n18a COMPLETE\n")

