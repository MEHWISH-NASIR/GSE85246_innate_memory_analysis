
rm(list = ls())

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(limma)
})

cat("\n=============================================\n")
cat("18b — JAK3 / EPHB2 H3K27ac RECONCILIATION\n")
cat("ALL-PROMOTER eBAYES DAY1 + DAY6 MODEL\n")
cat("=============================================\n\n")

outdir <- "results/18_JAK3_EPHB2_reconciliation/H3K27ac"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

day6_file <- paste0(
  "results/12_kinase_reevaluation/",
  "12d_day6_promoter_signal_raw.csv"
)

scale_file <- paste0(
  "results/12_kinase_reevaluation/",
  "12e_RPMI_rep2_scaling_factors.csv"
)

kinase_file <- paste0(
  "results/03_memory_kinases/",
  "03_Day6_KinHub_all_kinases.csv"
)

day1_files <- c(
  RPMI_d1_rep1 =
    "data/chip/day1/H3K27ac/RPMI_rep1_H3K27ac.bw",
  RPMI_d1_rep2 =
    "data/chip/day1/H3K27ac/RPMI_rep2_H3K27ac.bw",
  LPS_d1_rep1 =
    "data/chip/day1/H3K27ac/LPS_rep1_H3K27ac.bw",
  LPS_d1_rep2 =
    "data/chip/day1/H3K27ac/LPS_rep2_H3K27ac.bw"
)

required_files <- c(
  day6_file,
  scale_file,
  kinase_file,
  unname(day1_files)
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required files:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}

# ============================================================
# 1. LOAD EXACT DAY6 PROMOTER UNIVERSE
# ============================================================

pm <- read.csv(
  day6_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Full promoter rows loaded:", nrow(pm), "\n")

meta_cols <- c(
  "ENTREZID",
  "SYMBOL",
  "ENSEMBL",
  "chr",
  "start",
  "end",
  "strand",
  "width"
)

if (!all(meta_cols %in% names(pm))) {
  stop("Required promoter metadata columns are missing.")
}

gr <- GRanges(
  seqnames = pm$chr,
  ranges = IRanges(
    start = pm$start,
    end = pm$end
  ),
  strand = "*"
)

# ============================================================
# 2. EXACT BigWig mean0 EXTRACTION
#
# Import one BigWig at a time.
# Signal is overlap-width weighted.
# Uncovered promoter bases contribute zero.
# ============================================================

extract_mean0_all <- function(file, promoters) {

  cat("\nImporting:", basename(file), "\n")

  bw <- import(
    BigWigFile(file)
  )

  cat("BigWig intervals:", length(bw), "\n")

  hits <- findOverlaps(
    promoters,
    bw,
    ignore.strand = TRUE
  )

  q <- queryHits(hits)
  s <- subjectHits(hits)

  signal_sum <- numeric(
    length(promoters)
  )

  if (length(q) > 0) {

    ov_start <- pmax(
      start(promoters)[q],
      start(bw)[s]
    )

    ov_end <- pmin(
      end(promoters)[q],
      end(bw)[s]
    )

    ov_width <- pmax(
      0,
      ov_end - ov_start + 1
    )

    weighted_signal <-
      ov_width * bw$score[s]

    tmp <- rowsum(
      weighted_signal,
      group = q,
      reorder = FALSE
    )

    idx <- as.integer(
      rownames(tmp)
    )

    signal_sum[idx] <- tmp[, 1]
  }

  result <-
    signal_sum /
    width(promoters)

  rm(bw, hits)
  invisible(gc())

  result
}

# ============================================================
# 3. EXTRACT DAY1 FULL-PROMOTER SIGNALS
# ============================================================

day1_cache <- file.path(
  outdir,
  "18b_day1_H3K27ac_full_promoter_signal.csv"
)

if (file.exists(day1_cache)) {

  cat("\nUsing cached Day1 promoter matrix:\n")
  cat(day1_cache, "\n")

  d1 <- read.csv(
    day1_cache,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(d1) != nrow(pm)) {
    stop("Cached Day1 matrix row count differs from Day6 manifest.")
  }

} else {

  cat("\nExtracting Day1 promoter signals...\n")

  d1 <- pm[, meta_cols, drop = FALSE]

  for (nm in names(day1_files)) {

    d1[[nm]] <- extract_mean0_all(
      day1_files[[nm]],
      gr
    )

    cat(
      nm,
      "complete | mean =",
      mean(d1[[nm]]),
      "\n"
    )
  }

  write.csv(
    d1,
    day1_cache,
    row.names = FALSE
  )

  cat("\nDay1 promoter matrix saved.\n")
}

# ============================================================
# 4. DAY6 H3K27ac SIGNALS
# ============================================================

day6_cols <- c(
  "H3K27ac_RPMI_rep1",
  "H3K27ac_RPMI_rep2",
  "H3K27ac_LPS_rep1",
  "H3K27ac_LPS_rep2"
)

if (!all(day6_cols %in% names(pm))) {
  stop(
    paste(
      "Missing Day6 columns:",
      paste(
        setdiff(day6_cols, names(pm)),
        collapse = ", "
      )
    )
  )
}

sc <- read.csv(
  scale_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

h3 <- sc[
  sc$mark == "H3K27ac",
  ,
  drop = FALSE
]

if (nrow(h3) != 1) {
  stop("Could not uniquely resolve H3K27ac scaling row.")
}

primary_factor <-
  h3$primary_scaling_factor

sensitivity_factor <-
  h3$sensitivity_scaling_factor

cat("\nH3K27ac primary scaling factor:",
    primary_factor, "\n")

cat("H3K27ac sensitivity scaling factor:",
    sensitivity_factor, "\n")

# ============================================================
# 5. KINASE UNIVERSE
# ============================================================

kin <- read.csv(
  kinase_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kinase_symbols <- unique(
  kin$SYMBOL[
    !is.na(kin$SYMBOL) &
      kin$SYMBOL != "" &
      kin$Is_KinHub_Kinase
  ]
)

kinase_idx <-
  !is.na(pm$SYMBOL) &
  pm$SYMBOL %in% kinase_symbols

cat(
  "\nKinHub symbols:",
  length(kinase_symbols),
  "\n"
)

cat(
  "Kinase promoter rows:",
  sum(kinase_idx),
  "\n"
)

# ============================================================
# 6. BUILD + FIT FULL 8-SAMPLE MODEL
# ============================================================

fit_trajectory <- function(
  scaling_factor,
  scaling_method
) {

  cat("\n========================================\n")
  cat("SCALING:", scaling_method, "\n")
  cat("========================================\n")

  # Exact sample order:
  #
  # RPMI_d1_rep1
  # RPMI_d1_rep2
  # LPS_d1_rep1
  # LPS_d1_rep2
  # RPMI_d6_rep1
  # RPMI_d6_rep2 corrected
  # LPS_d6_rep1
  # LPS_d6_rep2

  mat_raw <- cbind(
    RPMI_d1_rep1 =
      d1$RPMI_d1_rep1,

    RPMI_d1_rep2 =
      d1$RPMI_d1_rep2,

    LPS_d1_rep1 =
      d1$LPS_d1_rep1,

    LPS_d1_rep2 =
      d1$LPS_d1_rep2,

    RPMI_d6_rep1 =
      pm$H3K27ac_RPMI_rep1,

    RPMI_d6_rep2 =
      pm$H3K27ac_RPMI_rep2 *
      scaling_factor,

    LPS_d6_rep1 =
      pm$H3K27ac_LPS_rep1,

    LPS_d6_rep2 =
      pm$H3K27ac_LPS_rep2
  )

  storage.mode(mat_raw) <- "numeric"

  if (anyNA(mat_raw)) {
    stop(
      paste(
        "NA signal encountered:",
        scaling_method
      )
    )
  }

  if (any(mat_raw < 0)) {
    stop(
      paste(
        "Negative signal encountered:",
        scaling_method
      )
    )
  }

  mat <- log2(
    mat_raw + 1
  )

  rownames(mat) <-
    seq_len(nrow(pm))

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

  design <- model.matrix(
    ~0 + group
  )

  colnames(design) <-
    levels(group)

  contrast_matrix <-
    makeContrasts(

      Day1_effect =
        LPS_d1 - RPMI_d1,

      Day6_effect =
        LPS_d6 - RPMI_d6,

      Trajectory_change =
        (LPS_d6 - RPMI_d6) -
        (LPS_d1 - RPMI_d1),

      levels = design
    )

  cat("\nDesign:\n")
  print(design)

  cat("\nContrasts:\n")
  print(contrast_matrix)

  # ----------------------------------------------------------
  # Critical point:
  # lmFit + eBayes on ALL promoter rows.
  # ----------------------------------------------------------

  fit <- lmFit(
    mat,
    design
  )

  fit2 <- contrasts.fit(
    fit,
    contrast_matrix
  )

  fit2 <- eBayes(
    fit2
  )

  all_results <- list()

  target_results <- list()

  kinase_results <- list()

  for (cc in colnames(contrast_matrix)) {

    tt <- topTable(
      fit2,
      coef = cc,
      number = Inf,
      sort.by = "none",
      adjust.method = "none",
      confint = TRUE
    )

    res <- data.frame(
      pm[, meta_cols, drop = FALSE],
      scaling_method =
        scaling_method,
      contrast =
        cc,
      logFC =
        tt$logFC,
      CI_low =
        tt$CI.L,
      CI_high =
        tt$CI.R,
      AveExpr =
        tt$AveExpr,
      t =
        tt$t,
      P_value =
        tt$P.Value,
      B =
        tt$B,
      stringsAsFactors = FALSE
    )

    # Kinase-restricted FDR
    kres <- res[
      kinase_idx,
      ,
      drop = FALSE
    ]

    kres$FDR_KinHub <-
      p.adjust(
        kres$P_value,
        method = "BH"
      )

    # Put kinase FDR back into full table
    res$FDR_KinHub <- NA_real_

    res$FDR_KinHub[
      which(kinase_idx)
    ] <- kres$FDR_KinHub

    # Prespecified two-gene target correction
    target_idx <-
      !is.na(res$SYMBOL) &
      res$SYMBOL %in%
      c("JAK3", "EPHB2")

    tres <- res[
      target_idx,
      ,
      drop = FALSE
    ]

    tres$FDR_target2 <-
      p.adjust(
        tres$P_value,
        method = "BH"
      )

    all_results[[cc]] <- res
    kinase_results[[cc]] <- kres
    target_results[[cc]] <- tres
  }

  # ----------------------------------------------------------
  # Individual replicate table for JAK3 + EPHB2
  # ----------------------------------------------------------

  target_rows <-
    which(
      pm$SYMBOL %in%
      c("JAK3", "EPHB2")
    )

  repl <- list()
  z <- 1

  sample_names <-
    colnames(mat_raw)

  for (rr in target_rows) {

    for (j in seq_along(sample_names)) {

      repl[[z]] <- data.frame(
        ENTREZID =
          pm$ENTREZID[rr],
        SYMBOL =
          pm$SYMBOL[rr],
        chr =
          pm$chr[rr],
        start =
          pm$start[rr],
        end =
          pm$end[rr],
        sample =
          sample_names[j],
        raw_or_corrected_mean0 =
          mat_raw[rr, j],
        log2_signal =
          mat[rr, j],
        scaling_method =
          scaling_method,
        stringsAsFactors = FALSE
      )

      z <- z + 1
    }
  }

  list(
    all =
      do.call(rbind, all_results),

    kinases =
      do.call(rbind, kinase_results),

    targets =
      do.call(rbind, target_results),

    replicates =
      do.call(rbind, repl)
  )
}

# ============================================================
# 7. PRIMARY MODEL
# ============================================================

primary <- fit_trajectory(
  scaling_factor =
    primary_factor,
  scaling_method =
    "PeerMedian_primary"
)

# ============================================================
# 8. SENSITIVITY MODEL
# ============================================================

sensitivity <- fit_trajectory(
  scaling_factor =
    sensitivity_factor,
  scaling_method =
    "RPMI_rep1_matching_sensitivity"
)

# ============================================================
# 9. SAVE RESULTS
# ============================================================

targets_combined <- rbind(
  primary$targets,
  sensitivity$targets
)

replicates_combined <- rbind(
  primary$replicates,
  sensitivity$replicates
)

kinases_combined <- rbind(
  primary$kinases,
  sensitivity$kinases
)

write.csv(
  targets_combined,
  file.path(
    outdir,
    "18b_JAK3_EPHB2_H3K27ac_allpromoter_trajectory.csv"
  ),
  row.names = FALSE
)

write.csv(
  replicates_combined,
  file.path(
    outdir,
    "18b_JAK3_EPHB2_H3K27ac_replicate_values.csv"
  ),
  row.names = FALSE
)

write.csv(
  kinases_combined,
  file.path(
    outdir,
    "18b_H3K27ac_all_KinHub_trajectory_results.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 10. PRIMARY TARGET SUMMARY
# ============================================================

target_primary <- primary$targets

target_primary$Direction <-
  ifelse(
    target_primary$logFC > 0,
    "POSITIVE",
    ifelse(
      target_primary$logFC < 0,
      "NEGATIVE",
      "ZERO"
    )
  )

write.csv(
  target_primary,
  file.path(
    outdir,
    "18b_PRIMARY_JAK3_EPHB2_H3K27ac_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 11. SENSITIVITY CONCORDANCE
# ============================================================

p <- primary$targets[
  ,
  c(
    "SYMBOL",
    "contrast",
    "logFC",
    "CI_low",
    "CI_high",
    "P_value"
  )
]

s <- sensitivity$targets[
  ,
  c(
    "SYMBOL",
    "contrast",
    "logFC",
    "CI_low",
    "CI_high",
    "P_value"
  )
]

names(p)[3:6] <- paste0(
  names(p)[3:6],
  "_primary"
)

names(s)[3:6] <- paste0(
  names(s)[3:6],
  "_sensitivity"
)

concordance <- merge(
  p,
  s,
  by = c(
    "SYMBOL",
    "contrast"
  ),
  all = TRUE,
  sort = FALSE
)

concordance$Same_direction <-
  sign(
    concordance$logFC_primary
  ) ==
  sign(
    concordance$logFC_sensitivity
  )

write.csv(
  concordance,
  file.path(
    outdir,
    "18b_primary_vs_sensitivity_concordance.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 12. OPERATIONAL DEFINITIONS
#
# IMPORTANT:
# Definitions are based on estimated effects + uncertainty,
# not on "significant at one time but not another".
# ============================================================

definitions <- data.frame(

  label = c(
    "persistent_positive",
    "washout_emergent",
    "transient_positive",
    "trajectory_unresolved"
  ),

  definition = c(

    paste(
      "Positive LPS-vs-RPMI effect at Day1 and Day6,",
      "with the Day6 effect supported by a 95% CI above zero.",
      "The trajectory contrast quantifies strengthening or weakening",
      "but persistence is not inferred from equality of p-values."
    ),

    paste(
      "Day6 LPS-vs-RPMI effect is positive and exceeds the Day1 effect,",
      "with a positive direct trajectory-change estimate;",
      "strong washout-emergent evidence requires the 95% CI for",
      "the trajectory contrast to exclude zero."
    ),

    paste(
      "Positive Day1 LPS-vs-RPMI effect followed by substantial",
      "attenuation, loss or reversal at Day6, evaluated using",
      "the direct trajectory contrast and its uncertainty."
    ),

    paste(
      "The estimated effects or trajectory CI are too uncertain",
      "to support a persistent, emergent or transient trajectory."
    )
  ),

  stringsAsFactors = FALSE
)

write.csv(
  definitions,
  file.path(
    outdir,
    "18b_operational_trajectory_definitions.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 13. CONSOLE OUTPUT
# ============================================================

cat("\n\n========================================\n")
cat("PRIMARY JAK3 / EPHB2 RESULTS\n")
cat("========================================\n")

show_cols <- c(
  "SYMBOL",
  "contrast",
  "logFC",
  "CI_low",
  "CI_high",
  "P_value",
  "FDR_KinHub",
  "FDR_target2"
)

print(
  target_primary[
    ,
    show_cols,
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\n========================================\n")
cat("PRIMARY VS SENSITIVITY\n")
cat("========================================\n")

print(
  concordance,
  row.names = FALSE
)

cat("\n========================================\n")
cat("MODEL INFORMATION\n")
cat("========================================\n")

cat(
  "Promoters used for lmFit/eBayes:",
  nrow(pm),
  "\n"
)

cat(
  "Kinase promoters used for KinHub BH:",
  sum(kinase_idx),
  "\n"
)

cat(
  "Transformation: log2(signal + 1)\n"
)

cat(
  "Design: unpaired four-group Day x Treatment model\n"
)

cat(
  "Primary scaling: peer-median correction of Day6 RPMI rep2\n"
)

cat(
  "Sensitivity scaling: Day6 RPMI rep2 matched to RPMI rep1\n"
)

cat("\n18b COMPLETE\n")

