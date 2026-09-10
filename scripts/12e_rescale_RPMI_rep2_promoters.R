# ============================================================
# 12e — RESCALE THREE NOT-NORMALIZED RPMI_d6_rep2 TRACKS
#
# Input:
#   raw promoter matrix from Step 12d
#   genome-wide global-signal QC from Step 12c
#
# Two correction schemes are produced:
#
# PRIMARY:
#   Match each NotNormalized RPMI_rep2 track to the
#   median genome-wide mean of all NORMALIZED Day-6
#   tracks for the same histone mark.
#
# SENSITIVITY:
#   Match RPMI_rep2 directly to RPMI_rep1 for the same mark.
#
# No statistical testing is done here.
# ============================================================

rm(list = ls())

outdir <- "results/12_kinase_reevaluation"

raw_file <- file.path(
  outdir,
  "12d_day6_promoter_signal_raw.rds"
)

qc_file <- file.path(
  outdir,
  "12c_bigwig_global_qc.csv"
)

if (!file.exists(raw_file)) {
  stop("Missing Step 12d raw promoter matrix.")
}

if (!file.exists(qc_file)) {
  stop("Missing Step 12c BigWig QC table.")
}

x <- readRDS(raw_file)

signal_raw <- x$signal_matrix_raw
promoter_manifest <- x$promoter_manifest
sample_manifest <- x$sample_manifest

qc <- read.csv(
  qc_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

marks <- c(
  "H3K27ac",
  "H3K4me1",
  "H3K4me3"
)

cat("\n========================================\n")
cat("12e — RPMI REP2 RESCALING\n")
cat("========================================\n\n")

signal_primary <- signal_raw
signal_sensitivity <- signal_raw

factor_table <- list()

for (m in marks) {

  q <- qc[
    qc$mark == m,
  ]

  normalized_means <- q$global_mean[
    q$normalization_label == "Normalized"
  ]

  bad_mean <- q$global_mean[
    q$normalization_label == "NotNormalized"
  ]

  if (
    length(normalized_means) != 5 ||
    length(bad_mean) != 1
  ) {
    stop(
      paste(
        "Unexpected QC structure for",
        m
      )
    )
  }

  # ----------------------------------------------------------
  # PRIMARY:
  # median of all normalized Day-6 peers of same mark
  # ----------------------------------------------------------

  peer_median <- median(
    normalized_means,
    na.rm = TRUE
  )

  factor_primary <-
    peer_median / bad_mean

  # ----------------------------------------------------------
  # SENSITIVITY:
  # RPMI_rep1 / RPMI_rep2 genome-wide mean
  # ----------------------------------------------------------

  rpm1 <- q$global_mean[
    grepl(
      "RPMI_rep1",
      q$sample
    )
  ]

  if (length(rpm1) != 1) {
    stop(
      paste(
        "Could not identify RPMI_rep1 for",
        m
      )
    )
  }

  factor_sensitivity <-
    rpm1 / bad_mean

  sample_id <- paste(
    m,
    "RPMI",
    "rep2",
    sep = "_"
  )

  if (!(sample_id %in% colnames(signal_raw))) {
    stop(
      paste(
        "Missing promoter column:",
        sample_id
      )
    )
  }

  signal_primary[, sample_id] <-
    signal_raw[, sample_id] *
    factor_primary

  signal_sensitivity[, sample_id] <-
    signal_raw[, sample_id] *
    factor_sensitivity

  factor_table[[m]] <- data.frame(
    mark = m,
    deposited_RPMI_rep2_global_mean =
      bad_mean,
    normalized_peer_median =
      peer_median,
    primary_scaling_factor =
      factor_primary,
    RPMI_rep1_global_mean =
      rpm1,
    sensitivity_scaling_factor =
      factor_sensitivity,
    stringsAsFactors = FALSE
  )
}

factor_table <- do.call(
  rbind,
  factor_table
)

cat("Scaling factors:\n\n")

print(
  factor_table,
  row.names = FALSE
)

write.csv(
  factor_table,
  file.path(
    outdir,
    "12e_RPMI_rep2_scaling_factors.csv"
  ),
  row.names = FALSE
)

# ============================================================
# SAVE BOTH CORRECTED MATRICES
# ============================================================

saveRDS(
  list(
    promoters = x$promoters,
    promoter_manifest =
      promoter_manifest,
    sample_manifest =
      sample_manifest,
    signal_matrix_raw =
      signal_raw,
    signal_matrix_primary =
      signal_primary,
    signal_matrix_sensitivity =
      signal_sensitivity,
    scaling_factors =
      factor_table
  ),
  file.path(
    outdir,
    "12e_day6_promoter_signal_rescaled.rds"
  )
)

# ============================================================
# TARGET-GENE AUDIT
# ============================================================

targets <- c(
  "JAK3",
  "EPHB2",
  "MET",
  "MAP3K8",
  "BMPR1A"
)

target_idx <- which(
  promoter_manifest$SYMBOL %in%
    targets
)

make_target_table <- function(mat, method) {

  z <- cbind(
    promoter_manifest[target_idx, ],
    as.data.frame(
      mat[target_idx, , drop = FALSE],
      check.names = FALSE
    )
  )

  z$Scaling_method <- method

  z
}

target_raw <- make_target_table(
  signal_raw,
  "As_deposited"
)

target_primary <- make_target_table(
  signal_primary,
  "Peer_median_corrected"
)

target_sensitivity <- make_target_table(
  signal_sensitivity,
  "RPMI_rep1_matched"
)

target_all <- rbind(
  target_raw,
  target_primary,
  target_sensitivity
)

write.csv(
  target_all,
  file.path(
    outdir,
    "12e_target_promoter_signal_scaling_audit.csv"
  ),
  row.names = FALSE
)

# ============================================================
# POST-CORRECTION COLUMN MEDIANS
# ============================================================

cat("\n========================================\n")
cat("RAW COLUMN MEDIANS\n")
cat("========================================\n\n")

print(
  apply(
    signal_raw,
    2,
    median,
    na.rm = TRUE
  )
)

cat("\n========================================\n")
cat("PRIMARY-CORRECTED COLUMN MEDIANS\n")
cat("========================================\n\n")

print(
  apply(
    signal_primary,
    2,
    median,
    na.rm = TRUE
  )
)

cat("\n========================================\n")
cat("TARGET SIGNALS — PRIMARY CORRECTION\n")
cat("========================================\n\n")

print(
  target_primary,
  row.names = FALSE
)

cat("\nSaved:\n")
cat(
  " ",
  file.path(
    outdir,
    "12e_RPMI_rep2_scaling_factors.csv"
  ),
  "\n"
)
cat(
  " ",
  file.path(
    outdir,
    "12e_day6_promoter_signal_rescaled.rds"
  ),
  "\n"
)
cat(
  " ",
  file.path(
    outdir,
    "12e_target_promoter_signal_scaling_audit.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("12e COMPLETE\n")
cat("========================================\n")
