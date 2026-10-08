
rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
})

setDTthreads(1)

indir <- "macrophage_transcriptomics_validation/GSE168468/data/processed"
outdir <- "macrophage_transcriptomics_validation/GSE168468/results"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  indir,
  pattern = "gene\\.mmseq\\.txt\\.gz$",
  full.names = TRUE
)

# Restrict to resting Day6 samples only:
keep <- grepl(
  "D6_RPMI|d6R|d6L|D6_Den|d6D|D6_Bul|d6B",
  basename(files)
)

# Exclude restimulated d7 / D6_4h samples
keep <- keep & !grepl(
  "d7|D6_4h",
  basename(files)
)

files <- files[keep]

cat("Day6 files found:", length(files), "\n")

if (length(files) != 20) {
  stop("Expected exactly 20 resting Day6 samples.")
}

res <- vector("list", length(files))

for (i in seq_along(files)) {

  f <- files[i]

  x <- fread(
    f,
    skip = 1,
    header = TRUE,
    nThread = 1
  )

  # Use genes that are observed and finite
  y <- x[
    observed == 1 &
    is.finite(log_mu),
    log_mu
  ]

  res[[i]] <- data.frame(
    sample = basename(f),
    n_observed = length(y),
    median_log_mu = median(y),
    mean_log_mu = mean(y),
    q25 = unname(quantile(y, 0.25)),
    q75 = unname(quantile(y, 0.75)),
    stringsAsFactors = FALSE
  )
}

audit <- rbindlist(res)

global_median <- median(audit$median_log_mu)

audit[, median_center_offset :=
        median_log_mu - global_median]

audit[, scaling_adjustment :=
        -median_center_offset]

fwrite(
  audit,
  file.path(
    outdir,
    "20b_GSE168468_Day6_MMSEQ_scaling_audit.csv"
  )
)

cat("\n===== DAY6 MMSEQ SCALING AUDIT =====\n")
print(audit)

cat("\nGlobal median of sample medians:",
    global_median, "\n")

cat("\nRange of median offsets:",
    range(audit$median_center_offset), "\n")

cat("\n20b COMPLETE\n")

