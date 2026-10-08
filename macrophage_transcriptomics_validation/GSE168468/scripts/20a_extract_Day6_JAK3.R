
rm(list = ls())

library(data.table)
setDTthreads(1)

indir <- "macrophage_transcriptomics_validation/GSE168468/data/processed"
outdir <- "macrophage_transcriptomics_validation/GSE168468/results"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

jak3 <- "ENSG00000105639"

samples <- data.frame(
  donor = rep(paste0("A", 1:5), 4),

  condition = rep(
    c("RPMI", "LPS", "BCG_Denmark", "BCG_Bulgaria"),
    each = 5
  ),

  file = c(
    # Day6 RPMI
    "GSM5141114_A1_D6_RPMI.gene.mmseq.txt.gz",
    "GSM5141115_A2_d6R.gene.mmseq.txt.gz",
    "GSM5141116_A3_d6R.gene.mmseq.txt.gz",
    "GSM5141117_A4_d6R.gene.mmseq.txt.gz",
    "GSM5141118_A5_d6R.gene.mmseq.txt.gz",

    # Day6 LPS
    "GSM5141119_A1_d6L.gene.mmseq.txt.gz",
    "GSM5141120_A2_d6L.gene.mmseq.txt.gz",
    "GSM5141121_A3_d6L.gene.mmseq.txt.gz",
    "GSM5141122_A4_d6L.gene.mmseq.txt.gz",
    "GSM5141123_A5_d6L.gene.mmseq.txt.gz",

    # Day6 BCG Denmark
    "GSM5141124_A1_D6_Den.gene.mmseq.txt.gz",
    "GSM5141125_A2_d6D.gene.mmseq.txt.gz",
    "GSM5141126_A3_d6D.gene.mmseq.txt.gz",
    "GSM5141127_A4_d6D.gene.mmseq.txt.gz",
    "GSM5141128_A5_d6D.gene.mmseq.txt.gz",

    # Day6 BCG Bulgaria
    "GSM5141129_A1_D6_Bul.gene.mmseq.txt.gz",
    "GSM5141130_A2_d6B.gene.mmseq.txt.gz",
    "GSM5141131_A3_d6B.gene.mmseq.txt.gz",
    "GSM5141132_A4_d6B.gene.mmseq.txt.gz",
    "GSM5141133_A5_d6B.gene.mmseq.txt.gz"
  ),

  stringsAsFactors = FALSE
)

res <- vector("list", nrow(samples))

for (i in seq_len(nrow(samples))) {

  f <- file.path(indir, samples$file[i])

  if (!file.exists(f)) {
    stop("Missing file: ", f)
  }

  x <- fread(
    f,
    skip = 1,
    header = TRUE,
    nThread = 1
  )

  z <- x[feature_id == jak3]

  if (nrow(z) != 1) {
    stop("Expected one JAK3 row in: ", samples$file[i])
  }

  res[[i]] <- data.frame(
    donor = samples$donor[i],
    condition = samples$condition[i],
    file = samples$file[i],
    log_mu = z$log_mu,
    sd = z$sd,
    unique_hits = z$unique_hits,
    observed = z$observed,
    stringsAsFactors = FALSE
  )
}

res <- rbindlist(res)

fwrite(
  res,
  file.path(outdir, "20a_GSE168468_Day6_JAK3_raw_logmu.csv")
)

cat("\n===== GSE168468 DAY6 JAK3 =====\n")
print(res[, .(
  donor,
  condition,
  log_mu,
  sd,
  unique_hits
)])

wide <- dcast(
  res,
  donor ~ condition,
  value.var = "log_mu"
)

wide[, LPS_minus_RPMI := LPS - RPMI]
wide[, Denmark_minus_RPMI := BCG_Denmark - RPMI]
wide[, Bulgaria_minus_RPMI := BCG_Bulgaria - RPMI]

fwrite(
  wide,
  file.path(outdir, "20a_GSE168468_Day6_JAK3_paired_raw_effects.csv")
)

cat("\n===== RAW log_mu PAIRED DIFFERENCES =====\n")
print(wide)

cat("\n20a COMPLETE\n")
