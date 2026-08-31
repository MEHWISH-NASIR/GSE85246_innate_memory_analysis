# ============================================================
# 01_setup_and_data.R
#
# GSE85246 / GSE85243
#
# FINAL MEMORY-SAFE PROJECT SETUP
#
# Purpose:
#   1. Define all 25 RNA samples
#   2. Verify every MMSEQ input file
#   3. Verify identical gene IDs/order
#   4. Collect raw sample-level QC sequentially
#   5. Load canonical KinHub/OpenKinome kinase reference
#   6. Save lightweight setup objects
#
# IMPORTANT:
#   NO normalization
#   NO differential expression
#   NO kinase candidate selection
#
# Next:
#   02_initial_LPS_response.R
# ============================================================

rm(list = ls())

# ============================================================
# 1. DIRECTORIES
# ============================================================

outdir <- "results/01_setup"
figdir <- "figures/01_setup"

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

# ============================================================
# 2. SAMPLE MANIFEST
# ============================================================

samples <- data.frame(
  
  sample = c(
    
    # Day 1
    "LPS_d1_6094",
    "RPMI_d1_6092",
    "LPS_d1_6716",
    "RPMI_d1_6715",
    
    # Day 6 state
    "LPS_d6_6097",
    "RPMI_d6_6095",
    "LPS_d6_6719",
    "RPMI_d6_6718",
    "LPS_d6_8747",
    "RPMI_d6_8746",
    "LPS_d6_9244",
    "RPMI_d6_9243",
    "LPS_d6_10649",
    "RPMI_d6_10647",
    
    # BG rescue, non-restim
    "RESCUE_d6_10645",
    "RESCUE_d6_10653",
    
    # LPS restimulation
    "LPS_Restim_8750",
    "LPS_Restim_9247",
    "LPS_Restim_10650",
    
    # RPMI restimulation
    "RPMI_Restim_8749",
    "RPMI_Restim_8754",
    "RPMI_Restim_9246",
    "RPMI_Restim_10648",
    
    # BG rescue + restimulation
    "RESCUE_Restim_10646",
    "RESCUE_Restim_10654"
  ),
  
  GSM = c(
    
    "GSM2262878",
    "GSM2262899",
    "GSM2262879",
    "GSM2262900",
    
    "GSM2262880",
    "GSM2262901",
    "GSM2262881",
    "GSM2262902",
    "GSM2262882",
    "GSM2262903",
    "GSM2262883",
    "GSM2262905",
    "GSM2262886",
    "GSM2262909",
    
    "GSM2262913",
    "GSM2262914",
    
    "GSM2262884",
    "GSM2262885",
    "GSM2262887",
    
    "GSM2262906",
    "GSM2262907",
    "GSM2262908",
    "GSM2262910",
    
    "GSM2262911",
    "GSM2262912"
  ),
  
  analysis_set = c(
    
    rep("Day1_state", 4),
    rep("Day6_state", 10),
    rep("Day6_rescue", 2),
    rep("Day6_restim", 9)
  ),
  
  condition = c(
    
    "LPS", "RPMI",
    "LPS", "RPMI",
    
    "LPS", "RPMI",
    "LPS", "RPMI",
    "LPS", "RPMI",
    "LPS", "RPMI",
    "LPS", "RPMI",
    
    "RESCUE", "RESCUE",
    
    "LPS", "LPS", "LPS",
    
    "RPMI", "RPMI", "RPMI", "RPMI",
    
    "RESCUE", "RESCUE"
  ),
  
  state_pair = c(
    
    "d1_609x", "d1_609x",
    "d1_671x", "d1_671x",
    
    "d6_609x", "d6_609x",
    "d6_671x", "d6_671x",
    "d6_874x", "d6_874x",
    "d6_924x", "d6_924x",
    "d6_don70", "d6_don70",
    
    NA, NA,
    
    rep(NA, 9)
  ),
  
  response_pair = c(
    
    rep(NA, 8),
    
    # 874
    "resp_874_LPS",
    "resp_874_RPMI",
    
    # 924
    "resp_924_LPS",
    "resp_924_RPMI",
    
    # donor70
    "resp_don70_LPS",
    "resp_don70_RPMI",
    
    # rescue
    "resp_rescue1",
    "resp_rescue2",
    
    # corresponding restim samples
    "resp_874_LPS",
    "resp_924_LPS",
    "resp_don70_LPS",
    
    "resp_874_RPMI",
    NA,
    "resp_924_RPMI",
    "resp_don70_RPMI",
    
    "resp_rescue1",
    "resp_rescue2"
  ),
  
  file = c(
    
    # Day 1
    "data/processed/mmseq_state/GSM2262878_LPS_d1_6094.gene.mmseq.txt",
    "data/processed/mmseq_state/GSM2262899_RPMI_d1_6092.gene.mmseq.txt",
    "data/processed/mmseq_state/GSM2262879_LPS_d1_6716.gene.mmseq.txt",
    "data/processed/mmseq_state/GSM2262900_RPMI_d1_6715.gene.mmseq.txt",
    
    # Day 6
    "data/processed/mmseq_state/GSM2262880_LPS_d6_6097.gene.mmseq.txt",
    "data/processed/mmseq_state/GSM2262901_RPMI_d6_6095.gene.mmseq.txt",
    
    "data/processed/mmseq_state/GSM2262881_LPS_d6_6719.gene.mmseq.txt",
    "data/processed/mmseq_state/GSM2262902_RPMI_d6_6718.gene.mmseq.txt",
    
    "data/processed/mmseq_state/GSM2262882_LPS_d6_8747.gene.mmseq.txt",
    "data/processed/mmseq_state/GSM2262903_RPMI_d6_8746.gene.mmseq.txt",
    
    "data/processed/mmseq_state/GSM2262883_LPS_d6_9244.gene.mmseq.txt",
    "data/processed/mmseq_state/GSM2262905_RPMI_d6_9243.gene.mmseq.txt",
    
    "data/processed/mmseq_state/GSM2262886_LPS_d6_don70_10649.gene.mmseq.txt",
    "data/processed/mmseq_state/GSM2262909_RPMI_d6_don70_10647.gene.mmseq.txt",
    
    # Rescue
    "data/processed/mmseq_rescue/GSM2262913_Rescue_LPS_BGd2_d6_10645.gene.mmseq.txt",
    "data/processed/mmseq_rescue/GSM2262914_Rescue_LPS_BGd2_d6_10653.gene.mmseq.txt",
    
    # LPS restim
    "data/processed/mmseq_restim/GSM2262884_LPS_d6_Restim_8750.gene.mmseq.txt",
    "data/processed/mmseq_restim/GSM2262885_LPS_d6_Restim_9247.gene.mmseq.txt",
    "data/processed/mmseq_restim/GSM2262887_LPS_d6_Restim_don70_10650.gene.mmseq.txt",
    
    # RPMI restim
    "data/processed/mmseq_restim/GSM2262906_RPMI_d6_Restim_8749.gene.mmseq.txt",
    "data/processed/mmseq_restim/GSM2262907_RPMI_d6_Restim_8754.gene.mmseq.txt",
    "data/processed/mmseq_restim/GSM2262908_RPMI_d6_Restim_9246.gene.mmseq.txt",
    "data/processed/mmseq_restim/GSM2262910_RPMI_d6_Restim_don70_10648.gene.mmseq.txt",
    
    # Rescue restim
    "data/processed/mmseq_restim/GSM2262911_Rescue_LPS_BGd2_Restim_10646.gene.mmseq.txt",
    "data/processed/mmseq_restim/GSM2262912_Rescue_LPS_BGd2_Restim_10654.gene.mmseq.txt"
  ),
  
  stringsAsFactors = FALSE
)

# ============================================================
# 3. MANIFEST VALIDATION
# ============================================================

cat("\n========================================\n")
cat("GSE85246 RNA SAMPLE MANIFEST\n")
cat("========================================\n\n")

cat("Total samples:", nrow(samples), "\n\n")

print(
  table(samples$analysis_set)
)

if (nrow(samples) != 25) {
  stop("Expected exactly 25 RNA samples.")
}

if (anyDuplicated(samples$sample)) {
  stop("Duplicate sample names found.")
}

if (anyDuplicated(samples$GSM)) {
  stop("Duplicate GSM accessions found.")
}

missing <- samples[
  !file.exists(samples$file),
  c("sample", "GSM", "file")
]

if (nrow(missing) > 0) {
  
  cat("\nMissing files:\n")
  
  print(
    missing,
    row.names = FALSE
  )
  
  stop("One or more RNA files are missing.")
}

cat("\nAll 25 RNA MMSEQ files found successfully.\n")

write.csv(
  samples,
  file.path(
    outdir,
    "01_RNA_sample_manifest.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 4. KINHUB / OPENKINOME REFERENCE
# ============================================================

kinhub_file <-
  "data/reference/KinHubKinaseList.csv"

if (!file.exists(kinhub_file)) {
  
  stop(
    "KinHub reference file not found."
  )
}

kinhub <- read.csv(
  kinhub_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("KINHUB / OPENKINOME REFERENCE\n")
cat("========================================\n\n")

cat(
  "Reference rows:",
  nrow(kinhub),
  "\n"
)

raw_names <- colnames(kinhub)

clean_names <- gsub(
  "\u00A0",
  " ",
  raw_names,
  fixed = TRUE
)

clean_names <- gsub(
  "[[:space:]]+",
  " ",
  clean_names
)

clean_names <- trimws(
  clean_names
)

symbol_index <- match(
  "HGNC Name",
  clean_names
)

if (is.na(symbol_index)) {
  
  cat("\nAvailable columns:\n")
  print(raw_names)
  
  stop(
    "Cannot identify HGNC Name column."
  )
}

symbol_column <-
  raw_names[symbol_index]

cat(
  "Using symbol column:",
  symbol_column,
  "\n"
)

kinase_symbols <- as.character(
  kinhub[[symbol_column]]
)

kinase_symbols <- gsub(
  "\u00A0",
  " ",
  kinase_symbols,
  fixed = TRUE
)

kinase_symbols <- trimws(
  kinase_symbols
)

kinase_symbols <- unique(
  kinase_symbols[
    !is.na(kinase_symbols) &
      kinase_symbols != ""
  ]
)

cat(
  "Unique KinHub/OpenKinome symbols:",
  length(kinase_symbols),
  "\n"
)

write.csv(
  data.frame(
    SYMBOL = sort(kinase_symbols),
    stringsAsFactors = FALSE
  ),
  file.path(
    outdir,
    "01_KinHub_OpenKinome_symbols.csv"
  ),
  row.names = FALSE
)

# We do not need the complete reference table
# in memory after extracting the symbol list.
rm(kinhub)
gc()

# ============================================================
# 5. SEQUENTIAL MMSEQ AUDIT
#
# IMPORTANT:
# Read ONE file at a time.
# Do NOT retain all 25 full MMSEQ tables.
# ============================================================

required_columns <- c(
  "feature_id",
  "log_mu",
  "sd",
  "unique_hits",
  "observed"
)

reference_genes <- NULL

qc_list <- vector(
  "list",
  nrow(samples)
)

cat("\n========================================\n")
cat("SEQUENTIAL MMSEQ AUDIT\n")
cat("========================================\n\n")

for (i in seq_len(nrow(samples))) {
  
  cat(
    sprintf(
      "[%02d/%02d] %s ... ",
      i,
      nrow(samples),
      samples$sample[i]
    )
  )
  
  path <- samples$file[i]
  
  # ----------------------------------------------------------
  # Mapped-fragment count
  # ----------------------------------------------------------
  
  first_line <- readLines(
    path,
    n = 1,
    warn = FALSE
  )
  
  mapped_fragments <- suppressWarnings(
    as.numeric(
      sub(
        "^# Mapped fragments:[[:space:]]*",
        "",
        first_line
      )
    )
  )
  
  # ----------------------------------------------------------
  # Read this ONE MMSEQ table
  # ----------------------------------------------------------
  
  dat <- read.delim(
    path,
    skip = 1,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(dat)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      paste(
        "Missing MMSEQ columns in",
        basename(path),
        ":",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
  
  if (anyDuplicated(dat$feature_id)) {
    
    stop(
      paste(
        "Duplicate gene IDs in",
        basename(path)
      )
    )
  }
  
  # ----------------------------------------------------------
  # Establish / verify gene universe
  # ----------------------------------------------------------
  
  if (is.null(reference_genes)) {
    
    reference_genes <-
      dat$feature_id
    
  } else {
    
    if (!identical(
      dat$feature_id,
      reference_genes
    )) {
      
      stop(
        paste(
          "Gene ID/order mismatch in",
          basename(path)
        )
      )
    }
  }
  
  # ----------------------------------------------------------
  # Sample-level QC
  # ----------------------------------------------------------
  
  detected <- (
    dat$unique_hits > 0 &
      is.finite(dat$sd)
  )
  
  median_sd_detected <- if (
    any(detected)
  ) {
    
    median(
      dat$sd[detected],
      na.rm = TRUE
    )
    
  } else {
    
    NA_real_
  }
  
  qc_list[[i]] <- data.frame(
    
    sample =
      samples$sample[i],
    
    GSM =
      samples$GSM[i],
    
    analysis_set =
      samples$analysis_set[i],
    
    condition =
      samples$condition[i],
    
    mapped_fragments =
      mapped_fragments,
    
    genes_total =
      nrow(dat),
    
    genes_unique_hits_gt0 =
      sum(
        dat$unique_hits > 0,
        na.rm = TRUE
      ),
    
    genes_unique_hits_ge10 =
      sum(
        dat$unique_hits >= 10,
        na.rm = TRUE
      ),
    
    genes_observed_gt0 =
      sum(
        dat$observed > 0,
        na.rm = TRUE
      ),
    
    median_SD_detected =
      median_sd_detected,
    
    nonfinite_log_mu =
      sum(
        !is.finite(dat$log_mu)
      ),
    
    stringsAsFactors = FALSE
  )
  
  cat("OK\n")
  
  # Critical memory cleanup
  rm(
    dat,
    detected,
    first_line
  )
  
  gc(
    verbose = FALSE
  )
}

qc <- do.call(
  rbind,
  qc_list
)

rm(qc_list)

cat(
  "\nIdentical gene identity/order across all 25 samples: YES\n"
)

cat(
  "Genes represented per file:",
  length(reference_genes),
  "\n"
)

write.csv(
  qc,
  file.path(
    outdir,
    "01_RNA_sample_QC.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 6. STATE PAIR VALIDATION
# ============================================================

cat("\n========================================\n")
cat("STATE PAIR VALIDATION\n")
cat("========================================\n\n")

day1 <- samples[
  samples$analysis_set == "Day1_state",
]

day6 <- samples[
  samples$analysis_set == "Day6_state",
]

cat("Day 1:\n")

print(
  table(
    day1$state_pair,
    day1$condition
  )
)

cat("\nDay 6:\n")

print(
  table(
    day6$state_pair,
    day6$condition
  )
)

day1_tab <- table(
  day1$state_pair,
  day1$condition
)

day6_tab <- table(
  day6$state_pair,
  day6$condition
)

if (!all(
  day1_tab[, c("LPS", "RPMI")] == 1
)) {
  
  stop(
    "Day-1 pair structure is incorrect."
  )
}

if (!all(
  day6_tab[, c("LPS", "RPMI")] == 1
)) {
  
  stop(
    "Day-6 pair structure is incorrect."
  )
}

cat(
  "\nDay-1 and Day-6 pair structure: VALID\n"
)

# ============================================================
# 7. BASELINE / RESTIMULATION PAIR VALIDATION
# ============================================================

cat("\n========================================\n")
cat("BASELINE / RESTIMULATION PAIRS\n")
cat("========================================\n\n")

paired <- samples[
  !is.na(samples$response_pair),
]

pair_counts <- table(
  paired$response_pair
)

print(
  pair_counts
)

if (!all(pair_counts == 2)) {
  
  cat(
    "\nIncorrect response-pair counts:\n"
  )
  
  print(
    pair_counts[
      pair_counts != 2
    ]
  )
  
  stop(
    "Response-pair structure is incorrect."
  )
}

cat(
  "\nAll defined baseline/restimulation pairs have 2 samples: YES\n"
)

# ============================================================
# 8. UNMATCHED RESTIMULATION SAMPLE
# ============================================================

unmatched <- samples[
  samples$analysis_set == "Day6_restim" &
    is.na(samples$response_pair),
  c(
    "sample",
    "GSM",
    "condition",
    "file"
  )
]

cat("\n========================================\n")
cat("UNMATCHED RESTIMULATION SAMPLES\n")
cat("========================================\n\n")

if (nrow(unmatched) == 0) {
  
  cat("NONE\n")
  
} else {
  
  print(
    unmatched[
      ,
      c(
        "sample",
        "GSM",
        "condition"
      )
    ],
    row.names = FALSE
  )
}

write.csv(
  unmatched,
  file.path(
    outdir,
    "01_unmatched_restimulation_samples.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 9. BASIC QC FIGURES
# ============================================================

png(
  file.path(
    figdir,
    "01_mapped_fragments.png"
  ),
  width = 2400,
  height = 1300,
  res = 180
)

par(
  mar = c(
    11,
    5,
    4,
    1
  )
)

barplot(
  qc$mapped_fragments / 1e6,
  names.arg = qc$sample,
  las = 2,
  ylab = "Mapped fragments (millions)",
  main = "GSE85246 RNA-seq mapped fragments"
)

dev.off()

png(
  file.path(
    figdir,
    "01_detected_genes.png"
  ),
  width = 2400,
  height = 1300,
  res = 180
)

par(
  mar = c(
    11,
    5,
    4,
    1
  )
)

barplot(
  qc$genes_unique_hits_gt0,
  names.arg = qc$sample,
  las = 2,
  ylab = "Genes with >=1 unique hit",
  main = "GSE85246 detected genes per sample"
)

dev.off()

# ============================================================
# 10. SAVE LIGHTWEIGHT SETUP OBJECT
#
# NO large expression matrices are saved here.
# ============================================================

setup_object <- list(
  
  samples =
    samples,
  
  kinase_symbols =
    kinase_symbols,
  
  kinase_symbol_column =
    symbol_column,
  
  gene_ids =
    reference_genes,
  
  sample_QC =
    qc,
  
  unmatched_restim =
    unmatched
)

saveRDS(
  setup_object,
  file.path(
    outdir,
    "01_setup_object.rds"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    outdir,
    "01_sessionInfo.txt"
  )
)

# ============================================================
# 11. FINAL SUMMARY
# ============================================================

cat("\n========================================\n")
cat("01 SETUP + DATA COMPLETE\n")
cat("========================================\n\n")

cat(
  "RNA samples:",
  nrow(samples),
  "\n"
)

cat(
  "Day-1 state samples:",
  sum(
    samples$analysis_set ==
      "Day1_state"
  ),
  "\n"
)

cat(
  "Day-6 state samples:",
  sum(
    samples$analysis_set ==
      "Day6_state"
  ),
  "\n"
)

cat(
  "Day-6 rescue samples:",
  sum(
    samples$analysis_set ==
      "Day6_rescue"
  ),
  "\n"
)

cat(
  "Day-6 restimulation samples:",
  sum(
    samples$analysis_set ==
      "Day6_restim"
  ),
  "\n"
)

cat(
  "MMSEQ genes:",
  length(reference_genes),
  "\n"
)

cat(
  "KinHub/OpenKinome symbols:",
  length(kinase_symbols),
  "\n"
)

cat(
  "Unmatched restimulation samples:",
  nrow(unmatched),
  "\n"
)

cat("\n----------------------------------------\n")

cat(
  "NO normalization performed.\n"
)

cat(
  "NO differential-expression analysis performed.\n"
)

cat(
  "NO kinase candidates selected.\n"
)

cat("----------------------------------------\n")

cat(
  "\nNext canonical step:\n"
)

cat(
  "scripts/02_initial_LPS_response.R\n"
)

cat("\n========================================\n")