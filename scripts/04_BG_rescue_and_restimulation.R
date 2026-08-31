# ============================================================
# 04_BG_rescue_and_restimulation.R
#
# GSE85246 / GSE85243
#
# BETA-GLUCAN RESCUE + LPS RESTIMULATION
#
# Questions:
#
# A. BASELINE RESCUE
#    Does beta-glucan move the persistent Day-6
#    LPS-conditioned state toward RPMI?
#
# B. RESTIMULATION / TOLERANCE
#    Does prior LPS exposure alter the response to
#    a second LPS stimulation?
#
# C. RESPONSE RESTORATION
#    Does beta-glucan rescue move the restimulation
#    response toward the RPMI/naive response?
#
#
# Formal paired tolerance design:
#
# donor 874:
#   RPMI_d6_8746  -> RPMI_Restim_8749
#   LPS_d6_8747   -> LPS_Restim_8750
#
# donor 924:
#   RPMI_d6_9243  -> RPMI_Restim_9246
#   LPS_d6_9244   -> LPS_Restim_9247
#
# donor 70:
#   RPMI_d6_10647 -> RPMI_Restim_10648
#   LPS_d6_10649  -> LPS_Restim_10650
#
#
# Formal interaction:
#
#   (LPS_Restim - LPS_baseline)
#                -
#   (RPMI_Restim - RPMI_baseline)
#
#
# Beta-glucan rescue pairs:
#
#   RESCUE_d6_10645 -> RESCUE_Restim_10646
#   RESCUE_d6_10653 -> RESCUE_Restim_10654
#
# Rescue n = 2 and is therefore treated primarily as
# descriptive/directional evidence.
#
#
# IMPORTANT:
#
# - All samples used here are normalized JOINTLY.
#
# - RPMI_Restim_8754 is excluded from formal analysis
#   because it lacks a matched non-restimulated baseline.
#
# - Primary candidates come from Step 03.
#
# Next:
#   05_integrated_kinase_trajectory.R
# ============================================================


rm(list = ls())


# ============================================================
# 1. REQUIRED PACKAGES / FILES
# ============================================================

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Required package 'limma' is not installed.")
}

if (!file.exists("scripts/mmseq.R")) {
  stop("Missing scripts/mmseq.R")
}

if (!file.exists(
  "results/01_setup/01_setup_object.rds"
)) {
  stop(
    "Step 01 output missing. Run 01_setup_and_data.R first."
  )
}

if (!file.exists(
  "results/03_memory_kinases/03_memory_kinases_FDR05.csv"
)) {
  stop(
    "Step 03 memory-kinase table is missing."
  )
}

if (!file.exists(
  "results/03_memory_kinases/03_Day6_LPS_vs_RPMI_all_genes.csv"
)) {
  stop(
    "Step 03 genome-wide Day-6 result is missing."
  )
}


source("scripts/mmseq.R")


# ============================================================
# 2. OUTPUT DIRECTORIES
# ============================================================

outdir <- "results/04_BG_rescue_restimulation"
figdir <- "figures/04_BG_rescue_restimulation"

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
# 3. LOAD STEP-01 / STEP-03 INFORMATION
# ============================================================

setup <- readRDS(
  "results/01_setup/01_setup_object.rds"
)

project_samples <- setup$samples


memory24 <- read.csv(
  "results/03_memory_kinases/03_memory_kinases_FDR05.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


day6_annotation <- read.csv(
  "results/03_memory_kinases/03_Day6_LPS_vs_RPMI_all_genes.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("\n========================================\n")
cat("04 — BG RESCUE + RESTIMULATION\n")
cat("========================================\n\n")


cat(
  "Step-03 primary memory kinase candidates:",
  nrow(memory24),
  "\n"
)


if (anyDuplicated(memory24$ENSEMBL)) {
  stop("Duplicate ENSEMBL IDs found in Step-03 memory candidates.")
}


# ============================================================
# 4. DEFINE THE 16 ANALYSIS SAMPLES
#
# 12 samples:
#   formal RPMI/LPS response interaction
#
# 4 samples:
#   beta-glucan rescue response
#
# The unmatched RPMI_Restim_8754 is intentionally omitted.
# ============================================================

analysis_info <- data.frame(
  
  sample = c(
    
    # --------------------------------------------------------
    # donor 874
    # --------------------------------------------------------
    "RPMI_d6_8746",
    "RPMI_Restim_8749",
    "LPS_d6_8747",
    "LPS_Restim_8750",
    
    # --------------------------------------------------------
    # donor 924
    # --------------------------------------------------------
    "RPMI_d6_9243",
    "RPMI_Restim_9246",
    "LPS_d6_9244",
    "LPS_Restim_9247",
    
    # --------------------------------------------------------
    # donor 70
    # --------------------------------------------------------
    "RPMI_d6_10647",
    "RPMI_Restim_10648",
    "LPS_d6_10649",
    "LPS_Restim_10650",
    
    # --------------------------------------------------------
    # beta-glucan rescue pairs
    # --------------------------------------------------------
    "RESCUE_d6_10645",
    "RESCUE_Restim_10646",
    "RESCUE_d6_10653",
    "RESCUE_Restim_10654"
  ),
  
  arm = c(
    
    "RPMI", "RPMI", "LPS", "LPS",
    "RPMI", "RPMI", "LPS", "LPS",
    "RPMI", "RPMI", "LPS", "LPS",
    
    "RESCUE", "RESCUE",
    "RESCUE", "RESCUE"
  ),
  
  phase = c(
    
    "Baseline", "Restim",
    "Baseline", "Restim",
    
    "Baseline", "Restim",
    "Baseline", "Restim",
    
    "Baseline", "Restim",
    "Baseline", "Restim",
    
    "Baseline", "Restim",
    "Baseline", "Restim"
  ),
  
  donor_block = c(
    
    "874", "874", "874", "874",
    "924", "924", "924", "924",
    "70",  "70",  "70",  "70",
    
    "rescue1", "rescue1",
    "rescue2", "rescue2"
  ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# 5. LINK ANALYSIS MANIFEST TO STEP-01 FILES
# ============================================================

manifest_idx <- match(
  analysis_info$sample,
  project_samples$sample
)


if (anyNA(manifest_idx)) {
  
  cat("\nMissing samples:\n")
  
  print(
    analysis_info$sample[
      is.na(manifest_idx)
    ]
  )
  
  stop(
    "One or more Step-04 samples are absent from Step-01 manifest."
  )
}


analysis_info$GSM <-
  project_samples$GSM[
    manifest_idx
  ]


analysis_info$file <-
  project_samples$file[
    manifest_idx
  ]


if (!all(file.exists(analysis_info$file))) {
  
  stop(
    "One or more Step-04 MMSEQ files are missing."
  )
}


write.csv(
  analysis_info,
  file.path(
    outdir,
    "04_analysis_sample_manifest.csv"
  ),
  row.names = FALSE
)


cat("\nAnalysis samples:\n\n")

print(
  analysis_info[
    ,
    c(
      "sample",
      "arm",
      "phase",
      "donor_block"
    )
  ],
  row.names = FALSE
)


cat(
  "\nUnmatched RPMI_Restim_8754 is NOT used in formal Step 04.\n"
)


# ============================================================
# 6. READ + JOINTLY NORMALIZE ALL 16 ANALYSIS SAMPLES
#
# Critical:
# Baseline and restimulation samples are normalized together
# before response differences are calculated.
# ============================================================

cat("\n========================================\n")
cat("JOINT MMSEQ NORMALIZATION\n")
cat("========================================\n\n")


mm <- readmmseq(
  
  mmseq_files =
    analysis_info$file,
  
  sample_names =
    analysis_info$sample,
  
  normalize = TRUE
)


required_mmseq_objects <- c(
  "log_mu",
  "unique_hits"
)


missing_mmseq <- setdiff(
  required_mmseq_objects,
  names(mm)
)


if (length(missing_mmseq) > 0) {
  
  stop(
    paste(
      "readmmseq output missing:",
      paste(
        missing_mmseq,
        collapse = ", "
      )
    )
  )
}


cat(
  "\nGenes in joint MMSEQ object:",
  nrow(mm$log_mu),
  "\n"
)

cat(
  "Samples in joint MMSEQ object:",
  ncol(mm$log_mu),
  "\n"
)


# ============================================================
# 7. NATURAL LOG -> LOG2
# ============================================================

expr <- mm$log_mu / log(2)

colnames(expr) <-
  analysis_info$sample


# ============================================================
# 8. EXPRESSION FILTER
#
# Same general philosophy as Steps 02/03:
# retain genes with sequencing support in >=2 samples.
# ============================================================

keep <- rowSums(
  mm$unique_hits >= 1,
  na.rm = TRUE
) >= 2


cat(
  "\nGenes before filtering:",
  nrow(expr),
  "\n"
)

cat(
  "Genes retained:",
  sum(keep),
  "\n"
)

cat(
  "Genes removed:",
  sum(!keep),
  "\n"
)


expr_f <- expr[
  keep,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# Zero-variance removal
# ------------------------------------------------------------

gene_var <- apply(
  expr_f,
  1,
  var,
  na.rm = TRUE
)


keep_var <-
  is.finite(gene_var) &
  gene_var > 0


expr_f <- expr_f[
  keep_var,
  ,
  drop = FALSE
]


cat(
  "Genes after zero-variance removal:",
  nrow(expr_f),
  "\n"
)


# ============================================================
# 9. VERIFY STEP-03 MEMORY CANDIDATES ARE AVAILABLE
# ============================================================

candidate_present <-
  memory24$ENSEMBL %in%
  rownames(expr_f)


cat(
  "\nStep-03 memory candidates retained in Step 04:",
  sum(candidate_present),
  "/",
  nrow(memory24),
  "\n"
)


if (any(!candidate_present)) {
  
  cat(
    "\nCandidates not evaluable after Step-04 filtering:\n"
  )
  
  print(
    memory24[
      !candidate_present,
      c(
        "ENSEMBL",
        "SYMBOL"
      )
    ],
    row.names = FALSE
  )
}


# ============================================================
# 10. FORMAL 3-DONOR TOLERANCE INTERACTION DATA
#
# Use only:
#   RPMI + LPS
#   donors 874, 924, 70
#
# Rescue samples are NOT part of this formal interaction.
# ============================================================

formal_info <- analysis_info[
  analysis_info$arm %in%
    c(
      "RPMI",
      "LPS"
    ),
  ,
  drop = FALSE
]


formal_samples <-
  formal_info$sample


expr_formal <- expr_f[
  ,
  formal_samples,
  drop = FALSE
]


if (ncol(expr_formal) != 12) {
  
  stop(
    "Expected exactly 12 samples in formal tolerance model."
  )
}


# ============================================================
# 11. FORMAL FACTORS
# ============================================================

formal_info$donor <- factor(
  formal_info$donor_block,
  levels = c(
    "874",
    "924",
    "70"
  )
)


formal_info$history <- factor(
  formal_info$arm,
  levels = c(
    "RPMI",
    "LPS"
  )
)


formal_info$restim <- factor(
  formal_info$phase,
  levels = c(
    "Baseline",
    "Restim"
  )
)


# ============================================================
# 12. FORMAL INTERACTION DESIGN
#
# Model:
#
# expression ~ donor + history * restim
#
# historyLPS:
#   baseline LPS-history vs RPMI-history difference
#
# restimRestim:
#   RPMI/naive restimulation response
#
# historyLPS:restimRestim:
#   difference in restimulation response caused by
#   prior LPS exposure
# ============================================================

design <- model.matrix(
  ~ donor + history * restim,
  data = formal_info
)


rownames(design) <-
  formal_info$sample


# Make coefficient names safe for contrasts
colnames(design) <-
  make.names(
    colnames(design)
  )


cat("\n========================================\n")
cat("FORMAL TOLERANCE INTERACTION DESIGN\n")
cat("========================================\n\n")


print(design)


cat(
  "\nDesign rank:",
  qr(design)$rank,
  "/",
  ncol(design),
  "\n"
)


if (qr(design)$rank != ncol(design)) {
  
  stop(
    "Formal tolerance design matrix is not full rank."
  )
}


required_coef <- c(
  "restimRestim",
  "historyLPS.restimRestim"
)


if (!all(
  required_coef %in%
  colnames(design)
)) {
  
  cat(
    "\nDesign coefficient names:\n"
  )
  
  print(
    colnames(design)
  )
  
  stop(
    "Required interaction coefficients not found."
  )
}


# ============================================================
# 13. FIT FORMAL LIMMA MODEL
# ============================================================

fit0 <- limma::lmFit(
  expr_formal,
  design
)


contrast_matrix <- limma::makeContrasts(
  
  Naive_response =
    restimRestim,
  
  Tolerant_response =
    restimRestim +
    historyLPS.restimRestim,
  
  Tolerance_interaction =
    historyLPS.restimRestim,
  
  levels = design
)


fit <- limma::contrasts.fit(
  fit0,
  contrast_matrix
)


fit <- limma::eBayes(
  fit,
  trend = TRUE,
  robust = TRUE
)


# ============================================================
# 14. BUILD COMPLETE FORMAL INTERACTION TABLE
# ============================================================

formal_results <- data.frame(
  
  ENSEMBL =
    rownames(expr_formal),
  
  Naive_response_log2 =
    fit$coefficients[
      ,
      "Naive_response"
    ],
  
  Tolerant_response_log2 =
    fit$coefficients[
      ,
      "Tolerant_response"
    ],
  
  Tolerance_interaction_log2 =
    fit$coefficients[
      ,
      "Tolerance_interaction"
    ],
  
  Interaction_t =
    fit$t[
      ,
      "Tolerance_interaction"
    ],
  
  Interaction_P =
    fit$p.value[
      ,
      "Tolerance_interaction"
    ],
  
  stringsAsFactors = FALSE
)


formal_results$Interaction_FDR_genome <-
  p.adjust(
    formal_results$Interaction_P,
    method = "BH"
  )


# ============================================================
# 15. ADD AVAILABLE GENE ANNOTATION
# ============================================================

annot_idx <- match(
  formal_results$ENSEMBL,
  day6_annotation$ENSEMBL
)


formal_results$SYMBOL <-
  day6_annotation$SYMBOL[
    annot_idx
  ]


formal_results$GENENAME <-
  day6_annotation$GENENAME[
    annot_idx
  ]


formal_results <- formal_results[
  ,
  c(
    "ENSEMBL",
    "SYMBOL",
    "GENENAME",
    "Naive_response_log2",
    "Tolerant_response_log2",
    "Tolerance_interaction_log2",
    "Interaction_t",
    "Interaction_P",
    "Interaction_FDR_genome"
  )
]


formal_results <- formal_results[
  order(
    formal_results$Interaction_P
  ),
  ,
  drop = FALSE
]


write.csv(
  formal_results,
  file.path(
    outdir,
    "04_formal_tolerance_interaction_all_genes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. MANUAL PER-DONOR RESPONSE CALCULATION
#
# Used for transparent directional summaries.
# ============================================================

response_matrix <- cbind(
  
  RPMI_874 =
    expr_f[, "RPMI_Restim_8749"] -
    expr_f[, "RPMI_d6_8746"],
  
  RPMI_924 =
    expr_f[, "RPMI_Restim_9246"] -
    expr_f[, "RPMI_d6_9243"],
  
  RPMI_70 =
    expr_f[, "RPMI_Restim_10648"] -
    expr_f[, "RPMI_d6_10647"],
  
  
  LPS_874 =
    expr_f[, "LPS_Restim_8750"] -
    expr_f[, "LPS_d6_8747"],
  
  LPS_924 =
    expr_f[, "LPS_Restim_9247"] -
    expr_f[, "LPS_d6_9244"],
  
  LPS_70 =
    expr_f[, "LPS_Restim_10650"] -
    expr_f[, "LPS_d6_10649"],
  
  
  RESCUE_1 =
    expr_f[, "RESCUE_Restim_10646"] -
    expr_f[, "RESCUE_d6_10645"],
  
  RESCUE_2 =
    expr_f[, "RESCUE_Restim_10654"] -
    expr_f[, "RESCUE_d6_10653"]
)


naive_response_mean <- rowMeans(
  response_matrix[
    ,
    c(
      "RPMI_874",
      "RPMI_924",
      "RPMI_70"
    ),
    drop = FALSE
  ],
  na.rm = TRUE
)


tolerant_response_mean <- rowMeans(
  response_matrix[
    ,
    c(
      "LPS_874",
      "LPS_924",
      "LPS_70"
    ),
    drop = FALSE
  ],
  na.rm = TRUE
)


rescue_response_mean <- rowMeans(
  response_matrix[
    ,
    c(
      "RESCUE_1",
      "RESCUE_2"
    ),
    drop = FALSE
  ],
  na.rm = TRUE
)


manual_interaction <-
  tolerant_response_mean -
  naive_response_mean


# ============================================================
# 17. VERIFY MANUAL INTERACTION == LIMMA COEFFICIENT
# ============================================================

formal_ordered <- formal_results[
  match(
    rownames(expr_f),
    formal_results$ENSEMBL
  ),
  ,
  drop = FALSE
]


difference_check <- max(
  abs(
    manual_interaction -
      formal_ordered$Tolerance_interaction_log2
  ),
  na.rm = TRUE
)


cat(
  "\nMaximum manual-vs-model interaction difference:",
  signif(
    difference_check,
    6
  ),
  "\n"
)


if (
  is.finite(difference_check) &&
  difference_check > 1e-6
) {
  
  warning(
    "Manual and model interaction estimates differ unexpectedly."
  )
}


# ============================================================
# 18. BASELINE BG RESCUE METRICS
#
# Descriptive only.
#
# Compare:
#   mean RPMI baseline
#   mean LPS baseline
#   mean BG-rescue baseline
# ============================================================

rpmi_baseline_mean <- rowMeans(
  expr_f[
    ,
    c(
      "RPMI_d6_8746",
      "RPMI_d6_9243",
      "RPMI_d6_10647"
    ),
    drop = FALSE
  ],
  na.rm = TRUE
)


lps_baseline_mean <- rowMeans(
  expr_f[
    ,
    c(
      "LPS_d6_8747",
      "LPS_d6_9244",
      "LPS_d6_10649"
    ),
    drop = FALSE
  ],
  na.rm = TRUE
)


rescue_baseline_mean <- rowMeans(
  expr_f[
    ,
    c(
      "RESCUE_d6_10645",
      "RESCUE_d6_10653"
    ),
    drop = FALSE
  ],
  na.rm = TRUE
)


baseline_memory_effect_joint <-
  lps_baseline_mean -
  rpmi_baseline_mean


BG_baseline_shift <-
  rescue_baseline_mean -
  lps_baseline_mean


# ============================================================
# 19. HELPER FOR SAFE RATIOS
# ============================================================

safe_ratio <- function(
    numerator,
    denominator,
    epsilon = 1e-8
) {

  out <- rep(
    NA_real_,
    length(numerator)
  )

  # Preserve ENSEMBL gene names
  names(out) <- names(numerator)

  ok <-
    is.finite(numerator) &
    is.finite(denominator) &
    abs(denominator) > epsilon

  out[ok] <-
    numerator[ok] /
    denominator[ok]

  out
}


# ============================================================
# 20. BASELINE RESCUE FRACTION
#
# 0 = no movement from LPS state
# 1 = reaches RPMI mean
# >1 = overshoots RPMI
# <0 = moves away from RPMI
# ============================================================

BG_baseline_rescue_fraction <- safe_ratio(
  
  lps_baseline_mean -
    rescue_baseline_mean,
  
  lps_baseline_mean -
    rpmi_baseline_mean
)


BG_baseline_closer_to_RPMI <-
  
  abs(
    rescue_baseline_mean -
      rpmi_baseline_mean
  ) <
  
  abs(
    lps_baseline_mean -
      rpmi_baseline_mean
  )


BG_baseline_opposes_memory <-
  
  sign(BG_baseline_shift) ==
  -sign(baseline_memory_effect_joint)


# ============================================================
# 21. RESPONSE RESTORATION METRICS
#
# Compare:
#
# naive response
# tolerant response
# BG-rescue response
#
#
# Restoration fraction:
#
# (rescue - tolerant)
# ---------------------
# (naive  - tolerant)
#
# 0 = no restoration
# 1 = reaches naive/RPMI response
# ============================================================

response_restoration_fraction <- safe_ratio(
  
  rescue_response_mean -
    tolerant_response_mean,
  
  naive_response_mean -
    tolerant_response_mean
)


BG_response_moves_toward_naive <-
  
  (
    rescue_response_mean -
      tolerant_response_mean
  ) *
  
  (
    naive_response_mean -
      tolerant_response_mean
  ) > 0


BG_response_closer_to_naive <-
  
  abs(
    rescue_response_mean -
      naive_response_mean
  ) <
  
  abs(
    tolerant_response_mean -
      naive_response_mean
  )


# ============================================================
# 22. DESCRIPTIVE RESPONSE PATTERN
#
# This avoids the incorrect assumption that a negative
# interaction always means tolerance.
#
# For a positively induced gene, blunting generally produces
# a negative interaction.
#
# For a negatively responsive gene, blunting can produce
# a positive interaction.
#
# Therefore classify using response magnitude/direction.
# ============================================================

classify_response <- function(
    naive,
    tolerant
) {

  out <- rep(
    NA_character_,
    length(naive)
  )

  # Preserve ENSEMBL gene names
  names(out) <- names(naive)

  ok <-
    is.finite(naive) &
    is.finite(tolerant)

  same_direction <-
    sign(naive) ==
    sign(tolerant)

  out[
    ok &
    !same_direction
  ] <- "direction_changed"

  out[
    ok &
    same_direction &
    abs(tolerant) <
      abs(naive)
  ] <- "blunted_same_direction"

  out[
    ok &
    same_direction &
    abs(tolerant) >=
      abs(naive)
  ] <- "maintained_or_enhanced"

  out
}


response_pattern <- classify_response(
  naive_response_mean,
  tolerant_response_mean
)


# ============================================================
# 23. BUILD 24-MEMORY-KINASE INTEGRATED STEP-04 TABLE
# ============================================================

candidate <- memory24


# Rename Step-03 columns explicitly so there is no confusion
names(candidate)[
  names(candidate) ==
    "log2FC_Day6_LPS_vs_RPMI"
] <- "Memory_Day6_log2FC"


names(candidate)[
  names(candidate) ==
    "P.Value"
] <- "Memory_Day6_P"


names(candidate)[
  names(candidate) ==
    "FDR_kinase"
] <- "Memory_Day6_FDR_kinase"


names(candidate)[
  names(candidate) ==
    "FDR_genome"
] <- "Memory_Day6_FDR_genome"


# ------------------------------------------------------------
# Match expression-level metrics
# ------------------------------------------------------------

candidate$RPMI_baseline_mean_log2 <-
  rpmi_baseline_mean[
    candidate$ENSEMBL
  ]


candidate$LPS_baseline_mean_log2 <-
  lps_baseline_mean[
    candidate$ENSEMBL
  ]


candidate$BG_rescue_baseline_mean_log2 <-
  rescue_baseline_mean[
    candidate$ENSEMBL
  ]


candidate$Baseline_memory_effect_joint <-
  baseline_memory_effect_joint[
    candidate$ENSEMBL
  ]


candidate$BG_baseline_shift_from_LPS <-
  BG_baseline_shift[
    candidate$ENSEMBL
  ]


candidate$BG_baseline_rescue_fraction <-
  BG_baseline_rescue_fraction[
    candidate$ENSEMBL
  ]


candidate$BG_baseline_closer_to_RPMI <-
  BG_baseline_closer_to_RPMI[
    candidate$ENSEMBL
  ]


candidate$BG_baseline_opposes_memory <-
  BG_baseline_opposes_memory[
    candidate$ENSEMBL
  ]


# ------------------------------------------------------------
# Restimulation response metrics
# ------------------------------------------------------------

candidate$Naive_response_log2 <-
  naive_response_mean[
    candidate$ENSEMBL
  ]


candidate$Tolerant_response_log2 <-
  tolerant_response_mean[
    candidate$ENSEMBL
  ]


candidate$BG_rescue_response_log2 <-
  rescue_response_mean[
    candidate$ENSEMBL
  ]


candidate$Tolerance_interaction_log2 <-
  manual_interaction[
    candidate$ENSEMBL
  ]


candidate$Response_pattern <-
  response_pattern[
    candidate$ENSEMBL
  ]


candidate$BG_response_restoration_fraction <-
  response_restoration_fraction[
    candidate$ENSEMBL
  ]


candidate$BG_response_moves_toward_naive <-
  BG_response_moves_toward_naive[
    candidate$ENSEMBL
  ]


candidate$BG_response_closer_to_naive <-
  BG_response_closer_to_naive[
    candidate$ENSEMBL
  ]


# ============================================================
# 24. ADD FORMAL INTERACTION STATISTICS
# ============================================================

formal_idx <- match(
  candidate$ENSEMBL,
  formal_results$ENSEMBL
)


candidate$Interaction_P <-
  formal_results$Interaction_P[
    formal_idx
  ]


candidate$Interaction_FDR_genome <-
  formal_results$Interaction_FDR_genome[
    formal_idx
  ]


# Candidate-family FDR:
#
# Secondary/exploratory because these 24 genes were selected
# previously from the Day-6 memory analysis.
candidate$Interaction_FDR_memory24_exploratory <-
  p.adjust(
    candidate$Interaction_P,
    method = "BH"
  )


candidate$Interaction_genomeFDR05 <-
  candidate$Interaction_FDR_genome <
  0.05


candidate$Interaction_memory24_FDR05_exploratory <-
  candidate$Interaction_FDR_memory24_exploratory <
  0.05


# ============================================================
# 25. CLASSICAL BLUNTING FLAGS
#
# Formal interaction significance + descriptive response
# magnitude pattern.
# ============================================================

candidate$Formal_blunting_genomeFDR05 <-
  
  candidate$Interaction_genomeFDR05 &
  candidate$Response_pattern ==
  "blunted_same_direction"


candidate$Exploratory_blunting_memory24_FDR05 <-
  
  candidate$Interaction_memory24_FDR05_exploratory &
  candidate$Response_pattern ==
  "blunted_same_direction"


# ============================================================
# 26. SAVE BASELINE BG RESCUE TABLE
# ============================================================

baseline_table <- candidate[
  ,
  c(
    "ENSEMBL",
    "SYMBOL",
    "Memory_Day6_log2FC",
    "Memory_Day6_FDR_kinase",
    "RPMI_baseline_mean_log2",
    "LPS_baseline_mean_log2",
    "BG_rescue_baseline_mean_log2",
    "Baseline_memory_effect_joint",
    "BG_baseline_shift_from_LPS",
    "BG_baseline_rescue_fraction",
    "BG_baseline_closer_to_RPMI",
    "BG_baseline_opposes_memory"
  )
]


write.csv(
  baseline_table,
  file.path(
    outdir,
    "04_memory24_baseline_BG_rescue.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 27. SAVE RESTIMULATION TABLE
# ============================================================

response_table <- candidate[
  ,
  c(
    "ENSEMBL",
    "SYMBOL",
    "Memory_Day6_log2FC",
    "Naive_response_log2",
    "Tolerant_response_log2",
    "Tolerance_interaction_log2",
    "Interaction_P",
    "Interaction_FDR_genome",
    "Interaction_FDR_memory24_exploratory",
    "Response_pattern",
    "BG_rescue_response_log2",
    "BG_response_restoration_fraction",
    "BG_response_moves_toward_naive",
    "BG_response_closer_to_naive",
    "Formal_blunting_genomeFDR05",
    "Exploratory_blunting_memory24_FDR05"
  )
]


write.csv(
  response_table,
  file.path(
    outdir,
    "04_memory24_restimulation_response.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 28. SAVE COMPLETE INTEGRATED STEP-04 TABLE
# ============================================================

candidate <- candidate[
  order(
    candidate$Interaction_P
  ),
  ,
  drop = FALSE
]


write.csv(
  candidate,
  file.path(
    outdir,
    "04_memory24_BG_rescue_restimulation_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 29. FIGURE 04A
#
# Naive vs tolerant response
# ============================================================

png(
  file.path(
    figdir,
    "04A_memory24_naive_vs_tolerant_response.png"
  ),
  width = 1800,
  height = 1500,
  res = 180
)


plot(
  candidate$Naive_response_log2,
  candidate$Tolerant_response_log2,
  pch = 16,
  cex = 0.9,
  xlab =
    "RPMI / naive restimulation response (log2)",
  ylab =
    "LPS-conditioned restimulation response (log2)",
  main =
    "Memory kinases: naive vs LPS-conditioned response"
)


abline(
  a = 0,
  b = 1,
  lty = 2
)


abline(
  h = 0,
  lty = 3
)


abline(
  v = 0,
  lty = 3
)


label_n <- min(
  12,
  nrow(candidate)
)


label_idx <- seq_len(
  label_n
)


text(
  candidate$Naive_response_log2[
    label_idx
  ],
  candidate$Tolerant_response_log2[
    label_idx
  ],
  labels =
    candidate$SYMBOL[
      label_idx
    ],
  pos = 3,
  cex = 0.65
)


dev.off()


# ============================================================
# 30. FIGURE 04B
#
# Tolerant vs BG-rescue response
# ============================================================

movement <- abs(
  candidate$BG_rescue_response_log2 -
    candidate$Tolerant_response_log2
)


movement_order <- order(
  movement,
  decreasing = TRUE,
  na.last = NA
)


label_idx_2 <- head(
  movement_order,
  12
)


png(
  file.path(
    figdir,
    "04B_memory24_tolerant_vs_BG_rescue_response.png"
  ),
  width = 1800,
  height = 1500,
  res = 180
)


plot(
  candidate$Tolerant_response_log2,
  candidate$BG_rescue_response_log2,
  pch = 16,
  cex = 0.9,
  xlab =
    "LPS-conditioned restimulation response (log2)",
  ylab =
    "BG-rescue restimulation response (log2)",
  main =
    "Memory kinases: tolerant vs beta-glucan rescue response"
)


abline(
  a = 0,
  b = 1,
  lty = 2
)


abline(
  h = 0,
  lty = 3
)


abline(
  v = 0,
  lty = 3
)


text(
  candidate$Tolerant_response_log2[
    label_idx_2
  ],
  candidate$BG_rescue_response_log2[
    label_idx_2
  ],
  labels =
    candidate$SYMBOL[
      label_idx_2
    ],
  pos = 3,
  cex = 0.65
)


dev.off()


# ============================================================
# 31. SUMMARY COUNTS
# ============================================================

n_interaction_genome_fdr05 <- sum(
  candidate$Interaction_genomeFDR05 %in%
    TRUE
)


n_interaction_candidate_fdr05 <- sum(
  candidate$Interaction_memory24_FDR05_exploratory %in%
    TRUE
)


n_blunted <- sum(
  candidate$Response_pattern ==
    "blunted_same_direction",
  na.rm = TRUE
)


n_baseline_closer <- sum(
  candidate$BG_baseline_closer_to_RPMI %in%
    TRUE
)


n_response_closer <- sum(
  candidate$BG_response_closer_to_naive %in%
    TRUE
)


n_response_toward <- sum(
  candidate$BG_response_moves_toward_naive %in%
    TRUE
)


n_formal_blunting <- sum(
  candidate$Formal_blunting_genomeFDR05 %in%
    TRUE
)


n_exploratory_blunting <- sum(
  candidate$Exploratory_blunting_memory24_FDR05 %in%
    TRUE
)


# ============================================================
# 32. PRINT TOP STEP-04 RESULTS
# ============================================================

display_cols <- c(
  "SYMBOL",
  "Memory_Day6_log2FC",
  "Naive_response_log2",
  "Tolerant_response_log2",
  "Tolerance_interaction_log2",
  "Interaction_P",
  "Interaction_FDR_genome",
  "Interaction_FDR_memory24_exploratory",
  "Response_pattern",
  "BG_rescue_response_log2",
  "BG_response_restoration_fraction",
  "BG_response_closer_to_naive"
)


cat("\n========================================\n")
cat("MEMORY-KINASE RESTIMULATION RESULTS\n")
cat("========================================\n\n")


print(
  candidate[
    ,
    display_cols
  ],
  row.names = FALSE,
  digits = 4
)


# ============================================================
# 33. SAVE RDS
# ============================================================

saveRDS(
  list(
    
    analysis_manifest =
      analysis_info,
    
    expression_log2 =
      expr_f,
    
    formal_manifest =
      formal_info,
    
    formal_design =
      design,
    
    formal_results =
      formal_results,
    
    response_matrix =
      response_matrix,
    
    memory24_summary =
      candidate
    
  ),
  file.path(
    outdir,
    "04_BG_rescue_and_restimulation.rds"
  )
)


# ============================================================
# 34. TEXT SUMMARY
# ============================================================

summary_lines <- c(
  
  "GSE85246 — Step 04 BG Rescue + Restimulation",
  
  "",
  
  paste(
    "Primary Step-03 memory kinases:",
    nrow(candidate)
  ),
  
  paste(
    "Jointly normalized Step-04 samples:",
    ncol(expr_f)
  ),
  
  paste(
    "Genes retained:",
    nrow(expr_f)
  ),
  
  "",
  
  "Formal tolerance interaction:",
  "3 matched donors x RPMI/LPS history x baseline/restimulation",
  
  "",
  
  paste(
    "Memory kinases with genome-wide interaction FDR < 0.05:",
    n_interaction_genome_fdr05
  ),
  
  paste(
    "Memory kinases with exploratory candidate-family interaction FDR < 0.05:",
    n_interaction_candidate_fdr05
  ),
  
  paste(
    "Memory kinases with descriptive blunted response:",
    n_blunted
  ),
  
  paste(
    "Formal blunting + genome-wide interaction FDR < 0.05:",
    n_formal_blunting
  ),
  
  paste(
    "Exploratory blunting + memory24 interaction FDR < 0.05:",
    n_exploratory_blunting
  ),
  
  "",
  
  paste(
    "Memory kinases with BG baseline closer to RPMI:",
    n_baseline_closer
  ),
  
  paste(
    "Memory kinases with BG response moving toward naive response:",
    n_response_toward
  ),
  
  paste(
    "Memory kinases with BG response closer to naive response:",
    n_response_closer
  ),
  
  "",
  
  "Important interpretation:",
  paste(
    "Beta-glucan rescue has only two matched rescue-response pairs;",
    "therefore rescue evidence is treated primarily as directional/descriptive."
  ),
  
  paste(
    "The unmatched RPMI_Restim_8754 sample is not used",
    "in the formal paired interaction."
  )
)


writeLines(
  summary_lines,
  file.path(
    outdir,
    "04_BG_rescue_restimulation_summary.txt"
  )
)


# ============================================================
# 35. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n========================================\n")
cat("04 BG RESCUE + RESTIMULATION COMPLETE\n")
cat("========================================\n\n")


cat(
  "Step-03 memory kinases evaluated:",
  nrow(candidate),
  "\n"
)


cat(
  "Formal matched donors:",
  3,
  "\n"
)


cat(
  "BG rescue response pairs:",
  2,
  "\n"
)


cat(
  "Genome-wide interaction FDR < 0.05:",
  n_interaction_genome_fdr05,
  "\n"
)


cat(
  "Exploratory memory24 interaction FDR < 0.05:",
  n_interaction_candidate_fdr05,
  "\n"
)


cat(
  "Descriptively blunted memory-kinase responses:",
  n_blunted,
  "\n"
)


cat(
  "Formal blunting + genome-wide FDR < 0.05:",
  n_formal_blunting,
  "\n"
)


cat(
  "BG baseline closer to RPMI:",
  n_baseline_closer,
  "/",
  nrow(candidate),
  "\n"
)


cat(
  "BG response moving toward naive response:",
  n_response_toward,
  "/",
  nrow(candidate),
  "\n"
)


cat(
  "BG response closer to naive response:",
  n_response_closer,
  "/",
  nrow(candidate),
  "\n"
)


cat(
  "\nIMPORTANT:\n",
  "BG rescue evidence is directional/descriptive because n = 2.\n",
  sep = ""
)


cat(
  "\nNext canonical step:\n"
)


cat(
  "scripts/05_integrated_kinase_trajectory.R\n"
)


cat("\n========================================\n")