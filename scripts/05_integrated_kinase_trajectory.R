# ============================================================
# 05_integrated_kinase_trajectory.R
#
# GSE85246 / GSE85243
#
# INTEGRATED RNA KINASE TRAJECTORY
#
# Integrates:
#
#   Step 02 — Initial Day-1 LPS response
#   Step 03 — Persistent Day-6 memory state
#   Step 04 — Restimulation/tolerance + beta-glucan rescue
#
#
# IMPORTANT:
#
# - No new differential-expression model is fitted here.
#
# - The primary candidate universe is the 24 canonical
#   KinHub/OpenKinome Day-6 memory kinases selected in Step 03.
#
# - Day-1 evidence is contextual because no Day-1 kinase
#   survived multiple-testing correction.
#
# - Restimulation interaction statistics are retained, but
#   none of the 24 candidates passed interaction FDR < 0.05.
#
# - Beta-glucan rescue evidence is descriptive/directional
#   because only two rescue pairs are available.
#
#
# Main output:
#
#   One integrated trajectory table for all 24 memory kinases.
#
#
# Next:
#   06_make_final_kinase_figure.R
# ============================================================


rm(list = ls())


# ============================================================
# 1. INPUT FILES
# ============================================================

day1_file <-
  "results/02_initial_LPS/02_Day1_LPS_vs_RPMI_KinHub_kinases.csv"


memory_file <-
  "results/03_memory_kinases/03_memory_kinases_FDR05.csv"


step04_file <-
  paste0(
    "results/04_BG_rescue_restimulation/",
    "04_memory24_BG_rescue_restimulation_summary.csv"
  )


required_files <- c(
  day1_file,
  memory_file,
  step04_file
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {
  
  stop(
    paste(
      "Missing required file(s):",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


# ============================================================
# 2. OUTPUT DIRECTORY
# ============================================================

outdir <-
  "results/05_integrated_trajectory"


dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. LOAD RESULTS
# ============================================================

day1 <- read.csv(
  day1_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


memory <- read.csv(
  memory_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


step04 <- read.csv(
  step04_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("\n========================================\n")
cat("05 — INTEGRATED KINASE TRAJECTORY\n")
cat("========================================\n\n")


cat(
  "Day-1 canonical kinases:",
  nrow(day1),
  "\n"
)


cat(
  "Step-03 primary memory kinases:",
  nrow(memory),
  "\n"
)


cat(
  "Step-04 memory kinases:",
  nrow(step04),
  "\n"
)


# ============================================================
# 4. BASIC VALIDATION
# ============================================================

if (anyDuplicated(day1$ENSEMBL)) {
  
  stop(
    "Duplicate ENSEMBL IDs in Day-1 kinase table."
  )
}


if (anyDuplicated(memory$ENSEMBL)) {
  
  stop(
    "Duplicate ENSEMBL IDs in Step-03 memory table."
  )
}


if (anyDuplicated(step04$ENSEMBL)) {
  
  stop(
    "Duplicate ENSEMBL IDs in Step-04 table."
  )
}


if (nrow(memory) != 24) {
  
  warning(
    paste(
      "Expected 24 Step-03 memory kinases but found",
      nrow(memory)
    )
  )
}


# ============================================================
# 5. VERIFY STEP-03 AND STEP-04 CANDIDATE SETS
# ============================================================

same_candidate_set <-
  setequal(
    memory$ENSEMBL,
    step04$ENSEMBL
  )


cat(
  "\nStep-03 and Step-04 candidate sets identical:",
  same_candidate_set,
  "\n"
)


if (!same_candidate_set) {
  
  cat(
    "\nPresent in Step 03 but absent from Step 04:\n"
  )
  
  print(
    setdiff(
      memory$ENSEMBL,
      step04$ENSEMBL
    )
  )
  
  
  cat(
    "\nPresent in Step 04 but absent from Step 03:\n"
  )
  
  print(
    setdiff(
      step04$ENSEMBL,
      memory$ENSEMBL
    )
  )
  
  
  stop(
    "Step-03 and Step-04 candidate sets do not match."
  )
}


# ============================================================
# 6. START FROM STEP-03 MEMORY CANDIDATES
# ============================================================

integrated <- memory


# Explicitly rename memory columns
names(integrated)[
  names(integrated) ==
    "log2FC_Day6_LPS_vs_RPMI"
] <- "Day6_memory_log2FC"


names(integrated)[
  names(integrated) ==
    "P.Value"
] <- "Day6_memory_P"


names(integrated)[
  names(integrated) ==
    "FDR_kinase"
] <- "Day6_memory_FDR_kinase"


names(integrated)[
  names(integrated) ==
    "FDR_genome"
] <- "Day6_memory_FDR_genome"


names(integrated)[
  names(integrated) ==
    "Direction"
] <- "Day6_memory_direction"


# ============================================================
# 7. ADD DAY-1 INITIAL RESPONSE
# ============================================================

d1_idx <- match(
  integrated$ENSEMBL,
  day1$ENSEMBL
)


if (anyNA(d1_idx)) {
  
  cat(
    "\nMemory candidates missing from Day-1 kinase table:\n"
  )
  
  print(
    integrated[
      is.na(d1_idx),
      c(
        "ENSEMBL",
        "SYMBOL"
      )
    ],
    row.names = FALSE
  )
  
  stop(
    "One or more memory kinases cannot be evaluated at Day 1."
  )
}


integrated$Day1_log2FC <-
  day1$log2FC_Day1_LPS_vs_RPMI[
    d1_idx
  ]


integrated$Day1_P <-
  day1$P.Value[
    d1_idx
  ]


integrated$Day1_FDR_kinase <-
  day1$FDR_kinase[
    d1_idx
  ]


integrated$Day1_FDR_genome <-
  day1$FDR_genome[
    d1_idx
  ]


integrated$Day1_direction <- ifelse(
  integrated$Day1_log2FC > 0,
  "UP",
  ifelse(
    integrated$Day1_log2FC < 0,
    "DOWN",
    "UNCHANGED"
  )
)


integrated$Day1_nominal_P05 <-
  integrated$Day1_P < 0.05


integrated$Day1_kinase_FDR05 <-
  integrated$Day1_FDR_kinase < 0.05


# ============================================================
# 8. INITIAL-RESPONSE / MEMORY RELATIONSHIP
#
# No Day-1 kinase survived FDR correction in Step 02.
#
# Therefore:
#
# acute_supported_same_direction
#   Day1 nominal P < 0.05 and Day1/Day6 directions agree.
#
# acute_nominal_opposite_direction
#   Day1 nominal P < 0.05 but direction differs.
#
# memory_emergent_small_Day1_effect
#   No nominal Day1 evidence AND |Day1 log2FC| < 0.25.
#
# same_direction_no_nominal_support
#   Direction agrees but Day1 nominal P >= 0.05.
#
# opposite_direction_no_nominal_support
#   Direction differs and Day1 nominal P >= 0.05.
# ============================================================

same_direction <-
  sign(
    integrated$Day1_log2FC
  ) ==
  sign(
    integrated$Day6_memory_log2FC
  )


small_day1_effect <-
  abs(
    integrated$Day1_log2FC
  ) < 0.25


integrated$Day1_Day6_relationship <-
  NA_character_


integrated$Day1_Day6_relationship[
  integrated$Day1_nominal_P05 &
    same_direction
] <- "acute_supported_same_direction"


integrated$Day1_Day6_relationship[
  integrated$Day1_nominal_P05 &
    !same_direction
] <- "acute_nominal_opposite_direction"


integrated$Day1_Day6_relationship[
  !integrated$Day1_nominal_P05 &
    small_day1_effect
] <- "memory_emergent_small_Day1_effect"


integrated$Day1_Day6_relationship[
  !integrated$Day1_nominal_P05 &
    !small_day1_effect &
    same_direction
] <- "same_direction_no_nominal_support"


integrated$Day1_Day6_relationship[
  !integrated$Day1_nominal_P05 &
    !small_day1_effect &
    !same_direction
] <- "opposite_direction_no_nominal_support"


# ============================================================
# 9. ADD STEP-04 RESTIMULATION + BG RESULTS
# ============================================================

s4_idx <- match(
  integrated$ENSEMBL,
  step04$ENSEMBL
)


if (anyNA(s4_idx)) {
  
  stop(
    "Failed to match one or more memory kinases to Step 04."
  )
}


copy_step04_column <- function(column_name) {
  
  if (!column_name %in%
      names(step04)) {
    
    stop(
      paste(
        "Required Step-04 column missing:",
        column_name
      )
    )
  }
  
  step04[
    s4_idx,
    column_name
  ]
}


integrated$Naive_response_log2 <-
  copy_step04_column(
    "Naive_response_log2"
  )


integrated$Tolerant_response_log2 <-
  copy_step04_column(
    "Tolerant_response_log2"
  )


integrated$Tolerance_interaction_log2 <-
  copy_step04_column(
    "Tolerance_interaction_log2"
  )


integrated$Interaction_P <-
  copy_step04_column(
    "Interaction_P"
  )


integrated$Interaction_FDR_genome <-
  copy_step04_column(
    "Interaction_FDR_genome"
  )


integrated$Interaction_FDR_memory24_exploratory <-
  copy_step04_column(
    "Interaction_FDR_memory24_exploratory"
  )


integrated$Response_pattern <-
  copy_step04_column(
    "Response_pattern"
  )


integrated$BG_rescue_response_log2 <-
  copy_step04_column(
    "BG_rescue_response_log2"
  )


integrated$BG_response_restoration_fraction <-
  copy_step04_column(
    "BG_response_restoration_fraction"
  )


integrated$BG_response_moves_toward_naive <-
  copy_step04_column(
    "BG_response_moves_toward_naive"
  )


integrated$BG_response_closer_to_naive <-
  copy_step04_column(
    "BG_response_closer_to_naive"
  )


integrated$RPMI_baseline_mean_log2 <-
  copy_step04_column(
    "RPMI_baseline_mean_log2"
  )


integrated$LPS_baseline_mean_log2 <-
  copy_step04_column(
    "LPS_baseline_mean_log2"
  )


integrated$BG_rescue_baseline_mean_log2 <-
  copy_step04_column(
    "BG_rescue_baseline_mean_log2"
  )


integrated$BG_baseline_rescue_fraction <-
  copy_step04_column(
    "BG_baseline_rescue_fraction"
  )


integrated$BG_baseline_closer_to_RPMI <-
  copy_step04_column(
    "BG_baseline_closer_to_RPMI"
  )


integrated$BG_baseline_opposes_memory <-
  copy_step04_column(
    "BG_baseline_opposes_memory"
  )


# ============================================================
# 10. ROBUST LOGICAL CONVERSION
# ============================================================

as_logical_safe <- function(x) {
  
  if (is.logical(x)) {
    
    return(x)
  }
  
  
  out <- rep(
    NA,
    length(x)
  )
  
  
  x_char <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  
  out[
    x_char == "TRUE"
  ] <- TRUE
  
  
  out[
    x_char == "FALSE"
  ] <- FALSE
  
  
  out
}


integrated$BG_response_moves_toward_naive <-
  as_logical_safe(
    integrated$BG_response_moves_toward_naive
  )


integrated$BG_response_closer_to_naive <-
  as_logical_safe(
    integrated$BG_response_closer_to_naive
  )


integrated$BG_baseline_closer_to_RPMI <-
  as_logical_safe(
    integrated$BG_baseline_closer_to_RPMI
  )


integrated$BG_baseline_opposes_memory <-
  as_logical_safe(
    integrated$BG_baseline_opposes_memory
  )


# ============================================================
# 11. FORMAL / EXPLORATORY INTERACTION FLAGS
# ============================================================

integrated$Interaction_nominal_P05 <-
  integrated$Interaction_P < 0.05


integrated$Interaction_genome_FDR05 <-
  integrated$Interaction_FDR_genome <
  0.05


integrated$Interaction_memory24_FDR05_exploratory <-
  integrated$
  Interaction_FDR_memory24_exploratory <
  0.05


# ============================================================
# 12. RESPONSE-PATTERN FLAGS
# ============================================================

integrated$Blunted_same_direction <-
  integrated$Response_pattern ==
  "blunted_same_direction"


integrated$Direction_changed_after_restim <-
  integrated$Response_pattern ==
  "direction_changed"


integrated$Maintained_or_enhanced_response <-
  integrated$Response_pattern ==
  "maintained_or_enhanced"


# ============================================================
# 13. BG BASELINE REVERSAL SUPPORT
#
# Require both:
#
# 1. BG baseline is closer to RPMI than LPS baseline.
# 2. BG shift is opposite to the Day-6 memory direction.
#
# This remains directional evidence, not a formal paired test.
# ============================================================

integrated$BG_baseline_reversal_supported <-
  
  integrated$BG_baseline_closer_to_RPMI %in%
  TRUE &
  
  integrated$BG_baseline_opposes_memory %in%
  TRUE


# ============================================================
# 14. BG RESPONSE RESTORATION SUPPORT
#
# Require both:
#
# 1. BG response moves from tolerant toward naive.
# 2. BG response finishes closer to naive than tolerant.
#
# Again, descriptive because rescue n = 2.
# ============================================================

integrated$BG_response_restoration_supported <-
  
  integrated$BG_response_moves_toward_naive %in%
  TRUE &
  
  integrated$BG_response_closer_to_naive %in%
  TRUE


# ============================================================
# 15. DIRECTIONAL RNA SUPPORT SCORE
#
# All 24 already satisfy the primary Day-6 memory criterion.
#
# Additional independent directional layers:
#
# +1  BG reverses baseline memory state
# +1  Restimulation response is blunted in same direction
# +1  BG restores response toward naive state
#
# Score ranges from 0 to 3.
#
# NOTE:
# Nominal interaction P values are NOT included in the score
# because they do not survive multiple-testing correction.
# ============================================================

integrated$Directional_support_score <-
  
  as.integer(
    integrated$BG_baseline_reversal_supported
  ) +
  
  as.integer(
    integrated$Blunted_same_direction
  ) +
  
  as.integer(
    integrated$BG_response_restoration_supported
  )


# ============================================================
# 16. RNA PRIORITY FOR CHROMATIN FOLLOW-UP
#
# HIGH:
#   all three directional layers supported
#
# INTERMEDIATE:
#   two of the three supported
#
# EXPLORATORY:
#   zero or one supported
#
# This is a prioritisation category, NOT a statistical
# significance category.
# ============================================================

integrated$RNA_chromatin_priority <- ifelse(
  
  integrated$Directional_support_score == 3,
  
  "HIGH",
  
  ifelse(
    
    integrated$Directional_support_score == 2,
    
    "INTERMEDIATE",
    
    "EXPLORATORY"
  )
)


# ============================================================
# 17. HUMAN-READABLE INTEGRATED INTERPRETATION
# ============================================================

integrated$Integrated_RNA_interpretation <-
  paste0(
    "Day1=",
    integrated$Day1_Day6_relationship,
    "; Day6=memory_FDR05",
    "; Restim=",
    integrated$Response_pattern,
    "; BG_baseline_reversal=",
    integrated$BG_baseline_reversal_supported,
    "; BG_response_restoration=",
    integrated$BG_response_restoration_supported
  )


# ============================================================
# 18. SANITY CHECKS
# ============================================================

cat("\n========================================\n")
cat("INTEGRATION VALIDATION\n")
cat("========================================\n\n")


cat(
  "Integrated kinase rows:",
  nrow(integrated),
  "\n"
)


cat(
  "Missing Day-1 log2FC:",
  sum(
    is.na(
      integrated$Day1_log2FC
    )
  ),
  "\n"
)


cat(
  "Missing Response_pattern:",
  sum(
    is.na(
      integrated$Response_pattern
    )
  ),
  "\n"
)


cat(
  "Missing BG response-restoration flag:",
  sum(
    is.na(
      integrated$BG_response_restoration_supported
    )
  ),
  "\n"
)


if (anyNA(
  integrated$Day1_Day6_relationship
)) {
  
  warning(
    "Some Day1/Day6 trajectory classes are missing."
  )
}


# ============================================================
# 19. SORT TABLE
# ============================================================

priority_rank <- c(
  "HIGH" = 1,
  "INTERMEDIATE" = 2,
  "EXPLORATORY" = 3
)


integrated$priority_rank <-
  unname(
    priority_rank[
      integrated$RNA_chromatin_priority
    ]
  )


integrated <- integrated[
  order(
    integrated$priority_rank,
    -integrated$Directional_support_score,
    integrated$Interaction_P,
    integrated$Day6_memory_FDR_kinase
  ),
  ,
  drop = FALSE
]


integrated$priority_rank <- NULL


# ============================================================
# 20. SAVE COMPLETE INTEGRATED TABLE
# ============================================================

write.csv(
  integrated,
  file.path(
    outdir,
    "05_integrated_kinase_trajectory.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. SAVE PRIORITY SUBSETS
# ============================================================

high_priority <- integrated[
  integrated$RNA_chromatin_priority ==
    "HIGH",
  ,
  drop = FALSE
]


intermediate_priority <- integrated[
  integrated$RNA_chromatin_priority ==
    "INTERMEDIATE",
  ,
  drop = FALSE
]


exploratory_priority <- integrated[
  integrated$RNA_chromatin_priority ==
    "EXPLORATORY",
  ,
  drop = FALSE
]


write.csv(
  high_priority,
  file.path(
    outdir,
    "05_HIGH_priority_for_chromatin.csv"
  ),
  row.names = FALSE
)


write.csv(
  intermediate_priority,
  file.path(
    outdir,
    "05_INTERMEDIATE_priority_for_chromatin.csv"
  ),
  row.names = FALSE
)


write.csv(
  exploratory_priority,
  file.path(
    outdir,
    "05_EXPLORATORY_priority_for_chromatin.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. SAVE COMPACT REVIEW TABLE
# ============================================================

review_cols <- c(
  "ENSEMBL",
  "SYMBOL",
  "Day1_log2FC",
  "Day1_P",
  "Day1_FDR_kinase",
  "Day1_Day6_relationship",
  "Day6_memory_log2FC",
  "Day6_memory_FDR_kinase",
  "Day6_memory_FDR_genome",
  "Response_pattern",
  "Tolerance_interaction_log2",
  "Interaction_P",
  "Interaction_FDR_genome",
  "BG_baseline_reversal_supported",
  "BG_response_restoration_supported",
  "Directional_support_score",
  "RNA_chromatin_priority"
)


review_table <- integrated[
  ,
  review_cols,
  drop = FALSE
]


write.csv(
  review_table,
  file.path(
    outdir,
    "05_integrated_kinase_review_table.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. SUMMARY COUNTS
# ============================================================

n_acute_same <- sum(
  integrated$Day1_Day6_relationship ==
    "acute_supported_same_direction",
  na.rm = TRUE
)


n_acute_opposite <- sum(
  integrated$Day1_Day6_relationship ==
    "acute_nominal_opposite_direction",
  na.rm = TRUE
)


n_memory_emergent <- sum(
  integrated$Day1_Day6_relationship ==
    "memory_emergent_small_Day1_effect",
  na.rm = TRUE
)


n_same_non_nominal <- sum(
  integrated$Day1_Day6_relationship ==
    "same_direction_no_nominal_support",
  na.rm = TRUE
)


n_opposite_non_nominal <- sum(
  integrated$Day1_Day6_relationship ==
    "opposite_direction_no_nominal_support",
  na.rm = TRUE
)


n_nominal_interaction <- sum(
  integrated$Interaction_nominal_P05,
  na.rm = TRUE
)


n_interaction_fdr <- sum(
  integrated$Interaction_genome_FDR05,
  na.rm = TRUE
)


n_bg_baseline <- sum(
  integrated$BG_baseline_reversal_supported,
  na.rm = TRUE
)


n_bg_response <- sum(
  integrated$BG_response_restoration_supported,
  na.rm = TRUE
)


n_blunted <- sum(
  integrated$Blunted_same_direction,
  na.rm = TRUE
)


# ============================================================
# 24. CONSOLE TABLE
# ============================================================

console_cols <- c(
  "SYMBOL",
  "Day1_log2FC",
  "Day1_P",
  "Day1_Day6_relationship",
  "Day6_memory_log2FC",
  "Day6_memory_FDR_kinase",
  "Response_pattern",
  "Interaction_P",
  "BG_baseline_reversal_supported",
  "BG_response_restoration_supported",
  "Directional_support_score",
  "RNA_chromatin_priority"
)


cat("\n========================================\n")
cat("INTEGRATED MEMORY-KINASE TRAJECTORIES\n")
cat("========================================\n\n")


print(
  integrated[
    ,
    console_cols,
    drop = FALSE
  ],
  row.names = FALSE,
  digits = 4
)


# ============================================================
# 25. SAVE RDS
# ============================================================

saveRDS(
  list(
    
    integrated =
      integrated,
    
    high_priority =
      high_priority,
    
    intermediate_priority =
      intermediate_priority,
    
    exploratory_priority =
      exploratory_priority
    
  ),
  file.path(
    outdir,
    "05_integrated_kinase_trajectory.rds"
  )
)


# ============================================================
# 26. TEXT SUMMARY
# ============================================================

summary_lines <- c(
  
  "GSE85246 — Step 05 Integrated Kinase Trajectory",
  
  "",
  
  paste(
    "Primary Day-6 memory kinases:",
    nrow(integrated)
  ),
  
  "",
  
  "Day-1 / Day-6 relationship:",
  
  paste(
    "Acute-supported same direction:",
    n_acute_same
  ),
  
  paste(
    "Acute nominal opposite direction:",
    n_acute_opposite
  ),
  
  paste(
    "Memory-emergent small Day-1 effect:",
    n_memory_emergent
  ),
  
  paste(
    "Same direction without nominal Day-1 support:",
    n_same_non_nominal
  ),
  
  paste(
    "Opposite direction without nominal Day-1 support:",
    n_opposite_non_nominal
  ),
  
  "",
  
  paste(
    "Blunted same-direction restimulation response:",
    n_blunted
  ),
  
  paste(
    "Nominal interaction P < 0.05:",
    n_nominal_interaction
  ),
  
  paste(
    "Genome-wide interaction FDR < 0.05:",
    n_interaction_fdr
  ),
  
  "",
  
  paste(
    "BG baseline reversal supported:",
    n_bg_baseline
  ),
  
  paste(
    "BG response restoration supported:",
    n_bg_response
  ),
  
  "",
  
  paste(
    "HIGH RNA/chromatin priority:",
    nrow(high_priority)
  ),
  
  paste(
    "INTERMEDIATE RNA/chromatin priority:",
    nrow(intermediate_priority)
  ),
  
  paste(
    "EXPLORATORY RNA/chromatin priority:",
    nrow(exploratory_priority)
  ),
  
  "",
  
  paste(
    "IMPORTANT:",
    "RNA_chromatin_priority is an evidence-integration",
    "category, not an FDR significance category."
  ),
  
  paste(
    "Beta-glucan rescue evidence is directional/descriptive",
    "because only two rescue pairs are available."
  )
)


writeLines(
  summary_lines,
  file.path(
    outdir,
    "05_integrated_kinase_trajectory_summary.txt"
  )
)


# ============================================================
# 27. FINAL REPORT
# ============================================================

cat("\n========================================\n")
cat("05 INTEGRATED KINASE TRAJECTORY COMPLETE\n")
cat("========================================\n\n")


cat(
  "Primary memory kinases:",
  nrow(integrated),
  "\n"
)


cat("\nDay-1 / Day-6 classification:\n")


cat(
  "  Acute-supported same direction:",
  n_acute_same,
  "\n"
)


cat(
  "  Acute nominal opposite direction:",
  n_acute_opposite,
  "\n"
)


cat(
  "  Memory-emergent small Day-1 effect:",
  n_memory_emergent,
  "\n"
)


cat(
  "  Same direction, no nominal Day-1 support:",
  n_same_non_nominal,
  "\n"
)


cat(
  "  Opposite direction, no nominal Day-1 support:",
  n_opposite_non_nominal,
  "\n"
)


cat("\nRestimulation:\n")


cat(
  "  Blunted same-direction:",
  n_blunted,
  "\n"
)


cat(
  "  Nominal interaction P < 0.05:",
  n_nominal_interaction,
  "\n"
)


cat(
  "  Genome-wide interaction FDR < 0.05:",
  n_interaction_fdr,
  "\n"
)


cat("\nBeta-glucan:\n")


cat(
  "  Baseline reversal supported:",
  n_bg_baseline,
  "/",
  nrow(integrated),
  "\n"
)


cat(
  "  Response restoration supported:",
  n_bg_response,
  "/",
  nrow(integrated),
  "\n"
)


cat("\nRNA priority for chromatin follow-up:\n")


cat(
  "  HIGH:",
  nrow(high_priority),
  "\n"
)


cat(
  "  INTERMEDIATE:",
  nrow(intermediate_priority),
  "\n"
)


cat(
  "  EXPLORATORY:",
  nrow(exploratory_priority),
  "\n"
)


if (nrow(high_priority) > 0) {
  
  cat(
    "\nHIGH-priority kinase symbols:\n"
  )
  
  cat(
    paste(
      high_priority$SYMBOL,
      collapse = ", "
    ),
    "\n"
  )
}


cat(
  "\nIMPORTANT:\n"
)


cat(
  paste0(
    "Priority is based on integrated directional RNA evidence; ",
    "it is not a new significance test.\n"
  )
)


cat(
  "\nNext canonical step:\n"
)


cat(
  "scripts/06_make_final_kinase_figure.R\n"
)


cat("\n========================================\n")