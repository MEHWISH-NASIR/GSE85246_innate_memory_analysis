# ============================================================
# 09f_chromatin_strength_audit.R
#
# Quantitative audit of the 09e chromatin result.
#
# Questions:
# 1. How large are the Day1 and Day6 chromatin differences?
# 2. Are the Day6 effects seen in both donors 73 and 74?
# 3. Does beta-glucan oppose the LPS state in both donors?
#
# No new loci are selected here.
# No BigWig files are reread.
# ============================================================

rm(list = ls())

signal_file <-
  "results/09_MAP3K8_validation/09e_target_all_region_BigWig_signal.csv"

selected_file <-
  "results/09_MAP3K8_validation/09e_condition_blind_selected_distal_windows.csv"

outdir <-
  "results/09_MAP3K8_validation"

if (!file.exists(signal_file)) {
  stop("Missing 09e target BigWig signal file.")
}

if (!file.exists(selected_file)) {
  stop("Missing 09e selected distal-window file.")
}

sig <- read.csv(
  signal_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

selected <- read.csv(
  selected_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

targets <- c(
  "MAP3K7CL",
  "TTC39B",
  "TNIP3"
)

marks <- c(
  "H3K27ac",
  "H3K4me1"
)

cat("\n========================================\n")
cat("09f — CHROMATIN STRENGTH AUDIT\n")
cat("========================================\n\n")

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------

analyse_region_mark <- function(
    symbol,
    region_id,
    region_type,
    mark
) {

  x <- sig[
    sig$region_id == region_id,
    ,
    drop = FALSE
  ]

  if (nrow(x) != 1) {
    stop(
      paste(
        "Expected one row for",
        region_id,
        "found",
        nrow(x)
      )
    )
  }

  # Day 1
  d1_rpmi_values <- c(
    x[[paste0("D1_RPMI_r1_", mark)]],
    x[[paste0("D1_RPMI_r2_", mark)]]
  )

  d1_lps_values <- c(
    x[[paste0("D1_LPS_r1_", mark)]],
    x[[paste0("D1_LPS_r2_", mark)]]
  )

  d1_rpmi <- mean(d1_rpmi_values)
  d1_lps  <- mean(d1_lps_values)

  d1_delta <- d1_lps - d1_rpmi

  # Day 6 donor 73
  d73_rpmi <-
    x[[paste0("D6_RPMI_d73_", mark)]]

  d73_lps <-
    x[[paste0("D6_LPS_d73_", mark)]]

  d73_bg <-
    x[[paste0("D6_RESCUE_d73_", mark)]]

  # Day 6 donor 74
  d74_rpmi <-
    x[[paste0("D6_RPMI_d74_", mark)]]

  d74_lps <-
    x[[paste0("D6_LPS_d74_", mark)]]

  d74_bg <-
    x[[paste0("D6_RESCUE_d74_", mark)]]

  # Day6 group means
  d6_rpmi <- mean(
    c(d73_rpmi, d74_rpmi)
  )

  d6_lps <- mean(
    c(d73_lps, d74_lps)
  )

  d6_bg <- mean(
    c(d73_bg, d74_bg)
  )

  d6_delta <-
    d6_lps - d6_rpmi

  bg_shift <-
    d6_bg - d6_lps

  # Donor-specific LPS effects
  d73_delta <-
    d73_lps - d73_rpmi

  d74_delta <-
    d74_lps - d74_rpmi

  # Donor-specific BG shifts
  d73_bg_shift <-
    d73_bg - d73_lps

  d74_bg_shift <-
    d74_bg - d74_lps

  # Persistence criterion
  cross_time_same_direction <-
    d1_delta * d6_delta > 0

  # Does each Day6 donor agree with the mean LPS direction?
  day6_both_donors_consistent <-
    d73_delta * d6_delta > 0 &&
    d74_delta * d6_delta > 0

  # Does BG oppose LPS in each donor?
  bg_opposes_d73 <-
    d73_bg_shift * d73_delta < 0

  bg_opposes_d74 <-
    d74_bg_shift * d74_delta < 0

  bg_opposes_both <-
    bg_opposes_d73 &&
    bg_opposes_d74

  # Is BG closer to RPMI for each donor?
  bg_closer_d73 <-
    abs(d73_bg - d73_rpmi) <
    abs(d73_lps - d73_rpmi)

  bg_closer_d74 <-
    abs(d74_bg - d74_rpmi) <
    abs(d74_lps - d74_rpmi)

  bg_closer_both <-
    bg_closer_d73 &&
    bg_closer_d74

  # Original 07b-style mean-based support
  mean_based_support <-
    cross_time_same_direction &&
    (bg_shift * d6_delta < 0) &&
    (
      abs(d6_bg - d6_rpmi) <
      abs(d6_lps - d6_rpmi)
    )

  data.frame(
    SYMBOL = symbol,
    region_type = region_type,
    region_id = region_id,
    mark = mark,

    Day1_RPMI_mean = d1_rpmi,
    Day1_LPS_mean = d1_lps,
    Day1_delta = d1_delta,

    Day6_RPMI_mean = d6_rpmi,
    Day6_LPS_mean = d6_lps,
    Day6_BG_mean = d6_bg,
    Day6_delta = d6_delta,
    BG_shift = bg_shift,

    D73_LPS_minus_RPMI = d73_delta,
    D74_LPS_minus_RPMI = d74_delta,

    D73_BG_minus_LPS = d73_bg_shift,
    D74_BG_minus_LPS = d74_bg_shift,

    cross_time_same_direction =
      cross_time_same_direction,

    Day6_both_donors_consistent =
      day6_both_donors_consistent,

    BG_opposes_D73 =
      bg_opposes_d73,

    BG_opposes_D74 =
      bg_opposes_d74,

    BG_opposes_both =
      bg_opposes_both,

    BG_closer_D73 =
      bg_closer_d73,

    BG_closer_D74 =
      bg_closer_d74,

    BG_closer_both =
      bg_closer_both,

    Mean_based_09e_support =
      mean_based_support,

    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Analyse promoter + frozen distal region
# ------------------------------------------------------------

results <- list()

k <- 1

for (g in targets) {

  promoter_id <-
    paste0(
      g,
      "_PROMOTER"
    )

  distal_id <-
    selected$region_id[
      selected$SYMBOL == g
    ]

  if (length(distal_id) != 1) {
    stop(
      paste(
        "Could not identify frozen distal region for",
        g
      )
    )
  }

  for (m in marks) {

    results[[k]] <-
      analyse_region_mark(
        g,
        promoter_id,
        "promoter",
        m
      )

    k <- k + 1

    results[[k]] <-
      analyse_region_mark(
        g,
        distal_id,
        "fixed_nonpromoter",
        m
      )

    k <- k + 1
  }
}

audit <-
  do.call(
    rbind,
    results
  )

write.csv(
  audit,
  file.path(
    outdir,
    "09f_chromatin_strength_and_donor_audit.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Print the most important information
# ------------------------------------------------------------

cat("All region/mark combinations:\n\n")

print(
  audit[
    ,
    c(
      "SYMBOL",
      "region_type",
      "mark",
      "Day1_delta",
      "Day6_delta",
      "BG_shift",
      "D73_LPS_minus_RPMI",
      "D74_LPS_minus_RPMI",
      "D73_BG_minus_LPS",
      "D74_BG_minus_LPS",
      "Mean_based_09e_support"
    )
  ],
  row.names = FALSE
)

cat("\n========================================\n")
cat("09e-SUPPORTED CHROMATIN ONLY\n")
cat("========================================\n\n")

supported <-
  audit[
    audit$Mean_based_09e_support,
    ,
    drop = FALSE
  ]

print(
  supported[
    ,
    c(
      "SYMBOL",
      "region_type",
      "mark",
      "Day1_delta",
      "Day6_delta",
      "BG_shift",
      "Day6_both_donors_consistent",
      "BG_opposes_both",
      "BG_closer_both"
    )
  ],
  row.names = FALSE
)

cat("\n========================================\n")
cat("STRICT DONOR-CONSISTENT SUPPORT\n")
cat("========================================\n\n")

strict <-
  supported[
    supported$Day6_both_donors_consistent &
    supported$BG_opposes_both &
    supported$BG_closer_both,
    ,
    drop = FALSE
  ]

if (nrow(strict) == 0) {

  cat(
    "No 09e-supported region passes all donor-level checks.\n"
  )

} else {

  print(
    strict[
      ,
      c(
        "SYMBOL",
        "region_type",
        "mark",
        "Day1_delta",
        "Day6_delta",
        "BG_shift"
      )
    ],
    row.names = FALSE
  )
}

cat(
  "\n09e-supported region/mark combinations:",
  nrow(supported),
  "\n"
)

cat(
  "Strict donor-consistent combinations:",
  nrow(strict),
  "\n"
)

cat("\n========================================\n")
cat("09f COMPLETE\n")
cat("========================================\n")
