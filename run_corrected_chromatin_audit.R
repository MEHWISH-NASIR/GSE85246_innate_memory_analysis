# ============================================================
# run_corrected_chromatin_audit.R
#
# GSE85246 / GSE85245
# Corrected Day-6 promoter chromatin audit
#
# Final reproducible audit chain:
#   12a  audit Day-6 chromatin files
#   12c  BigWig global QC
#   12d  quantify TSS +/- 2 kb promoters
#   12e  correct RPMI_d6_rep2 NotNormalized tracks
#   12j  local Tobias-style kinase limma reproduction
#   12m  final transparent integration
#
# Diagnostic sensitivity scripts 12f-12l are intentionally
# NOT part of this final audit runner.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1)

pipeline <- data.frame(
  step = c("12a", "12c", "12d", "12e", "12j", "12m"),
  label = c(
    "Audit Day-6 chromatin files",
    "BigWig global QC",
    "Quantify Day-6 promoters",
    "Rescale RPMI rep2 promoter signals",
    "Local Tobias-style kinase limma",
    "Finalize corrected chromatin audit"
  ),
  script = c(
    "scripts/12a_audit_day6_chromatin_files.R",
    "scripts/12c_bigwig_global_qc.R",
    "scripts/12d_quantify_day6_promoters.R",
    "scripts/12e_rescale_RPMI_rep2_promoters.R",
    "scripts/12j_reproduce_Tobias_corrected_limma.R",
    "scripts/12m_finalize_chromatin_audit.R"
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Preflight
# ------------------------------------------------------------

required_root <- c(
  "scripts",
  "data",
  "results",
  "results/03_memory_kinases"
)

missing_root <- required_root[!file.exists(required_root)]

if (length(missing_root) > 0) {
  stop(
    paste(
      "Missing required project item(s):",
      paste(missing_root, collapse = ", ")
    ),
    call. = FALSE
  )
}

missing_scripts <- pipeline$script[
  !file.exists(pipeline$script)
]

if (length(missing_scripts) > 0) {
  stop(
    paste(
      "Missing audit script(s):",
      paste(missing_scripts, collapse = ", ")
    ),
    call. = FALSE
  )
}

# BigWigs are intentionally local / excluded from Git.
bw_files <- list.files(
  "data/chip",
  pattern = "\\.(bw|bigWig)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(bw_files) == 0) {
  stop(
    paste0(
      "No BigWig files found under data/chip/. ",
      "The corrected chromatin audit requires the local GEO tracks."
    ),
    call. = FALSE
  )
}

cat("\n========================================\n")
cat("CORRECTED DAY-6 CHROMATIN AUDIT\n")
cat("========================================\n\n")

cat(
  "BigWig files detected: ",
  length(bw_files),
  "\n\n",
  sep = ""
)

cat("Preflight syntax check:\n")

for (i in seq_len(nrow(pipeline))) {

  cat(
    sprintf(
      "  [%s] %-42s ... ",
      pipeline$step[i],
      basename(pipeline$script[i])
    )
  )

  z <- try(
    parse(
      file = pipeline$script[i],
      keep.source = FALSE
    ),
    silent = TRUE
  )

  if (inherits(z, "try-error")) {
    cat("FAILED\n")
    stop(
      paste(
        "Parse failed:",
        pipeline$script[i]
      ),
      call. = FALSE
    )
  }

  cat("OK\n")
}

# ------------------------------------------------------------
# Run in clean R processes
# ------------------------------------------------------------

rscript_bin <- Sys.which("Rscript")

if (!nzchar(rscript_bin)) {
  stop("Could not locate Rscript.", call. = FALSE)
}

log_dir <- "logs/corrected_chromatin_audit"
dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

status <- data.frame(
  step = pipeline$step,
  label = pipeline$label,
  script = pipeline$script,
  status = "NOT_STARTED",
  exit_status = NA_integer_,
  stringsAsFactors = FALSE
)

status_file <- file.path(
  "results/12_kinase_reevaluation",
  "12_run_status.csv"
)

write.csv(
  status,
  status_file,
  row.names = FALSE
)

for (i in seq_len(nrow(pipeline))) {

  step_id <- pipeline$step[i]
  script_file <- pipeline$script[i]

  cat("\n========================================\n")
  cat(
    "STEP ",
    step_id,
    " — ",
    pipeline$label[i],
    "\n",
    sep = ""
  )
  cat("========================================\n")

  log_file <- file.path(
    log_dir,
    paste0(step_id, ".log")
  )

  status$status[i] <- "RUNNING"

  write.csv(
    status,
    status_file,
    row.names = FALSE
  )

  exit_code <- system2(
    rscript_bin,
    args = c("--vanilla", script_file),
    stdout = log_file,
    stderr = log_file,
    wait = TRUE
  )

  status$exit_status[i] <- as.integer(exit_code)

  if (identical(as.integer(exit_code), 0L)) {

    status$status[i] <- "SUCCESS"

    cat("SUCCESS\n")

  } else {

    status$status[i] <- "FAILED"

    write.csv(
      status,
      status_file,
      row.names = FALSE
    )

    cat("\nFAILED\n")
    cat("Log:\n  ", log_file, "\n", sep = "")

    if (file.exists(log_file)) {
      x <- readLines(log_file, warn = FALSE)

      cat(
        paste(
          tail(x, 40),
          collapse = "\n"
        ),
        "\n"
      )
    }

    stop(
      paste(
        "Audit stopped at",
        script_file
      ),
      call. = FALSE
    )
  }

  write.csv(
    status,
    status_file,
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# Verify final products
# ------------------------------------------------------------

required_final <- c(
  "results/12_kinase_reevaluation/12j_all_kinase_corrected_limma.csv",
  "results/12_kinase_reevaluation/12j_target_results.csv",
  "results/12_kinase_reevaluation/12m_local_formal_promoter_results.csv",
  "results/12_kinase_reevaluation/12m_independent_corrected_benchmark.csv",
  "results/12_kinase_reevaluation/12m_FINAL_candidate_status.csv",
  "results/12_kinase_reevaluation/12m_reproducibility_summary.csv",
  "figures/12_kinase_reevaluation/12m_FINAL_candidate_summary.png"
)

missing_final <- required_final[
  !file.exists(required_final)
]

if (length(missing_final) > 0) {
  stop(
    paste(
      "Audit ran but final output(s) are missing:",
      paste(missing_final, collapse = ", ")
    ),
    call. = FALSE
  )
}

final_status <- read.csv(
  "results/12_kinase_reevaluation/12m_FINAL_candidate_status.csv",
  stringsAsFactors = FALSE
)

priority <- final_status$SYMBOL[
  final_status$Downstream_status == "PRIORITIZE"
]

cat("\n========================================\n")
cat("CORRECTED AUDIT COMPLETE\n")
cat("========================================\n\n")

cat(
  "Successful steps: ",
  sum(status$status == "SUCCESS"),
  " / ",
  nrow(status),
  "\n",
  sep = ""
)

cat(
  "Downstream priority: ",
  paste(priority, collapse = ", "),
  "\n",
  sep = ""
)

cat(
  "\nImportant:\n",
  "The downstream decision integrates local RNA evidence ",
  "with independently corrected Tobias promoter results.\n",
  "The exact Tobias FDR values are not claimed as locally ",
  "reproduced statistics.\n",
  sep = ""
)
