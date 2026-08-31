# ============================================================
# run_all.R
#
# GSE85246 innate-memory kinase/chromatin analysis
#
# Runs the canonical pipeline in order:
#   01 setup/data
#   02 initial LPS response
#   03 persistent Day-6 kinase memory
#   04 beta-glucan rescue + restimulation
#   05 integrated kinase trajectory
#   06 final kinase figures
#   07 epigenetic persistence / BigWig quantification
#   07b corrected condition-blind distal analysis
#   08 final chromatin figures
#
# Each step runs in its OWN clean R process using:
#   Rscript --vanilla
#
# Run from project root:
#   Rscript --vanilla run_all.R
#
# Logs:
#   logs/run_all/
#
# Status:
#   results/run_all_status.csv
# ============================================================

rm(list = ls())

options(stringsAsFactors = FALSE, warn = 1)

pipeline <- data.frame(
  step = c("01","02","03","04","05","06","07","07b","08"),
  label = c(
    "Setup and data audit",
    "Initial Day-1 LPS response",
    "Persistent Day-6 kinase memory",
    "Beta-glucan rescue and restimulation",
    "Integrated kinase trajectory",
    "Final kinase figures",
    "Epigenetic persistence",
    "Corrected unbiased distal chromatin",
    "Final chromatin figures"
  ),
  script = c(
    "scripts/01_setup_and_data.R",
    "scripts/02_initial_LPS_response.R",
    "scripts/03_persistent_kinase_memory.R",
    "scripts/04_BG_rescue_and_restimulation.R",
    "scripts/05_integrated_kinase_trajectory.R",
    "scripts/06_make_final_kinase_figure.R",
    "scripts/07_epigenetic_persistence.R",
    "scripts/07b_unbiased_distal_chromatin.R",
    "scripts/08_make_final_chromatin_figures.R"
  ),
  stringsAsFactors = FALSE
)

required_root_items <- c("scripts", "data", "results")
missing_root_items <- required_root_items[!file.exists(required_root_items)]

if (length(missing_root_items) > 0) {
  stop(
    paste0(
      "\nrun_all.R must be executed from the project root.\n",
      "Current working directory:\n  ",
      normalizePath(".", winslash = "/", mustWork = FALSE),
      "\n\nMissing expected project item(s):\n  ",
      paste(missing_root_items, collapse = "\n  ")
    ),
    call. = FALSE
  )
}

rscript_bin <- Sys.which("Rscript")

if (!nzchar(rscript_bin)) {
  fallback_candidates <- c(
    file.path(R.home("bin"), "Rscript.exe"),
    file.path(R.home("bin"), "Rscript")
  )
  fallback_candidates <- fallback_candidates[file.exists(fallback_candidates)]
  if (length(fallback_candidates) == 0) {
    stop("Could not locate Rscript executable.", call. = FALSE)
  }
  rscript_bin <- fallback_candidates[1]
}

cat("\n========================================\n")
cat("GSE85246 — CANONICAL RUN-ALL PIPELINE\n")
cat("========================================\n\n")

cat("Project root:\n  ",
    normalizePath(".", winslash = "/", mustWork = FALSE),
    "\n", sep = "")

cat("Rscript:\n  ",
    normalizePath(rscript_bin, winslash = "/", mustWork = FALSE),
    "\n\n", sep = "")

missing_scripts <- pipeline$script[!file.exists(pipeline$script)]

if (length(missing_scripts) > 0) {
  stop(
    paste0(
      "Missing canonical pipeline script(s):\n  ",
      paste(missing_scripts, collapse = "\n  ")
    ),
    call. = FALSE
  )
}

cat("All canonical scripts found.\n")

cat("\nPreflight syntax check:\n")

for (i in seq_len(nrow(pipeline))) {
  cat(sprintf(
    "  [%s] %-42s ... ",
    pipeline$step[i],
    basename(pipeline$script[i])
  ))

  parse_result <- try(
    parse(file = pipeline$script[i], keep.source = FALSE),
    silent = TRUE
  )

  if (inherits(parse_result, "try-error")) {
    cat("FAILED\n\n")
    cat(as.character(parse_result), "\n")
    stop(
      paste("Preflight parse failed:", pipeline$script[i]),
      call. = FALSE
    )
  } else {
    cat("PARSE OK\n")
  }
}

cat("\nPreflight complete: all scripts parse successfully.\n")

log_dir <- "logs/run_all"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
this_log_dir <- file.path(log_dir, run_id)
dir.create(this_log_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nRun ID: ", run_id, "\n", sep = "")
cat("Logs: ", this_log_dir, "\n\n", sep = "")

status <- data.frame(
  step = pipeline$step,
  label = pipeline$label,
  script = pipeline$script,
  status = rep("NOT_STARTED", nrow(pipeline)),
  started = rep(NA_character_, nrow(pipeline)),
  finished = rep(NA_character_, nrow(pipeline)),
  elapsed_minutes = rep(NA_real_, nrow(pipeline)),
  exit_status = rep(NA_integer_, nrow(pipeline)),
  log_file = rep(NA_character_, nrow(pipeline)),
  stringsAsFactors = FALSE
)

write_status <- function() {
  write.csv(status, "results/run_all_status.csv", row.names = FALSE)
}

write_status()

print_log_tail <- function(log_file, n = 25) {
  if (!file.exists(log_file)) return(invisible(NULL))
  x <- readLines(log_file, warn = FALSE)
  if (length(x) == 0) return(invisible(NULL))
  cat(paste(tail(x, n), collapse = "\n"), "\n")
  invisible(NULL)
}

pipeline_start <- Sys.time()

for (i in seq_len(nrow(pipeline))) {
  step_id <- pipeline$step[i]
  step_label <- pipeline$label[i]
  script_file <- pipeline$script[i]

  log_file <- file.path(
    this_log_dir,
    paste0(
      step_id, "_",
      tools::file_path_sans_ext(basename(script_file)),
      ".log"
    )
  )

  status$log_file[i] <- log_file
  status$status[i] <- "RUNNING"
  status$started[i] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  write_status()

  cat("\n========================================\n")
  cat("STEP ", step_id, " — ", step_label, "\n", sep = "")
  cat("========================================\n")
  cat("Script: ", script_file, "\n", sep = "")
  cat("Log:    ", log_file, "\n\n", sep = "")

  step_start <- Sys.time()

  exit_code <- system2(
    command = rscript_bin,
    args = c("--vanilla", script_file),
    stdout = log_file,
    stderr = log_file,
    wait = TRUE
  )

  step_end <- Sys.time()

  elapsed_minutes <- as.numeric(
    difftime(step_end, step_start, units = "mins")
  )

  status$finished[i] <- format(step_end, "%Y-%m-%d %H:%M:%S")
  status$elapsed_minutes[i] <- round(elapsed_minutes, 3)
  status$exit_status[i] <- as.integer(exit_code)

  if (identical(as.integer(exit_code), 0L)) {
    status$status[i] <- "SUCCESS"
    write_status()

    cat(sprintf("SUCCESS — %.2f minutes\n", elapsed_minutes))
    cat("\nLast lines from this step:\n")
    cat("----------------------------------------\n")
    print_log_tail(log_file, n = 12)
    cat("----------------------------------------\n")

  } else {
    status$status[i] <- "FAILED"
    write_status()

    cat(sprintf("\nFAILED — exit status %s\n", exit_code))
    cat("\nLast lines from failed step:\n")
    cat("----------------------------------------\n")
    print_log_tail(log_file, n = 50)
    cat("----------------------------------------\n")

    cat(
      "\nPipeline stopped at:\n  ", script_file,
      "\n\nFull log:\n  ", log_file, "\n",
      sep = ""
    )

    stop(
      paste("run_all.R stopped because", script_file, "failed."),
      call. = FALSE
    )
  }
}

pipeline_end <- Sys.time()
total_minutes <- as.numeric(
  difftime(pipeline_end, pipeline_start, units = "mins")
)

write_status()

cat("\n========================================\n")
cat("RUN-ALL PIPELINE COMPLETE\n")
cat("========================================\n\n")

cat(
  "Successful steps: ",
  sum(status$status == "SUCCESS"),
  " / ",
  nrow(status),
  "\n",
  sep = ""
)

cat(sprintf("Total elapsed time: %.2f minutes\n", total_minutes))
cat("\nStatus table:\n  results/run_all_status.csv\n")
cat("\nRun logs:\n  ", this_log_dir, "\n", sep = "")

key_outputs <- c(
  "results/05_integrated_trajectory/05_integrated_kinase_trajectory.csv",
  "figures/06_final_kinase/06D_FINAL_integrated_kinase_summary.png",
  "results/07b_unbiased_distal_chromatin/07b_corrected_gene_level_chromatin_summary.csv",
  "results/08_final_chromatin/08_TOP_integrated_candidates.csv",
  "figures/08_final_chromatin/08D_FINAL_RNA_chromatin_integration.png",
  "results/08_final_chromatin/08_TOP_candidate_IGV_loci_hg19.csv"
)

output_check <- data.frame(
  file = key_outputs,
  exists = file.exists(key_outputs),
  stringsAsFactors = FALSE
)

write.csv(
  output_check,
  "results/run_all_final_output_check.csv",
  row.names = FALSE
)

cat("\nKey final-output check:\n")

for (i in seq_len(nrow(output_check))) {
  cat(
    ifelse(output_check$exists[i], "  [OK]      ", "  [MISSING] "),
    output_check$file[i],
    "\n",
    sep = ""
  )
}

if (all(output_check$exists)) {
  cat("\nAll key final outputs are present.\n")
} else {
  warning(
    paste(
      "Pipeline finished, but one or more",
      "key final outputs are missing."
    )
  )
}

top_file <- "results/08_final_chromatin/08_TOP_integrated_candidates.csv"

if (file.exists(top_file)) {
  top <- read.csv(
    top_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  cat("\nFinal top integrated candidates: ", nrow(top), "\n", sep = "")

  if (nrow(top) > 0 && "SYMBOL" %in% names(top)) {
    cat("  ", paste(top$SYMBOL, collapse = ", "), "\n", sep = "")
  }
}

cat(
  "\nNext project stage:\n",
  "  final report / README / GitHub packaging\n",
  sep = ""
)

cat("\n========================================\n")
