# ============================================================
# render_report.R
#
# Render the final GSE85246 analysis report.
#
# Run from project root:
#   Rscript --vanilla render_report.R
#
# Output:
#   docs/analysis_report.html
# ============================================================

rm(list = ls())

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop(
    "Package 'rmarkdown' is required.",
    call. = FALSE
  )
}

if (!rmarkdown::pandoc_available()) {
  stop(
    "Pandoc is not available to rmarkdown.",
    call. = FALSE
  )
}

input_file <- "analysis_report.Rmd"
output_dir <- "docs"
output_file <- "analysis_report.html"

if (!file.exists(input_file)) {
  stop(
    paste("Missing report source:", input_file),
    call. = FALSE
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n========================================\n")
cat("RENDERING FINAL ANALYSIS REPORT\n")
cat("========================================\n\n")

cat("Input:  ", input_file, "\n", sep = "")
cat(
  "Output: ",
  file.path(output_dir, output_file),
  "\n\n",
  sep = ""
)

result <- rmarkdown::render(
  input = input_file,
  output_format = "html_document",
  output_file = output_file,
  output_dir = output_dir,
  envir = new.env(parent = globalenv()),
  clean = TRUE,
  quiet = FALSE
)

cat("\n========================================\n")
cat("REPORT RENDER COMPLETE\n")
cat("========================================\n\n")

cat(
  "Created:\n  ",
  normalizePath(
    result,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "\nOpen from Git Bash with:\n",
  '  start "" "docs/analysis_report.html"\n',
  sep = ""
)

cat("\n========================================\n")
