library(TCGAbiolinks)
library(dplyr)

expr <- read.csv(
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_final_TPM_metadata_complete.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

projects <- sort(unique(expr$project))

cat("\n========================================\n")
cat("COMPLETE TCGA CLINICAL FETCH\n")
cat("========================================\n")

cat("\nProjects to fetch:", length(projects), "\n")
print(projects)

if (length(projects) != 33) {
  stop("Expected 33 TCGA projects; found ", length(projects))
}

clinical_all <- list()
failed_projects <- character()

for (p in projects) {

  cat("\nProcessing:", p, "\n")

  x <- tryCatch(
    {
      GDCquery_clinic(
        project = p,
        type = "clinical"
      )
    },
    error = function(e) {
      cat("FAILED:", p, "|", conditionMessage(e), "\n")
      failed_projects <<- c(failed_projects, p)
      NULL
    }
  )

  if (!is.null(x)) {

    # Remove list columns before writing/combining
    x <- x[, !sapply(x, is.list), drop = FALSE]

    x$project <- p

    clinical_all[[p]] <- x

    cat(
      "Clinical rows:",
      nrow(x),
      "\n"
    )
  }
}

clinical <- bind_rows(clinical_all)

cat("\n========================================\n")
cat("CLINICAL FETCH SUMMARY\n")
cat("========================================\n")

cat("\nSuccessful projects:",
    length(clinical_all), "\n")

cat("Failed projects:",
    length(failed_projects), "\n")

if (length(failed_projects) > 0) {
  cat("Failed project IDs:\n")
  print(failed_projects)
}

cat("\nTotal clinical rows:",
    nrow(clinical), "\n")

cat(
  "Unique project + patient combinations:",
  n_distinct(
    paste(
      clinical$project,
      clinical$submitter_id
    )
  ),
  "\n"
)

cat(
  "Projects represented:",
  n_distinct(clinical$project),
  "\n"
)

write.csv(
  clinical,
  "TCGA_JAK3_pan_cancer/results/JAK3_TCGA_clinical_complete.csv",
  row.names = FALSE
)

cat("\nSaved:\n")
cat(
  "TCGA_JAK3_pan_cancer/results/",
  "JAK3_TCGA_clinical_complete.csv\n",
  sep = ""
)

cat("\nCOMPLETE CLINICAL FETCH FINISHED\n")
