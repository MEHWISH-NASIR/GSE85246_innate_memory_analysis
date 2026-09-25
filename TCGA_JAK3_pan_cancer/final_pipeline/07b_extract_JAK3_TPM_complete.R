library(data.table)

files <- list.files(
  "GDCdata",
  pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("\nExact STAR-count files found:", length(files), "\n")

if (length(files) != 11505) {
  stop("Expected 11505 exact STAR-count files, found ", length(files))
}

results <- vector("list", length(files))

for (i in seq_along(files)) {

  f <- files[i]

  x <- fread(
    f,
    skip = 6,
    header = FALSE,
    showProgress = FALSE
  )

  hit <- x[V2 == "JAK3"]

  if (nrow(hit) != 1) {
    stop(
      "JAK3 row problem in file: ",
      f,
      " | rows found = ",
      nrow(hit)
    )
  }

  results[[i]] <- data.frame(
    sample = basename(f),
    gene = hit$V2,
    counts = hit$V4,
    TPM = hit$V7,
    stringsAsFactors = FALSE
  )

  if (i %% 500 == 0 || i == length(files)) {
    cat("Processed:", i, "/", length(files), "\n")
  }
}

JAK3_TPM <- rbindlist(results)

write.csv(
  JAK3_TPM,
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_TPM_complete.csv",
  row.names = FALSE
)

cat("\n========================================\n")
cat("JAK3 COMPLETE TPM EXTRACTION\n")
cat("========================================\n")
cat("Input STAR files:", length(files), "\n")
cat("Output rows:", nrow(JAK3_TPM), "\n")
cat("Unique samples:", uniqueN(JAK3_TPM$sample), "\n")
cat("Missing TPM:", sum(is.na(JAK3_TPM$TPM)), "\n")
cat("\nJAK3 TPM EXTRACTION COMPLETE\n")
