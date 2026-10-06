
library(data.table)

# ==========================================
# Configuration
# ==========================================

root <- "TCGA_multikinase_analysis"
outdir <- file.path(root, "results/expression")
checkpoint_dir <- file.path(outdir, "checkpoints")

dir.create(
  checkpoint_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

genes <- c(
  "JAK3", "MET", "HCK", "LYN", "DDR1",
  "IRAK2", "FJX1", "EPHB1", "RPS6KA2",
  "EPHB2"
)

meta <- fread(
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_GDC_metadata.csv"
)

files <- sort(list.files(
  "D:/TCGA_GDCdata",
  pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))

file_ids <- basename(dirname(files))

stopifnot(
  length(files) == 11505,
  !anyDuplicated(file_ids),
  setequal(file_ids, unique(meta$file_id))
)

# Test mode processes 10 files only.
test_mode <- "--test" %in% commandArgs(
  trailingOnly = TRUE
)

if (test_mode) {
  files <- files[1:10]
  file_ids <- file_ids[1:10]
}

# ==========================================
# Extraction function
# ==========================================

extract_one <- function(i) {

  x <- fread(
    files[i],
    skip = 6,
    header = FALSE,
    showProgress = FALSE
  )

  hit <- x[V2 %in% genes]

  if (
    nrow(hit) != length(genes) ||
    anyDuplicated(hit$V2) ||
    !setequal(hit$V2, genes)
  ) {
    stop(
      "Gene annotation problem: ",
      files[i]
    )
  }

  if (anyNA(hit$V7)) {
    stop("Missing TPM: ", files[i])
  }

  data.table(
    file_id = file_ids[i],
    gene_id = hit$V1,
    gene = hit$V2,
    counts = hit$V4,
    TPM = hit$V7
  )
}

# ==========================================
# Processing
# ==========================================

if (test_mode) {

  result <- rbindlist(
    lapply(seq_along(files), extract_one)
  )

  print(result[, .(
    Files = uniqueN(file_id),
    Missing_TPM = sum(is.na(TPM))
  ), by = gene])

  cat("\nTEST COMPLETE\n")

} else {

  chunks <- split(
    seq_along(files),
    ceiling(seq_along(files) / 100)
  )

  for (k in seq_along(chunks)) {

    output <- file.path(
      checkpoint_dir,
      sprintf("chunk_%04d.rds", k)
    )

    if (file.exists(output)) {

      saved <- readRDS(output)

      if (
        nrow(saved) !=
          length(chunks[[k]]) * length(genes) ||
        !setequal(
          unique(saved$file_id),
          file_ids[chunks[[k]]]
        )
      ) {
        stop("Invalid checkpoint: ", output)
      }

      cat("Skipping completed chunk:", k, "\n")
      next
    }

    result <- rbindlist(
      lapply(chunks[[k]], extract_one)
    )

    saveRDS(result, output)

    cat(
      "Completed chunk:",
      k, "/", length(chunks),
      "| Files processed:",
      max(chunks[[k]]),
      "\n"
    )
  }

  outputs <- sort(list.files(
    checkpoint_dir,
    pattern = "^chunk_[0-9]+\\.rds$",
    full.names = TRUE
  ))

  stopifnot(length(outputs) == length(chunks))

  final <- rbindlist(lapply(outputs, readRDS))

  stopifnot(
    nrow(final) == 11505 * length(genes),
    uniqueN(final$file_id) == 11505
  )

  fwrite(
    final,
    file.path(
      outdir,
      "01_multikinase_TPM_complete.csv"
    )
  )

  cat("\n================================\n")
  cat("MULTIKINASE EXTRACTION COMPLETE\n")
  cat("================================\n")
  cat("Samples:", uniqueN(final$file_id), "\n")
  cat("Genes:", uniqueN(final$gene), "\n")
  cat("Rows:", nrow(final), "\n")
  cat("Missing TPM:", sum(is.na(final$TPM)), "\n")
}

