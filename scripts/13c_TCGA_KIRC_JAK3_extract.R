# ============================================================
# 13c — TCGA-KIRC JAK3 EXPRESSION EXTRACTION
#
# Downloads one STAR-count file at a time, extracts JAK3,
# then deletes the temporary file.
#
# Samples:
#   Primary Tumor
#   Solid Tissue Normal
#
# Measures retained:
#   raw unstranded count
#   TPM
#   FPKM
#   FPKM-UQ
#
# Resumable: previously completed file IDs are skipped.
# ============================================================

rm(list = ls())

metadata_file <-
  "results/13_TCGA_JAK3/13a_TCGA_STARCounts_metadata.csv"

outfile <-
  "results/13_TCGA_JAK3/13c_KIRC_JAK3_expression.csv"

if (!file.exists(metadata_file)) {
  stop("Missing 13a metadata file.")
}

meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# Select KIRC tumor + normal only
# ------------------------------------------------------------

kirc <- meta[
  meta$TCGA_project == "TCGA-KIRC" &
  meta$sample_type %in%
    c("Primary Tumor", "Solid Tissue Normal"),
]

kirc <- kirc[
  !duplicated(kirc$file_id),
]

cat("\n========================================\n")
cat("13c — TCGA-KIRC JAK3 EXTRACTION\n")
cat("========================================\n\n")

cat("Files selected:", nrow(kirc), "\n\n")

print(
  table(kirc$sample_type)
)

if (nrow(kirc) != 613) {
  warning(
    paste(
      "Expected 613 KIRC tumor/normal files;",
      "found",
      nrow(kirc)
    )
  )
}

# ------------------------------------------------------------
# Resume support
# ------------------------------------------------------------

if (file.exists(outfile)) {

  completed <- read.csv(
    outfile,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  done_ids <- unique(
    completed$file_id
  )

  cat(
    "\nExisting completed records:",
    length(done_ids),
    "\n"
  )

} else {

  completed <- NULL
  done_ids <- character(0)
}

todo <- kirc[
  !(kirc$file_id %in% done_ids),
]

cat(
  "Files remaining:",
  nrow(todo),
  "\n\n"
)

if (nrow(todo) == 0) {
  cat("Nothing to do — all files already extracted.\n")
  quit(save = "no")
}

new_results <- list()

# ------------------------------------------------------------
# Helper: download with retries
# ------------------------------------------------------------

download_with_retry <- function(
    url,
    dest,
    retries = 3
) {

  for (attempt in seq_len(retries)) {

    ok <- tryCatch({

      download.file(
        url,
        destfile = dest,
        mode = "wb",
        quiet = TRUE
      )

      file.exists(dest) &&
        file.info(dest)$size > 0

    }, error = function(e) {

      cat(
        "    download error:",
        conditionMessage(e),
        "\n"
      )

      FALSE
    })

    if (isTRUE(ok)) {
      return(TRUE)
    }

    if (file.exists(dest)) {
      unlink(dest)
    }

    cat(
      "    retry",
      attempt,
      "of",
      retries,
      "\n"
    )

    Sys.sleep(3)
  }

  FALSE
}

# ------------------------------------------------------------
# Extract JAK3 one sample at a time
# ------------------------------------------------------------

for (i in seq_len(nrow(todo))) {

  row <- todo[i, ]

  cat(
    sprintf(
      "[%03d/%03d] %s  %s\n",
      i,
      nrow(todo),
      row$sample_type,
      row$cases
    )
  )

  tmp <- tempfile(
    fileext = ".tsv"
  )

  url <- paste0(
    "https://api.gdc.cancer.gov/data/",
    row$file_id
  )

  ok <- download_with_retry(
    url,
    tmp,
    retries = 3
  )

  if (!ok) {

    cat("    FAILED after retries.\n")

    new_results[[length(new_results) + 1]] <-
      data.frame(
        project = "TCGA-KIRC",
        sample_barcode = row$cases,
        patient_id =
          substr(row$cases, 1, 12),
        sample_type = row$sample_type,
        file_id = row$file_id,
        gene_id = NA_character_,
        gene_name = "JAK3",
        unstranded_count = NA_real_,
        TPM = NA_real_,
        FPKM = NA_real_,
        FPKM_UQ = NA_real_,
        extraction_status = "DOWNLOAD_FAILED",
        stringsAsFactors = FALSE
      )

    next
  }

  dat <- tryCatch(

    read.delim(
      tmp,
      header = TRUE,
      sep = "\t",
      comment.char = "#",
      fill = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),

    error = function(e) NULL
  )

  unlink(tmp)

  if (is.null(dat)) {

    cat("    TABLE READ FAILED\n")

    next
  }

  # Remove Ensembl version suffix
  ens_base <- sub(
    "\\..*$",
    "",
    dat$gene_id
  )

  hit <- which(
    dat$gene_name == "JAK3" &
    ens_base == "ENSG00000105639"
  )

  if (length(hit) != 1) {

    cat(
      "    JAK3 match count:",
      length(hit),
      "\n"
    )

    next
  }

  j <- dat[hit, ]

  result <- data.frame(

    project = "TCGA-KIRC",

    sample_barcode =
      row$cases,

    patient_id =
      substr(
        row$cases,
        1,
        12
      ),

    sample_type =
      row$sample_type,

    file_id =
      row$file_id,

    gene_id =
      j$gene_id,

    gene_name =
      j$gene_name,

    unstranded_count =
      as.numeric(
        j$unstranded
      ),

    TPM =
      as.numeric(
        j$tpm_unstranded
      ),

    FPKM =
      as.numeric(
        j$fpkm_unstranded
      ),

    FPKM_UQ =
      as.numeric(
        j$fpkm_uq_unstranded
      ),

    extraction_status =
      "OK",

    stringsAsFactors = FALSE
  )

  new_results[[length(new_results) + 1]] <-
    result

  cat(
    sprintf(
      "    JAK3 TPM = %.4f\n",
      result$TPM
    )
  )

  # ----------------------------------------------------------
  # Checkpoint every 25 successful/attempted records
  # ----------------------------------------------------------

  if (
    i %% 25 == 0 ||
    i == nrow(todo)
  ) {

    current_new <- do.call(
      rbind,
      new_results
    )

    if (is.null(completed)) {

      current_all <- current_new

    } else {

      current_all <- rbind(
        completed,
        current_new
      )
    }

    current_all <- current_all[
      !duplicated(
        current_all$file_id,
        fromLast = TRUE
      ),
    ]

    write.csv(
      current_all,
      outfile,
      row.names = FALSE
    )

    cat(
      "    CHECKPOINT:",
      nrow(current_all),
      "records saved\n"
    )
  }

  # Light pause between GDC requests
  Sys.sleep(0.1)
}

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

final <- read.csv(
  outfile,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("FINAL EXTRACTION SUMMARY\n")
cat("========================================\n\n")

cat(
  "Records:",
  nrow(final),
  "\n"
)

cat(
  "Successful:",
  sum(final$extraction_status == "OK"),
  "\n"
)

cat(
  "Failed:",
  sum(final$extraction_status != "OK"),
  "\n\n"
)

cat("Sample types:\n")

print(
  table(
    final$sample_type,
    useNA = "ifany"
  )
)

cat("\nJAK3 TPM summary:\n")

print(
  tapply(
    final$TPM,
    final$sample_type,
    summary
  )
)

cat("\nSaved:\n")
cat(" ", outfile, "\n")

cat("\n========================================\n")
cat("13c COMPLETE\n")
cat("========================================\n")
