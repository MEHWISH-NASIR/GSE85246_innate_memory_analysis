# ============================================================
# 13n — TCGA-KIRC immune/inflammatory marker extraction
#
# Purpose:
# Control bulk RNA-seq JAK3/EPHB2 signals for immune-cell content.
#
# Extract all markers from each STAR-count file in one pass.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
options(timeout = 90)

metadata_file <-
  "results/13_TCGA_JAK3/13a_TCGA_STARCounts_metadata.csv"

outdir <-
  "results/13_TCGA_priority_genes"

outfile <-
  file.path(
    outdir,
    "13n_KIRC_immune_markers_expression.csv"
  )

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Target immune markers
# ------------------------------------------------------------

markers <- c(
  "PTPRC",
  "CD3D",
  "CD3E",
  "MS4A1",
  "CD79A",
  "NKG7",
  "LST1",
  "TYROBP",
  "FCER1G",
  "CD68"
)

# ------------------------------------------------------------
# Metadata
# ------------------------------------------------------------

meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kirc <- meta[
  meta$TCGA_project == "TCGA-KIRC" &
  meta$sample_type %in%
    c(
      "Primary Tumor",
      "Solid Tissue Normal"
    ),
]

kirc <- kirc[
  !duplicated(kirc$file_id),
]

kirc$sample_barcode <- kirc$cases

kirc$patient_id <- substr(
  kirc$cases,
  1,
  12
)

cat("\n========================================\n")
cat("KIRC IMMUNE MARKER EXTRACTION\n")
cat("========================================\n\n")

cat("Files:", nrow(kirc), "\n")
cat("Markers:", length(markers), "\n\n")

print(
  table(kirc$sample_type)
)

# ------------------------------------------------------------
# Resume support
#
# A file is complete only if all marker genes were extracted OK.
# ------------------------------------------------------------

if (file.exists(outfile)) {

  old <- read.csv(
    outfile,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

} else {

  old <- data.frame()
}

done <- character(0)

if (nrow(old) > 0) {

  ok <- old[
    old$extraction_status == "OK",
  ]

  counts <- table(
    ok$file_id
  )

  done <- names(
    counts[counts >= length(markers)]
  )
}

todo <- kirc[
  !(kirc$file_id %in% done),
]

cat(
  "\nPreviously completed files:",
  length(done),
  "\n"
)

cat(
  "Files remaining:",
  nrow(todo),
  "\n\n"
)

if (nrow(todo) == 0) {

  cat("Nothing left to extract.\n")
  quit(save = "no")
}

# ------------------------------------------------------------
# Download helper
# ------------------------------------------------------------

download_with_retry <- function(
  url,
  dest,
  attempts = 3
) {

  for (a in seq_len(attempts)) {

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
      a,
      "of",
      attempts,
      "\n"
    )

    Sys.sleep(2 * a)
  }

  FALSE
}

# ------------------------------------------------------------
# Save checkpoint
# ------------------------------------------------------------

save_checkpoint <- function(
  old,
  new_rows
) {

  if (length(new_rows) == 0) {
    return(old)
  }

  ndf <- do.call(
    rbind,
    new_rows
  )

  if (nrow(old) > 0) {

    all <- rbind(
      old,
      ndf
    )

  } else {

    all <- ndf
  }

  key <- paste(
    all$file_id,
    all$gene_name,
    sep = "__"
  )

  all <- all[
    !duplicated(
      key,
      fromLast = TRUE
    ),
  ]

  write.csv(
    all,
    outfile,
    row.names = FALSE
  )

  all
}

# ------------------------------------------------------------
# Extraction
# ------------------------------------------------------------

new_rows <- list()

for (i in seq_len(nrow(todo))) {

  row <- todo[i, ]

  cat(
    sprintf(
      "[%03d/%03d] %s  %s\n",
      i,
      nrow(todo),
      row$sample_type,
      row$sample_barcode
    )
  )

  url <- paste0(
    "https://api.gdc.cancer.gov/data/",
    row$file_id
  )

  tmp <- tempfile(
    fileext = ".tsv"
  )

  downloaded <- download_with_retry(
    url,
    tmp,
    attempts = 3
  )

  # ----------------------------------------------------------
  # Download failure
  # ----------------------------------------------------------

  if (!downloaded) {

    cat("    DOWNLOAD FAILED\n")

    for (gene in markers) {

      new_rows[[
        length(new_rows) + 1
      ]] <- data.frame(

        project =
          "TCGA-KIRC",

        file_id =
          row$file_id,

        sample_barcode =
          row$sample_barcode,

        patient_id =
          row$patient_id,

        sample_type =
          row$sample_type,

        gene_name =
          gene,

        gene_id =
          NA_character_,

        TPM =
          NA_real_,

        extraction_status =
          "DOWNLOAD_FAILED",

        stringsAsFactors = FALSE
      )
    }

    next
  }

  # ----------------------------------------------------------
  # Read STAR-count table
  # ----------------------------------------------------------

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

  # ----------------------------------------------------------
  # Extract all markers
  # ----------------------------------------------------------

  for (gene in markers) {

    hit <- which(
      dat$gene_name == gene
    )

    if (length(hit) != 1) {

      cat(
        "    ",
        gene,
        " match count = ",
        length(hit),
        "\n",
        sep = ""
      )

      new_rows[[
        length(new_rows) + 1
      ]] <- data.frame(

        project =
          "TCGA-KIRC",

        file_id =
          row$file_id,

        sample_barcode =
          row$sample_barcode,

        patient_id =
          row$patient_id,

        sample_type =
          row$sample_type,

        gene_name =
          gene,

        gene_id =
          NA_character_,

        TPM =
          NA_real_,

        extraction_status =
          "GENE_NOT_UNIQUE",

        stringsAsFactors = FALSE
      )

      next
    }

    h <- hit[1]

    val <- suppressWarnings(
      as.numeric(
        dat$tpm_unstranded[h]
      )
    )

    new_rows[[
      length(new_rows) + 1
    ]] <- data.frame(

      project =
        "TCGA-KIRC",

      file_id =
        row$file_id,

      sample_barcode =
        row$sample_barcode,

      patient_id =
        row$patient_id,

      sample_type =
        row$sample_type,

      gene_name =
        gene,

      gene_id =
        dat$gene_id[h],

      TPM =
        val,

      extraction_status =
        "OK",

      stringsAsFactors = FALSE
    )
  }

  cat("    immune markers extracted\n")

  # ----------------------------------------------------------
  # Checkpoint every 10 files
  # ----------------------------------------------------------

  if (
    i %% 10 == 0 ||
    i == nrow(todo)
  ) {

    old <- save_checkpoint(
      old,
      new_rows
    )

    new_rows <- list()

    cat(
      "    CHECKPOINT:",
      length(unique(old$file_id)),
      "files saved\n"
    )
  }
}

# ------------------------------------------------------------
# Final audit
# ------------------------------------------------------------

res <- read.csv(
  outfile,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("FINAL KIRC IMMUNE-MARKER AUDIT\n")
cat("========================================\n\n")

cat(
  "Unique files:",
  length(unique(res$file_id)),
  "\n"
)

cat("\nExtraction status by marker:\n")

print(
  table(
    res$gene_name,
    res$extraction_status,
    useNA = "ifany"
  )
)

cat("\nSample counts by marker:\n")

print(
  table(
    res$gene_name,
    res$sample_type
  )
)

cat("\nSaved:\n ", outfile, "\n")

cat("\n========================================\n")
cat("13n COMPLETE\n")
cat("========================================\n")
