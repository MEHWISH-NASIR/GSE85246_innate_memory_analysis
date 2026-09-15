# ============================================================
# 14a — TCGA-BLCA combined extraction
# JAK3 + EPHB2 from the same STAR-count files
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
options(timeout = 90)

manifest_file <-
  "results/13_TCGA_JAK3/13e_pancancer_JAK3_expression.csv"

outdir <- "results/14_TCGA_BLCA"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <-
  file.path(
    outdir,
    "14a_BLCA_JAK3_EPHB2_expression.csv"
  )

if (!file.exists(manifest_file)) {
  stop("Missing existing TCGA manifest.")
}

# ------------------------------------------------------------
# Existing BLCA file manifest
# ------------------------------------------------------------

m <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

m <- m[
  m$project == "TCGA-BLCA",
]

if (nrow(m) == 0) {
  stop("No TCGA-BLCA files found.")
}

required <- c(
  "file_id",
  "sample_barcode",
  "sample_type"
)

miss <- setdiff(required, names(m))

if (length(miss) > 0) {
  stop(
    "Missing columns: ",
    paste(miss, collapse = ", ")
  )
}

# one manifest row per file
m <- m[
  !duplicated(m$file_id),
]

m$patient_id <-
  substr(m$sample_barcode, 1, 12)

cat("========================================\n")
cat("TCGA-BLCA JAK3 + EPHB2 EXTRACTION\n")
cat("========================================\n\n")

cat("Files:", nrow(m), "\n")

cat("\nSample types:\n")
print(table(m$sample_type))

# ------------------------------------------------------------
# Targets
# ------------------------------------------------------------

targets <- data.frame(
  gene_name = c(
    "JAK3",
    "EPHB2"
  ),

  ensembl_id = c(
    "ENSG00000105639",
    "ENSG00000133216"
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Resume support
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

  counts <- table(ok$file_id)

  done <- names(
    counts[counts >= 2]
  )
}

todo <- m[
  !(m$file_id %in% done),
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
# Helpers
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


read_star_file <- function(path) {

  con <- gzfile(
    path,
    open = "rt"
  )

  first <- readLines(
    con,
    n = 20,
    warn = FALSE
  )

  close(con)

  header_line <- which(
    grepl(
      "^gene_id\\tgene_name\\t",
      first
    )
  )[1]

  if (is.na(header_line)) {

    stop(
      "Could not locate STAR-count header."
    )
  }

  read.delim(
    gzfile(path),
    skip = header_line - 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


save_checkpoint <- function(old, new_rows) {

  if (length(new_rows) == 0) {
    return(invisible(NULL))
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

  # keep latest result for file × gene
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

  invisible(all)
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

  dest <- tempfile(
    fileext = ".tsv.gz"
  )

  downloaded <-
    download_with_retry(
      url,
      dest,
      attempts = 3
    )

  # ----------------------------------------------------------
  # Download failure
  # ----------------------------------------------------------

  if (!downloaded) {

    cat("    DOWNLOAD FAILED\n")

    for (g in seq_len(nrow(targets))) {

      new_rows[[length(new_rows) + 1]] <-
        data.frame(
          project = "TCGA-BLCA",
          file_id = row$file_id,
          sample_barcode =
            row$sample_barcode,
          patient_id =
            row$patient_id,
          sample_type =
            row$sample_type,
          gene_name =
            targets$gene_name[g],
          ensembl_id =
            targets$ensembl_id[g],
          unstranded =
            NA_real_,
          TPM =
            NA_real_,
          FPKM =
            NA_real_,
          FPKM_UQ =
            NA_real_,
          extraction_status =
            "DOWNLOAD_FAILED",
          stringsAsFactors = FALSE
        )
    }

    next
  }

  # ----------------------------------------------------------
  # Read STAR-count file
  # ----------------------------------------------------------

  dat <- tryCatch(
    read_star_file(dest),
    error = function(e) NULL
  )

  unlink(dest)

  if (is.null(dat)) {

    cat("    READ FAILED\n")

    for (g in seq_len(nrow(targets))) {

      new_rows[[length(new_rows) + 1]] <-
        data.frame(
          project = "TCGA-BLCA",
          file_id = row$file_id,
          sample_barcode =
            row$sample_barcode,
          patient_id =
            row$patient_id,
          sample_type =
            row$sample_type,
          gene_name =
            targets$gene_name[g],
          ensembl_id =
            targets$ensembl_id[g],
          unstranded =
            NA_real_,
          TPM =
            NA_real_,
          FPKM =
            NA_real_,
          FPKM_UQ =
            NA_real_,
          extraction_status =
            "READ_FAILED",
          stringsAsFactors = FALSE
        )
    }

    next
  }

  ens_base <- sub(
    "\\..*$",
    "",
    dat$gene_id
  )

  # ----------------------------------------------------------
  # Extract BOTH genes
  # ----------------------------------------------------------

  for (g in seq_len(nrow(targets))) {

    gene <-
      targets$gene_name[g]

    ensid <-
      targets$ensembl_id[g]

    hit <- which(
      dat$gene_name == gene &
      ens_base == ensid
    )

    if (length(hit) != 1) {

      cat(
        "    ",
        gene,
        " match count: ",
        length(hit),
        "\n",
        sep = ""
      )

      new_rows[[length(new_rows) + 1]] <-
        data.frame(
          project = "TCGA-BLCA",
          file_id = row$file_id,
          sample_barcode =
            row$sample_barcode,
          patient_id =
            row$patient_id,
          sample_type =
            row$sample_type,
          gene_name =
            gene,
          ensembl_id =
            ensid,
          unstranded =
            NA_real_,
          TPM =
            NA_real_,
          FPKM =
            NA_real_,
          FPKM_UQ =
            NA_real_,
          extraction_status =
            "GENE_NOT_UNIQUE",
          stringsAsFactors = FALSE
        )

      next
    }

    h <- hit[1]

    new_rows[[length(new_rows) + 1]] <-
      data.frame(
        project = "TCGA-BLCA",
        file_id = row$file_id,
        sample_barcode =
          row$sample_barcode,
        patient_id =
          row$patient_id,
        sample_type =
          row$sample_type,
        gene_name =
          gene,
        ensembl_id =
          ensid,

        unstranded =
          as.numeric(
            dat$unstranded[h]
          ),

        TPM =
          as.numeric(
            dat$tpm_unstranded[h]
          ),

        FPKM =
          as.numeric(
            dat$fpkm_unstranded[h]
          ),

        FPKM_UQ =
          as.numeric(
            dat$fpkm_uq_unstranded[h]
          ),

        extraction_status =
          "OK",

        stringsAsFactors = FALSE
      )

    cat(
      sprintf(
        "    %s TPM = %.4f\n",
        gene,
        as.numeric(
          dat$tpm_unstranded[h]
        )
      )
    )
  }

  # ----------------------------------------------------------
  # Checkpoint
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
# Final summary
# ------------------------------------------------------------

res <- read.csv(
  outfile,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("FINAL BLCA EXTRACTION SUMMARY\n")
cat("========================================\n\n")

cat(
  "Unique files:",
  length(unique(res$file_id)),
  "\n"
)

cat("\nExtraction status:\n")
print(
  table(
    res$gene_name,
    res$extraction_status,
    useNA = "ifany"
  )
)

cat("\nSample counts:\n")
print(
  table(
    res$gene_name,
    res$sample_type
  )
)

cat("\nTPM summaries:\n")

for (g in c("JAK3", "EPHB2")) {

  cat("\n---", g, "---\n")

  tmp <- res[
    res$gene_name == g &
    res$extraction_status == "OK",
  ]

  print(
    tapply(
      tmp$TPM,
      tmp$sample_type,
      summary
    )
  )
}

cat("\nSaved:\n ", outfile, "\n")

cat("\n========================================\n")
cat("14a COMPLETE\n")
cat("========================================\n")
