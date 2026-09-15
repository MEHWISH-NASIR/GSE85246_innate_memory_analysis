# ============================================================
# 13e — PAN-CANCER TCGA JAK3 EXTRACTION
#
# Extract JAK3 from all TCGA projects having:
#   Primary Tumor + Solid Tissue Normal
#
# One GDC file is downloaded temporarily, JAK3 extracted,
# then deleted. Script is resumable.
# ============================================================

rm(list = ls())

options(timeout = 300)

meta_file <-
  "results/13_TCGA_JAK3/13a_TCGA_STARCounts_metadata.csv"

audit_file <-
  "results/13_TCGA_JAK3/13a_TCGA_STARCounts_sample_audit.csv"

outfile <-
  "results/13_TCGA_JAK3/13e_pancancer_JAK3_expression.csv"

if (!file.exists(meta_file) || !file.exists(audit_file)) {
  stop("Missing Step 13a files.")
}

meta <- read.csv(
  meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

audit <- read.csv(
  audit_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

eligible <- audit$Project[
  audit$Tumor_Normal_available
]

cat("\n========================================\n")
cat("13e — PAN-CANCER JAK3 EXTRACTION\n")
cat("========================================\n\n")

cat("Eligible TCGA projects:", length(eligible), "\n")
print(eligible)

x <- meta[
  meta$TCGA_project %in% eligible &
  meta$sample_type %in%
    c("Primary Tumor", "Solid Tissue Normal"),
]

x <- x[
  !duplicated(x$file_id),
]

cat("\nTotal files selected:", nrow(x), "\n")

cat("\nFiles by project:\n")
print(
  sort(
    table(x$TCGA_project),
    decreasing = TRUE
  )
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

  done <- unique(old$file_id)

  cat(
    "\nExisting records:",
    length(done),
    "\n"
  )

} else {

  old <- NULL
  done <- character(0)
}

todo <- x[
  !(x$file_id %in% done),
]

cat(
  "Remaining:",
  nrow(todo),
  "\n\n"
)

if (nrow(todo) == 0) {
  cat("All files already complete.\n")
  quit(save = "no")
}


# ------------------------------------------------------------
# Download helper
# ------------------------------------------------------------

download_retry <- function(
    url,
    dest,
    attempts = 4
) {

  for (a in seq_len(attempts)) {

    ok <- tryCatch({

      suppressWarnings(
        download.file(
          url,
          destfile = dest,
          mode = "wb",
          quiet = TRUE
        )
      )

      file.exists(dest) &&
        file.info(dest)$size > 1000

    }, error = function(e) FALSE)

    if (isTRUE(ok)) {
      return(TRUE)
    }

    if (file.exists(dest)) {
      unlink(dest)
    }

    cat(
      "    retry",
      a,
      "/",
      attempts,
      "\n"
    )

    Sys.sleep(3)
  }

  FALSE
}


# ------------------------------------------------------------
# Extraction
# ------------------------------------------------------------

new <- list()

for (i in seq_len(nrow(todo))) {

  r <- todo[i, ]

  cat(
    sprintf(
      "[%05d/%05d] %s | %s | %s\n",
      i,
      nrow(todo),
      r$TCGA_project,
      r$sample_type,
      r$cases
    )
  )

  tmp <- tempfile(
    fileext = ".tsv"
  )

  url <- paste0(
    "https://api.gdc.cancer.gov/data/",
    r$file_id
  )

  ok <- download_retry(
    url,
    tmp
  )

  if (!ok) {

    result <- data.frame(
      project = r$TCGA_project,
      sample_barcode = r$cases,
      patient_id = substr(r$cases, 1, 12),
      sample_type = r$sample_type,
      file_id = r$file_id,
      gene_id = NA_character_,
      gene_name = "JAK3",
      unstranded_count = NA_real_,
      TPM = NA_real_,
      FPKM = NA_real_,
      FPKM_UQ = NA_real_,
      extraction_status = "DOWNLOAD_FAILED",
      stringsAsFactors = FALSE
    )

  } else {

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

      result <- data.frame(
        project = r$TCGA_project,
        sample_barcode = r$cases,
        patient_id = substr(r$cases, 1, 12),
        sample_type = r$sample_type,
        file_id = r$file_id,
        gene_id = NA_character_,
        gene_name = "JAK3",
        unstranded_count = NA_real_,
        TPM = NA_real_,
        FPKM = NA_real_,
        FPKM_UQ = NA_real_,
        extraction_status = "READ_FAILED",
        stringsAsFactors = FALSE
      )

    } else {

      ens <- sub(
        "\\..*$",
        "",
        dat$gene_id
      )

      hit <- which(
        dat$gene_name == "JAK3" &
        ens == "ENSG00000105639"
      )

      if (length(hit) == 1) {

        j <- dat[hit, ]

        result <- data.frame(
          project = r$TCGA_project,
          sample_barcode = r$cases,
          patient_id = substr(r$cases, 1, 12),
          sample_type = r$sample_type,
          file_id = r$file_id,
          gene_id = j$gene_id,
          gene_name = j$gene_name,
          unstranded_count =
            as.numeric(j$unstranded),
          TPM =
            as.numeric(j$tpm_unstranded),
          FPKM =
            as.numeric(j$fpkm_unstranded),
          FPKM_UQ =
            as.numeric(j$fpkm_uq_unstranded),
          extraction_status = "OK",
          stringsAsFactors = FALSE
        )

      } else {

        result <- data.frame(
          project = r$TCGA_project,
          sample_barcode = r$cases,
          patient_id = substr(r$cases, 1, 12),
          sample_type = r$sample_type,
          file_id = r$file_id,
          gene_id = NA_character_,
          gene_name = "JAK3",
          unstranded_count = NA_real_,
          TPM = NA_real_,
          FPKM = NA_real_,
          FPKM_UQ = NA_real_,
          extraction_status = "JAK3_NOT_UNIQUE",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  new[[length(new) + 1]] <- result

  if (result$extraction_status == "OK") {
    cat(
      sprintf(
        "    TPM = %.4f\n",
        result$TPM
      )
    )
  }

  # ----------------------------------------------------------
  # checkpoint
  # ----------------------------------------------------------

  if (
    i %% 25 == 0 ||
    i == nrow(todo)
  ) {

    ndf <- do.call(
      rbind,
      new
    )

    if (is.null(old)) {
      combined <- ndf
    } else {
      combined <- rbind(
        old,
        ndf
      )
    }

    combined <- combined[
      !duplicated(
        combined$file_id,
        fromLast = TRUE
      ),
    ]

    write.csv(
      combined,
      outfile,
      row.names = FALSE
    )

    cat(
      "    CHECKPOINT:",
      nrow(combined),
      "saved\n"
    )
  }

  Sys.sleep(0.1)
}


# ============================================================
# FINAL SUMMARY
# ============================================================

final <- read.csv(
  outfile,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("FINAL PAN-CANCER EXTRACTION SUMMARY\n")
cat("========================================\n\n")

cat("Records:", nrow(final), "\n")

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

cat("Successful samples by project/type:\n\n")

print(
  with(
    final[
      final$extraction_status == "OK",
    ],
    table(
      project,
      sample_type
    )
  )
)

cat("\n========================================\n")
cat("13e COMPLETE\n")
cat("========================================\n")
