library(dplyr)

# ============================================================
# 53a — PAN-CANCER IMMUNE MARKER EXTRACTION
#
# One-time extraction for reuse across all candidate kinases.
#
# Markers:
# PTPRC, CD3D, CD3E, MS4A1, CD79A,
# NKG7, LST1, TYROBP, FCER1G, CD68
#
# Expression:
# GDC STAR tpm_unstranded
# ============================================================

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

metadata_file <-
  "TCGA_JAK3_pan_cancer/results/expression/JAK3_GDC_metadata.csv"

outdir <-
  "TCGA_JAK3_pan_cancer/results/immune_markers"

outfile <-
  file.path(
    outdir,
    "TCGA_pan_cancer_10_immune_markers_TPM.csv"
  )

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# METADATA
# ============================================================

meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("PAN-CANCER IMMUNE MARKER EXTRACTION\n")
cat("========================================\n")

cat("\nMetadata rows:", nrow(meta), "\n")
cat("Target markers:", length(markers), "\n")

# ============================================================
# FIND EXACT STAR FILES
# ============================================================

star_files <- list.files(
  "GDCdata",
  pattern =
    "rna_seq\\.augmented_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Exact STAR files found:",
    length(star_files), "\n")

if (length(star_files) != 11505) {

  stop(
    paste(
      "Expected 11,505 STAR files but found",
      length(star_files)
    )
  )
}

# Map basename -> full path
star_index <- data.frame(
  file_name = basename(star_files),
  path = star_files,
  stringsAsFactors = FALSE
)

if (anyDuplicated(star_index$file_name)) {
  stop("Duplicated STAR filenames detected.")
}

meta <- meta %>%
  left_join(
    star_index,
    by = "file_name"
  )

cat(
  "Metadata rows with local STAR path:",
  sum(!is.na(meta$path)),
  "/",
  nrow(meta),
  "\n"
)

if (any(is.na(meta$path))) {
  stop("Some metadata rows have no local STAR file.")
}

# ============================================================
# RESUME SUPPORT
# ============================================================

if (file.exists(outfile)) {

  old <- read.csv(
    outfile,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  completed <- old %>%
    filter(extraction_status == "OK") %>%
    count(file_name) %>%
    filter(n == length(markers)) %>%
    pull(file_name)

} else {

  old <- data.frame()
  completed <- character(0)
}

todo <- meta %>%
  filter(!(file_name %in% completed))

cat(
  "Previously completed files:",
  length(completed),
  "\n"
)

cat(
  "Files remaining:",
  nrow(todo),
  "\n"
)

if (nrow(todo) == 0) {
  cat("\nNothing left to extract.\n")
  quit(save = "no")
}

# ============================================================
# CHECKPOINT FUNCTION
# ============================================================

save_checkpoint <- function(old, new_rows) {

  if (length(new_rows) == 0) {
    return(old)
  }

  new_df <- bind_rows(new_rows)

  if (nrow(old) > 0) {
    all <- bind_rows(old, new_df)
  } else {
    all <- new_df
  }

  key <- paste(
    all$file_name,
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

# ============================================================
# EXTRACTION
# ============================================================

new_rows <- list()

for (i in seq_len(nrow(todo))) {

  m <- todo[i, ]

  if (
    i == 1 ||
    i %% 50 == 0 ||
    i == nrow(todo)
  ) {

    cat(
      sprintf(
        "[%05d/%05d] %s\n",
        i,
        nrow(todo),
        m$file_name
      )
    )
  }

  dat <- tryCatch(

    read.delim(
      m$path,
      skip = 1,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),

    error = function(e) NULL
  )

  # ----------------------------------------------------------
  # FILE READ FAILURE
  # ----------------------------------------------------------

  if (is.null(dat)) {

    for (g in markers) {

      new_rows[[
        length(new_rows) + 1
      ]] <- data.frame(

        project = m$project,
        file_id = m$file_id,
        file_name = m$file_name,

        cases.submitter_id =
          m$cases.submitter_id,

        sample.submitter_id =
          m$sample.submitter_id,

        sample_type =
          m$sample_type,

        gene_name = g,
        gene_id = NA_character_,
        TPM = NA_real_,
        extraction_status =
          "FILE_READ_ERROR",

        stringsAsFactors = FALSE
      )
    }

    next
  }

  required_cols <- c(
    "gene_id",
    "gene_name",
    "tpm_unstranded"
  )

  if (!all(required_cols %in% names(dat))) {

    stop(
      paste(
        "Required STAR columns missing in",
        m$file_name
      )
    )
  }

  # ----------------------------------------------------------
  # MARKER EXTRACTION
  # ----------------------------------------------------------

  for (g in markers) {

    h <- which(
      dat$gene_name == g
    )

    if (length(h) == 1) {

      new_rows[[
        length(new_rows) + 1
      ]] <- data.frame(

        project = m$project,
        file_id = m$file_id,
        file_name = m$file_name,

        cases.submitter_id =
          m$cases.submitter_id,

        sample.submitter_id =
          m$sample.submitter_id,

        sample_type =
          m$sample_type,

        gene_name = g,

        gene_id =
          dat$gene_id[h],

        TPM =
          as.numeric(
            dat$tpm_unstranded[h]
          ),

        extraction_status =
          "OK",

        stringsAsFactors = FALSE
      )

    } else {

      new_rows[[
        length(new_rows) + 1
      ]] <- data.frame(

        project = m$project,
        file_id = m$file_id,
        file_name = m$file_name,

        cases.submitter_id =
          m$cases.submitter_id,

        sample.submitter_id =
          m$sample.submitter_id,

        sample_type =
          m$sample_type,

        gene_name = g,
        gene_id = NA_character_,
        TPM = NA_real_,

        extraction_status =
          ifelse(
            length(h) == 0,
            "GENE_NOT_FOUND",
            "MULTIPLE_GENE_MATCHES"
          ),

        stringsAsFactors = FALSE
      )
    }
  }

  # ----------------------------------------------------------
  # CHECKPOINT EVERY 100 FILES
  # ----------------------------------------------------------

  if (
    i %% 100 == 0 ||
    i == nrow(todo)
  ) {

    old <- save_checkpoint(
      old,
      new_rows
    )

    new_rows <- list()

    cat(
      "  CHECKPOINT:",
      length(unique(old$file_name)),
      "files saved\n"
    )
  }
}

# ============================================================
# FINAL AUDIT
# ============================================================

res <- read.csv(
  outfile,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========================================\n")
cat("FINAL PAN-CANCER IMMUNE-MARKER AUDIT\n")
cat("========================================\n")

cat(
  "\nRows:",
  nrow(res),
  "\n"
)

cat(
  "Unique files:",
  length(unique(res$file_name)),
  "\n"
)

cat(
  "Projects:",
  length(unique(res$project)),
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

cat("\nExpected rows:",
    11505 * length(markers),
    "\n")

cat(
  "Observed rows:",
  nrow(res),
  "\n"
)

cat(
  "Missing TPM:",
  sum(is.na(res$TPM)),
  "\n"
)

cat(
  "\nSaved:\n",
  outfile,
  "\n"
)

cat("\n53a COMPLETE\n")
