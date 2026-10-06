
library(data.table)
library(readxl)

setDTthreads(1)

root <- "TCGA_multikinase_analysis"

outdir <- file.path(
  root,
  "results/ATAC"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n========================================\n")
cat("MULTIKINASE LINKED ATAC EXTRACTION\n")
cat("SUMMIT-AWARE VERSION\n")
cat("========================================\n")

genes <- c(
  "JAK3",
  "MET",
  "HCK",
  "LYN",
  "DDR1",
  "IRAK2",
  "FJX1",
  "EPHB1",
  "RPS6KA2",
  "EPHB2"
)

# ============================================================
# Input files
# ============================================================

s7_file <- paste0(
  "TCGA_JAK3_pan_cancer/ATAC_followup/data/",
  "TCGA-ATAC_DataS7_PeakToGeneLinks_v2.xlsx"
)

rds_file <- paste0(
  "TCGA_JAK3_pan_cancer/ATAC_followup/data/",
  "TCGA-ATAC_PanCan_Log2Norm_Counts.rds"
)

mapping_file <- paste0(
  "TCGA_JAK3_pan_cancer/ATAC_followup/data/",
  "TCGA_identifier_mapping.txt"
)

# ============================================================
# 1. Read official Data S7 links
# ============================================================

cat("\n1. Reading official Data S7 links...\n")

links <- read_excel(
  s7_file,
  sheet = "All_Links_Merge",
  skip = 23
)

links <- as.data.table(
  links
)

candidate_links <- links[
  Linked_Gene %in% genes
]

cat(
  "Official candidate enhancer units:",
  nrow(candidate_links),
  "\n"
)

if (nrow(candidate_links) != 51) {
  stop("Expected 51 official candidate enhancer units.")
}

# ============================================================
# 2. Parse official summits
# ============================================================

cat("\n2. Parsing constituent ATAC summits...\n")

summit_rows <- list()

counter <- 1L

for (i in seq_len(nrow(candidate_links))) {

  summit_string <-
    as.character(
      candidate_links$Summits[i]
    )

  summits <- strsplit(
    summit_string,
    ";",
    fixed = TRUE
  )[[1]]

  summits <- suppressWarnings(
    as.numeric(
      trimws(summits)
    )
  )

  summits <- summits[
    !is.na(summits)
  ]

  if (length(summits) == 0) {
    stop(
      paste(
        "No valid summit found for",
        candidate_links$Enhancer_ID[i]
      )
    )
  }

  for (s in summits) {

    summit_rows[[counter]] <-
      data.table(

        gene =
          candidate_links$Linked_Gene[i],

        Enhancer_ID =
          candidate_links$Enhancer_ID[i],

        Chromosome =
          candidate_links$Chromosome[i],

        Enhancer_Start =
          candidate_links$Start[i],

        Enhancer_End =
          candidate_links$End[i],

        Enhancer_Width =
          candidate_links$Width[i],

        Summit =
          s,

        Enhancer_Distance =
          suppressWarnings(
            as.numeric(
              candidate_links$Enhancer_Distance[i]
            )
          ),

        Enhancer_Correlation =
          suppressWarnings(
            as.numeric(
              candidate_links$Enhancer_Correlation[i]
            )
          ),

        Enhancer_FDR =
          suppressWarnings(
            as.numeric(
              candidate_links$Enhancer_FDR[i]
            )
          )
      )

    counter <- counter + 1L
  }
}

summit_map <- rbindlist(
  summit_rows,
  fill = TRUE
)

cat(
  "Enhancer units:",
  uniqueN(summit_map$Enhancer_ID),
  "\n"
)

cat(
  "Constituent summit peaks expected:",
  nrow(summit_map),
  "\n"
)

cat("\nConstituent peaks per enhancer:\n")

print(
  summit_map[
    ,
    .N,
    by = .(
      gene,
      Enhancer_ID
    )
  ][
    N > 1
  ]
)

# ============================================================
# 3. Load normalized ATAC matrix
# ============================================================

cat("\n3. Loading normalized ATAC matrix...\n")

x <- readRDS(
  rds_file
)

cat(
  "ATAC matrix dimensions:",
  nrow(x),
  "x",
  ncol(x),
  "\n"
)

meta_cols <- c(
  "seqnames",
  "start",
  "end",
  "name",
  "score",
  "annotation",
  "GC"
)

missing_meta <- setdiff(
  meta_cols,
  names(x)
)

if (length(missing_meta) > 0) {

  stop(
    paste(
      "Missing ATAC matrix metadata:",
      paste(
        missing_meta,
        collapse = ", "
      )
    )
  )
}

# ============================================================
# 4. Map each official summit to matrix ATAC peak
# ============================================================

cat("\n4. Mapping official summits to matrix peaks...\n")

summit_map[, matrix_row := NA_integer_]
summit_map[, matrix_peak_name := NA_character_]
summit_map[, matrix_peak_start := NA_real_]
summit_map[, matrix_peak_end := NA_real_]

for (i in seq_len(nrow(summit_map))) {

  chr_i <- summit_map$Chromosome[i]
  summit_i <- summit_map$Summit[i]

  hits <- which(
    x$seqnames == chr_i &
    x$start <= summit_i &
    x$end >= summit_i
  )

  if (length(hits) == 0) {

    cat(
      "\nNO PEAK FOR SUMMIT:",
      summit_map$gene[i],
      summit_map$Enhancer_ID[i],
      chr_i,
      summit_i,
      "\n"
    )

    next
  }

  # If overlapping intervals exist, select the peak
  # whose midpoint is closest to the official summit.
  if (length(hits) > 1) {

    midpoint <-
      (
        as.numeric(x$start[hits]) +
        as.numeric(x$end[hits])
      ) / 2

    hits <- hits[
      which.min(
        abs(
          midpoint - summit_i
        )
      )
    ]
  }

  hit <- hits[1]

  summit_map$matrix_row[i] <-
    hit

  summit_map$matrix_peak_name[i] <-
    as.character(
      x$name[hit]
    )

  summit_map$matrix_peak_start[i] <-
    as.numeric(
      x$start[hit]
    )

  summit_map$matrix_peak_end[i] <-
    as.numeric(
      x$end[hit]
    )
}

cat(
  "Mapped summits:",
  sum(!is.na(summit_map$matrix_row)),
  "/",
  nrow(summit_map),
  "\n"
)

if (any(is.na(summit_map$matrix_row))) {

  cat("\nUnmapped summit records:\n")

  print(
    summit_map[
      is.na(matrix_row),
      .(
        gene,
        Enhancer_ID,
        Chromosome,
        Summit
      )
    ]
  )

  stop(
    "Not all official Data S7 summits mapped to ATAC matrix peaks."
  )
}

# Ensure each enhancer retained
if (
  uniqueN(summit_map$Enhancer_ID) !=
  nrow(candidate_links)
) {

  stop(
    "Not all official enhancer units retained."
  )
}

# ============================================================
# 5. Save summit-to-peak mapping audit
# ============================================================

fwrite(
  summit_map,
  file.path(
    outdir,
    "11_DataS7_summit_to_ATAC_peak_mapping.csv"
  )
)

mapping_summary <- summit_map[
  ,
  .(
    Constituent_peaks = .N,

    Matrix_peaks =
      uniqueN(
        matrix_peak_name
      ),

    All_summits_mapped =
      all(
        !is.na(matrix_row)
      )
  ),
  by = .(
    gene,
    Enhancer_ID
  )
]

fwrite(
  mapping_summary,
  file.path(
    outdir,
    "11_enhancer_constituent_peak_summary.csv"
  )
)

# ============================================================
# 6. ATAC sample columns
# ============================================================

sample_cols <- setdiff(
  names(x),
  meta_cols
)

cat(
  "\nATAC technical-replicate columns:",
  length(sample_cols),
  "\n"
)

# ============================================================
# 7. Identifier mapping
# ============================================================

cat("\n5. Loading ATAC identifier mapping...\n")

map <- fread(
  mapping_file
)

map[, matrix_column :=
  gsub(
    "-",
    "_",
    bam_prefix
  )
]

sample_map <- map[
  match(
    sample_cols,
    matrix_column
  )
]

sample_map[, matrix_column_actual :=
  sample_cols
]

unmapped_n <- sum(
  is.na(
    sample_map$stanfordUUID
  )
)

cat(
  "Unmapped ATAC sample columns:",
  unmapped_n,
  "\n"
)

if (unmapped_n > 0) {

  stop(
    "ATAC sample identifier mapping incomplete."
  )
}

sample_map[, Cancer :=
  sub(
    "^([A-Za-z0-9]+)[_-].*$",
    "\\1",
    matrix_column_actual
  )
]

sample_map[, Cancer :=
  sub(
    "x$",
    "",
    Cancer
  )
]

sample_map[, patient_id :=
  substr(
    aliquot_id,
    1,
    12
  )
]

cat(
  "ATAC cancers:",
  uniqueN(
    sample_map$Cancer
  ),
  "\n"
)

fwrite(
  sample_map,
  file.path(
    outdir,
    "11_ATAC_sample_identifier_mapping.csv"
  )
)

# ============================================================
# 8. Technical replicate groups
# ============================================================

groups <- split(
  seq_along(sample_cols),
  sample_map$stanfordUUID
)

cat(
  "Unique ATAC tissue fragments:",
  length(groups),
  "\n"
)

# ============================================================
# 9. Extract enhancer accessibility
#
# For multi-peak Data S7 enhancer units:
# first average constituent ATAC peaks within each technical
# replicate, then average technical replicates belonging to
# the same Stanford tissue fragment.
#
# Matrix values are already log2-normalized.
# ============================================================

cat(
  "\n6. Extracting official enhancer-unit accessibility...\n"
)

enhancers <- unique(
  summit_map[
    ,
    .(
      gene,
      Enhancer_ID
    )
  ]
)

long_rows <- vector(
  "list",
  nrow(enhancers) *
  length(groups)
)

counter <- 1L

for (e in seq_len(nrow(enhancers))) {

  gname <-
    enhancers$gene[e]

  enh <-
    enhancers$Enhancer_ID[e]

  sm <- summit_map[
    gene == gname &
    Enhancer_ID == enh
  ]

  peak_rows <-
    unique(
      sm$matrix_row
    )

  for (uuid in names(groups)) {

    technical_idx <-
      groups[[uuid]]

    cols <-
      sample_cols[
        technical_idx
      ]

    info <-
      sample_map[
        technical_idx
      ]

    # Matrix values:
    # rows = constituent peaks
    # cols = technical replicates
    vals <- as.matrix(
      x[
        peak_rows,
        cols,
        drop = FALSE
      ]
    )

    storage.mode(vals) <- "numeric"

    # Average constituent peaks for each technical replicate
    replicate_accessibility <-
      colMeans(
        vals,
        na.rm = TRUE
      )

    # Collapse technical replicates to tissue fragment
    enhancer_accessibility <-
      mean(
        replicate_accessibility,
        na.rm = TRUE
      )

    long_rows[[counter]] <-
      data.table(

        gene =
          gname,

        Enhancer_ID =
          enh,

        Chromosome =
          sm$Chromosome[1],

        Enhancer_Start =
          sm$Enhancer_Start[1],

        Enhancer_End =
          sm$Enhancer_End[1],

        Enhancer_Distance =
          sm$Enhancer_Distance[1],

        Enhancer_Correlation =
          sm$Enhancer_Correlation[1],

        Enhancer_FDR =
          sm$Enhancer_FDR[1],

        Constituent_peak_count =
          length(
            peak_rows
          ),

        Cancer =
          unique(
            info$Cancer
          )[1],

        stanfordUUID =
          uuid,

        patient_id =
          unique(
            info$patient_id
          )[1],

        n_technical_replicates =
          length(
            technical_idx
          ),

        accessibility =
          enhancer_accessibility
      )

    counter <- counter + 1L
  }
}

long_rows <- long_rows[
  seq_len(
    counter - 1L
  )
]

long <- rbindlist(
  long_rows,
  fill = TRUE
)

if (
  uniqueN(long$Enhancer_ID) != 51
) {
  stop(
    "Final tissue-level table does not contain all 51 enhancer units."
  )
}

cat(
  "Tissue-level enhancer observations:",
  nrow(long),
  "\n"
)

cat(
  "Genes retained:",
  uniqueN(long$gene),
  "\n"
)

cat(
  "Official enhancer units retained:",
  uniqueN(long$Enhancer_ID),
  "\n"
)

cat(
  "Unique ATAC patients:",
  uniqueN(long$patient_id),
  "\n"
)

fwrite(
  long,
  file.path(
    outdir,
    "11_multikinase_linked_ATAC_tissue_values.csv"
  )
)

# ============================================================
# 10. Collapse tissue fragments to patient level
# ============================================================

patient_atac <- long[
  is.finite(accessibility),
  .(
    accessibility =
      mean(
        accessibility,
        na.rm = TRUE
      ),

    tissue_fragments =
      uniqueN(
        stanfordUUID
      ),

    Constituent_peak_count =
      first(
        Constituent_peak_count
      ),

    Enhancer_Distance =
      first(
        Enhancer_Distance
      ),

    Enhancer_Correlation =
      first(
        Enhancer_Correlation
      ),

    Enhancer_FDR =
      first(
        Enhancer_FDR
      )
  ),
  by = .(
    gene,
    Cancer,
    Enhancer_ID,
    patient_id
  )
]

fwrite(
  patient_atac,
  file.path(
    outdir,
    "11_multikinase_linked_ATAC_patient_values.csv"
  )
)

# ============================================================
# 11. Availability by enhancer/cancer
# ============================================================

availability <- patient_atac[
  ,
  .(
    ATAC_patients =
      uniqueN(
        patient_id
      ),

    Median_accessibility =
      median(
        accessibility,
        na.rm = TRUE
      ),

    Mean_accessibility =
      mean(
        accessibility,
        na.rm = TRUE
      )
  ),
  by = .(
    gene,
    Cancer,
    Enhancer_ID,
    Constituent_peak_count
  )
]

fwrite(
  availability,
  file.path(
    outdir,
    "11_multikinase_ATAC_availability.csv"
  )
)

# ============================================================
# 12. Gene-level availability
# ============================================================

gene_summary <- patient_atac[
  ,
  .(
    Official_elements =
      uniqueN(
        Enhancer_ID
      ),

    Constituent_matrix_peaks =
      uniqueN(
        summit_map[
          gene == .BY$gene,
          matrix_peak_name
        ]
      ),

    ATAC_cancers =
      uniqueN(
        Cancer
      ),

    ATAC_patients =
      uniqueN(
        patient_id
      ),

    Patient_element_observations =
      .N
  ),
  by = gene
]

setorder(
  gene_summary,
  -Official_elements,
  gene
)

fwrite(
  gene_summary,
  file.path(
    outdir,
    "11_ATAC_gene_availability_summary.csv"
  )
)

# ============================================================
# 13. JAK3 validation
# ============================================================

cat("\n===== JAK3 EXTRACTION CHECK =====\n")

print(
  gene_summary[
    gene == "JAK3"
  ]
)

cat("\nJAK3 constituent peak mapping:\n")

print(
  mapping_summary[
    gene == "JAK3"
  ]
)

cat("\n===== ALL-GENE ATAC AVAILABILITY =====\n")

print(
  gene_summary
)

cat("\n===== MULTI-PEAK ENHANCER UNITS =====\n")

print(
  mapping_summary[
    Constituent_peaks > 1
  ]
)

cat("\n11 COMPLETE\n")

