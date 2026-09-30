# ============================================================
# 58d — Collapse TCGA ATAC technical replicates to tissue level
#       and summarize JAK3 promoter accessibility
#
# Input:
#   46 JAK3-region peaks x technical replicates
#
# Technical replicates from the same Stanford tissue UUID
# are averaged on the normalized log2 accessibility scale.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

matrix_file <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results/JAK3_ATAC_normalized_46peaks_9cancers.rds"

mapping_file <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/data/TCGA_identifier_mapping.txt"

promoter_file <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results/JAK3_ATAC_promoter_peaks.csv"

outdir <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results"

priority_cancers <- c(
  "CHOL","KIRC","KIRP","THCA","HNSC",
  "STAD","LUAD","COAD","LIHC"
)

cat("\n========================================\n")
cat("58d JAK3 ATAC TISSUE-LEVEL COLLAPSE\n")
cat("========================================\n\n")

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------

x <- readRDS(matrix_file)

map <- read.delim(
  mapping_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

promoter <- read.csv(
  promoter_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
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

sample_cols <- setdiff(
  names(x),
  meta_cols
)

cat("JAK3 peaks:", nrow(x), "\n")
cat("Technical-replicate columns:", length(sample_cols), "\n")
cat("Identifier-mapping rows:", nrow(map), "\n\n")

# ------------------------------------------------------------
# Convert bam_prefix to matrix-column format
#
# Example:
# BRCA-UUID-X017... -> BRCA_UUID_X017...
# ------------------------------------------------------------

map$matrix_column <- gsub(
  "-",
  "_",
  map$bam_prefix
)

# Match every ATAC column to identifier table
m <- map[
  match(
    sample_cols,
    map$matrix_column
  ),
]

m$matrix_column_actual <- sample_cols

cat(
  "Mapped technical replicates:",
  sum(!is.na(m$stanfordUUID)),
  "of",
  length(sample_cols),
  "\n"
)

if (any(is.na(m$stanfordUUID))) {

  cat("\nUnmapped matrix columns:\n")

  print(
    sample_cols[
      is.na(m$stanfordUUID)
    ]
  )

  stop("Some ATAC matrix columns could not be mapped.")
}

# ------------------------------------------------------------
# Cancer type
# ------------------------------------------------------------

m$Cancer <- sub(
  "^([^_-]+)[_-].*$",
  "\\1",
  sample_cols
)

m$Cancer <- sub(
  "x$",
  "",
  m$Cancer
)

if (!all(m$Cancer %in% priority_cancers)) {
  stop("Unexpected cancer type detected.")
}

# ------------------------------------------------------------
# TCGA patient ID
# ------------------------------------------------------------

m$patient_id <- substr(
  m$aliquot_id,
  1,
  12
)

# ------------------------------------------------------------
# Audit technical replicates per tissue fragment
# ------------------------------------------------------------

rep_counts <- table(
  m$stanfordUUID
)

cat(
  "Unique tissue fragments:",
  length(rep_counts),
  "\n"
)

cat("\nTechnical replicates per tissue fragment:\n")
print(
  table(rep_counts)
)

# ------------------------------------------------------------
# Collapse technical replicates
#
# Average normalized log2 ATAC values across technical
# replicates belonging to the same tissue fragment.
# ------------------------------------------------------------

groups <- split(
  seq_along(sample_cols),
  m$stanfordUUID
)

collapsed <- x[
  ,
  meta_cols,
  drop = FALSE
]

tissue_rows <- list()

for (uuid in names(groups)) {

  idx <- groups[[uuid]]

  cols <- sample_cols[idx]

  info <- m[idx, , drop = FALSE]

  cancer <- unique(info$Cancer)

  if (length(cancer) != 1) {
    stop(
      paste(
        "Multiple cancer labels for tissue:",
        uuid
      )
    )
  }

  tissue_name <- paste0(
    cancer,
    "__",
    gsub("-", "_", uuid)
  )

  vals <- as.matrix(
    x[
      ,
      cols,
      drop = FALSE
    ]
  )

  storage.mode(vals) <- "numeric"

  collapsed[[tissue_name]] <-
    rowMeans(
      vals,
      na.rm = TRUE
    )

  tissue_rows[[
    length(tissue_rows) + 1
  ]] <- data.frame(

    tissue_column =
      tissue_name,

    Cancer =
      cancer,

    stanfordUUID =
      uuid,

    patient_id =
      unique(info$patient_id)[1],

    aliquot_id =
      unique(info$aliquot_id)[1],

    Case_UUID =
      unique(info$Case_UUID)[1],

    Case_ID =
      unique(info$Case_ID)[1],

    n_technical_replicates =
      length(idx),

    stringsAsFactors = FALSE
  )
}

tissue_map <- do.call(
  rbind,
  tissue_rows
)

# ------------------------------------------------------------
# Expected tissue counts
# ------------------------------------------------------------

cat("\n========================================\n")
cat("TISSUE-LEVEL COUNTS\n")
cat("========================================\n\n")

tissue_counts <- as.data.frame(
  table(
    factor(
      tissue_map$Cancer,
      levels = priority_cancers
    )
  )
)

names(tissue_counts) <- c(
  "Cancer",
  "Tissue_samples"
)

print(
  tissue_counts,
  row.names = FALSE
)

cat(
  "\nTotal tissue samples:",
  nrow(tissue_map),
  "\n"
)

cat(
  "Unique TCGA patients:",
  length(
    unique(
      tissue_map$patient_id
    )
  ),
  "\n"
)

# ------------------------------------------------------------
# Save collapsed matrix + tissue map
# ------------------------------------------------------------

saveRDS(
  collapsed,
  file.path(
    outdir,
    "JAK3_ATAC_46peaks_tissue_level.rds"
  )
)

write.csv(
  tissue_map,
  file.path(
    outdir,
    "JAK3_ATAC_tissue_level_mapping.csv"
  ),
  row.names = FALSE
)

# ============================================================
# JAK3 PROMOTER ACCESSIBILITY
# ============================================================

if (nrow(promoter) != 1) {
  stop(
    paste(
      "Expected exactly 1 JAK3 promoter peak; found",
      nrow(promoter)
    )
  )
}

promoter_name <- promoter$name[1]

pidx <- match(
  promoter_name,
  collapsed$name
)

if (is.na(pidx)) {
  stop("JAK3 promoter peak not found in collapsed matrix.")
}

tissue_cols <- tissue_map$tissue_column

promoter_long <- data.frame(

  tissue_column =
    tissue_cols,

  accessibility =
    as.numeric(
      collapsed[
        pidx,
        tissue_cols
      ]
    ),

  stringsAsFactors = FALSE
)

promoter_long <- merge(
  tissue_map,
  promoter_long,
  by = "tissue_column",
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------------------
# Cancer-level promoter summary
# ------------------------------------------------------------

summary_rows <- list()

for (ca in priority_cancers) {

  d <- promoter_long[
    promoter_long$Cancer == ca &
    !is.na(promoter_long$accessibility),
  ]

  summary_rows[[
    length(summary_rows) + 1
  ]] <- data.frame(

    Cancer = ca,

    Tissue_samples =
      nrow(d),

    Patients =
      length(
        unique(
          d$patient_id
        )
      ),

    Median_promoter_accessibility =
      median(
        d$accessibility
      ),

    Mean_promoter_accessibility =
      mean(
        d$accessibility
      ),

    Q1 =
      as.numeric(
        quantile(
          d$accessibility,
          0.25
        )
      ),

    Q3 =
      as.numeric(
        quantile(
          d$accessibility,
          0.75
        )
      ),

    stringsAsFactors = FALSE
  )
}

promoter_summary <- do.call(
  rbind,
  summary_rows
)

promoter_summary$Promoter_peak <-
  promoter_name

promoter_summary$Promoter_coordinates <-
  paste0(
    promoter$seqnames[1],
    ":",
    promoter$start[1],
    "-",
    promoter$end[1]
  )

promoter_summary <- promoter_summary[
  order(
    -promoter_summary$Median_promoter_accessibility
  ),
]

promoter_summary$ATAC_accessibility_rank <-
  seq_len(
    nrow(promoter_summary)
  )

# ------------------------------------------------------------
# Save promoter results
# ------------------------------------------------------------

write.csv(
  promoter_long,
  file.path(
    outdir,
    "JAK3_ATAC_promoter_tissue_values.csv"
  ),
  row.names = FALSE
)

write.csv(
  promoter_summary,
  file.path(
    outdir,
    "JAK3_ATAC_promoter_summary_9cancers.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console report
# ------------------------------------------------------------

cat("\n========================================\n")
cat("JAK3 PROMOTER ACCESSIBILITY\n")
cat("========================================\n\n")

cat(
  "Promoter peak:",
  promoter_name,
  "\n"
)

cat(
  "Coordinates:",
  promoter_summary$Promoter_coordinates[1],
  "\n\n"
)

print(
  promoter_summary[
    ,
    c(
      "ATAC_accessibility_rank",
      "Cancer",
      "Tissue_samples",
      "Patients",
      "Median_promoter_accessibility",
      "Mean_promoter_accessibility",
      "Q1",
      "Q3"
    )
  ],
  row.names = FALSE
)

cat("\nSaved tissue-level matrix and promoter summaries.\n")

cat("\n========================================\n")
cat("58d COMPLETE\n")
cat("========================================\n")
