# ============================================================
# 58f — Extract official JAK3-linked TCGA ATAC elements
#       and summarize accessibility across priority tumors
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(
  library(readxl)
)

s7_file <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/data/TCGA-ATAC_DataS7_PeakToGeneLinks_v2.xlsx"

rds_file <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/data/TCGA-ATAC_PanCan_Log2Norm_Counts.rds"

mapping_file <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/data/TCGA_identifier_mapping.txt"

outdir <-
  "TCGA_JAK3_pan_cancer/ATAC_followup/results"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

priority_cancers <- c(
  "CHOL","KIRC","KIRP","THCA","HNSC",
  "STAD","LUAD","COAD","LIHC"
)

cat("\n========================================\n")
cat("58f JAK3 OFFICIAL LINKED ATAC ELEMENTS\n")
cat("========================================\n\n")

# ------------------------------------------------------------
# Read Data S7 using actual header row
# ------------------------------------------------------------

links <- read_excel(
  s7_file,
  sheet = "All_Links_Merge",
  skip = 23
)

links <- as.data.frame(
  links,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

jak3_links <- links[
  links$Linked_Gene == "JAK3",
]

cat("Official JAK3-linked enhancer units:", nrow(jak3_links), "\n\n")

if (nrow(jak3_links) == 0) {
  stop("No JAK3 links found.")
}

print(
  jak3_links[
    ,
    c(
      "Chromosome",
      "Start",
      "End",
      "Width",
      "Enhancer_ID",
      "Enhancer_Distance",
      "Linked_Gene",
      "Linked_Gene_Start",
      "Summits",
      "Enhancer_Correlation",
      "Enhancer_FDR"
    )
  ],
  row.names = FALSE
)

write.csv(
  jak3_links,
  file.path(
    outdir,
    "JAK3_DataS7_official_links.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Load normalized pan-cancer ATAC matrix
# ------------------------------------------------------------

cat("\nLoading normalized ATAC matrix...\n")

x <- readRDS(rds_file)

meta_cols <- c(
  "seqnames",
  "start",
  "end",
  "name",
  "score",
  "annotation",
  "GC"
)

# ------------------------------------------------------------
# Match the 3 JAK3 enhancer intervals to ATAC peaks
#
# All 3 JAK3 enhancer units are 501 bp wide, so each
# corresponds to a single fixed-width TCGA ATAC peak.
# ------------------------------------------------------------

key_matrix <- paste(
  x$seqnames,
  x$start,
  x$end,
  sep = ":"
)

key_links <- paste(
  jak3_links$Chromosome,
  jak3_links$Start,
  jak3_links$End,
  sep = ":"
)

idx <- match(
  key_links,
  key_matrix
)

cat(
  "\nMatched JAK3 linked elements:",
  sum(!is.na(idx)),
  "of",
  nrow(jak3_links),
  "\n"
)

if (any(is.na(idx))) {
  cat("\nUnmatched links:\n")
  print(jak3_links[is.na(idx), ])
  stop("Not all JAK3-linked elements matched ATAC matrix.")
}

linked_matrix <- x[
  idx,
  ,
  drop = FALSE
]

linked_matrix$Enhancer_ID <-
  jak3_links$Enhancer_ID

linked_matrix$Enhancer_Distance <-
  jak3_links$Enhancer_Distance

linked_matrix$Enhancer_Correlation <-
  jak3_links$Enhancer_Correlation

linked_matrix$Enhancer_FDR <-
  jak3_links$Enhancer_FDR

# ------------------------------------------------------------
# Priority-cancer technical replicate columns
# ------------------------------------------------------------

sample_cols <- setdiff(
  names(x),
  meta_cols
)

cancer_from_col <- sub(
  "^([A-Za-z0-9]+)[_-].*$",
  "\\1",
  sample_cols
)

cancer_from_col <- sub(
  "x$",
  "",
  cancer_from_col
)

keep <- cancer_from_col %in% priority_cancers

priority_cols <- sample_cols[keep]

# ------------------------------------------------------------
# Identifier mapping
# ------------------------------------------------------------

map <- read.delim(
  mapping_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

map$matrix_column <- gsub(
  "-",
  "_",
  map$bam_prefix
)

m <- map[
  match(
    priority_cols,
    map$matrix_column
  ),
]

m$matrix_column_actual <- priority_cols

if (any(is.na(m$stanfordUUID))) {
  stop("Some priority ATAC columns failed identifier mapping.")
}

m$Cancer <- sub(
  "^([A-Za-z0-9]+)[_-].*$",
  "\\1",
  priority_cols
)

m$Cancer <- sub(
  "x$",
  "",
  m$Cancer
)

m$patient_id <- substr(
  m$aliquot_id,
  1,
  12
)

# ------------------------------------------------------------
# Collapse technical replicates to tissue fragment
# ------------------------------------------------------------

groups <- split(
  seq_along(priority_cols),
  m$stanfordUUID
)

long_rows <- list()

for (enh_i in seq_len(nrow(linked_matrix))) {

  enh_id <- jak3_links$Enhancer_ID[enh_i]

  for (uuid in names(groups)) {

    g <- groups[[uuid]]

    cols <- priority_cols[g]

    info <- m[g, , drop = FALSE]

    vals <- as.numeric(
      linked_matrix[
        enh_i,
        cols
      ]
    )

    long_rows[[
      length(long_rows) + 1
    ]] <- data.frame(

      Enhancer_ID = enh_id,

      Chromosome =
        jak3_links$Chromosome[enh_i],

      Start =
        jak3_links$Start[enh_i],

      End =
        jak3_links$End[enh_i],

      Enhancer_Distance =
        jak3_links$Enhancer_Distance[enh_i],

      Enhancer_Correlation =
        jak3_links$Enhancer_Correlation[enh_i],

      Enhancer_FDR =
        jak3_links$Enhancer_FDR[enh_i],

      Cancer =
        unique(info$Cancer)[1],

      stanfordUUID =
        uuid,

      patient_id =
        unique(info$patient_id)[1],

      n_technical_replicates =
        length(g),

      accessibility =
        mean(
          vals,
          na.rm = TRUE
        ),

      stringsAsFactors = FALSE
    )
  }
}

long <- do.call(
  rbind,
  long_rows
)

write.csv(
  long,
  file.path(
    outdir,
    "JAK3_linked_ATAC_tissue_values.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Cancer-level summary
# ------------------------------------------------------------

summary_rows <- list()

for (enh in unique(long$Enhancer_ID)) {

  for (ca in priority_cancers) {

    d <- long[
      long$Enhancer_ID == enh &
      long$Cancer == ca &
      !is.na(long$accessibility),
    ]

    if (nrow(d) == 0) {
      next
    }

    info <- d[1, ]

    summary_rows[[
      length(summary_rows) + 1
    ]] <- data.frame(

      Enhancer_ID = enh,

      Cancer = ca,

      Tissue_samples =
        nrow(d),

      Patients =
        length(
          unique(
            d$patient_id
          )
        ),

      Enhancer_Distance =
        info$Enhancer_Distance,

      Enhancer_Correlation =
        info$Enhancer_Correlation,

      Enhancer_FDR =
        info$Enhancer_FDR,

      Median_accessibility =
        median(
          d$accessibility
        ),

      Mean_accessibility =
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
}

summary_tab <- do.call(
  rbind,
  summary_rows
)

write.csv(
  summary_tab,
  file.path(
    outdir,
    "JAK3_linked_ATAC_summary_9cancers.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("OFFICIAL JAK3 LINKS\n")
cat("========================================\n\n")

print(
  unique(
    summary_tab[
      ,
      c(
        "Enhancer_ID",
        "Enhancer_Distance",
        "Enhancer_Correlation",
        "Enhancer_FDR"
      )
    ]
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("MEDIAN ACCESSIBILITY BY CANCER\n")
cat("========================================\n\n")

wide <- reshape(
  summary_tab[
    ,
    c(
      "Cancer",
      "Enhancer_ID",
      "Median_accessibility"
    )
  ],
  idvar = "Cancer",
  timevar = "Enhancer_ID",
  direction = "wide"
)

wide <- wide[
  match(
    priority_cancers,
    wide$Cancer
  ),
]

print(
  wide,
  row.names = FALSE
)

cat("\nSaved:\n")
cat(
  " ",
  file.path(
    outdir,
    "JAK3_DataS7_official_links.csv"
  ),
  "\n"
)

cat(
  " ",
  file.path(
    outdir,
    "JAK3_linked_ATAC_summary_9cancers.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("58f COMPLETE\n")
cat("========================================\n")
