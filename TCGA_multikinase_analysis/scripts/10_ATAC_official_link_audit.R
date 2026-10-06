
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
cat("MULTIKINASE TCGA ATAC LINK AUDIT\n")
cat("========================================\n")

# ============================================================
# Candidate genes
# ============================================================

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
# Integrated RNA ranking
# ============================================================

ranking <- fread(
  file.path(
    root,
    "results/integrated",
    "09_integrated_multikinase_evidence.csv"
  )
)

ranking <- ranking[
  gene %in% genes
]

# ============================================================
# Official TCGA ATAC Data S7
# ============================================================

s7_file <- paste0(
  "TCGA_JAK3_pan_cancer/ATAC_followup/data/",
  "TCGA-ATAC_DataS7_PeakToGeneLinks_v2.xlsx"
)

cat("\nReading official Data S7 peak-to-gene links...\n")

links <- read_excel(
  s7_file,
  sheet = "All_Links_Merge",
  skip = 23
)

links <- as.data.table(
  links
)

cat(
  "Total official links:",
  nrow(links),
  "\n"
)

cat(
  "Unique linked genes:",
  uniqueN(links$Linked_Gene),
  "\n"
)

# ============================================================
# Extract all candidate-gene links
# ============================================================

candidate_links <- links[
  Linked_Gene %in% genes
]

cat(
  "\nCandidate-linked elements:",
  nrow(candidate_links),
  "\n"
)

# ============================================================
# Basic interval validation
# ============================================================

required <- c(
  "Chromosome",
  "Start",
  "End",
  "Enhancer_ID",
  "Enhancer_Distance",
  "Linked_Gene",
  "Linked_Gene_Start",
  "Enhancer_Correlation",
  "Enhancer_FDR"
)

missing_cols <- setdiff(
  required,
  names(candidate_links)
)

if (length(missing_cols) > 0) {

  stop(
    paste(
      "Missing expected Data S7 columns:",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}

# ============================================================
# Audit per gene
# ============================================================

audit <- rbindlist(
  lapply(
    genes,
    function(g) {

      x <- candidate_links[
        Linked_Gene == g
      ]

      if (nrow(x) == 0) {

        return(
          data.table(
            gene = g,

            Official_linked_elements = 0L,

            Unique_enhancer_IDs = 0L,

            Median_abs_distance_bp =
              NA_real_,

            Min_abs_distance_bp =
              NA_real_,

            Max_abs_distance_bp =
              NA_real_,

            Median_link_correlation =
              NA_real_,

            Positive_link_correlations =
              0L,

            Negative_link_correlations =
              0L,

            Links_FDR_lt_005 =
              0L,

            Has_official_ATAC_links =
              FALSE
          )
        )
      }

      distance <-
        suppressWarnings(
          as.numeric(
            x$Enhancer_Distance
          )
        )

      correlation <-
        suppressWarnings(
          as.numeric(
            x$Enhancer_Correlation
          )
        )

      fdr <-
        suppressWarnings(
          as.numeric(
            x$Enhancer_FDR
          )
        )

      data.table(
        gene = g,

        Official_linked_elements =
          nrow(x),

        Unique_enhancer_IDs =
          uniqueN(x$Enhancer_ID),

        Median_abs_distance_bp =
          median(
            abs(distance),
            na.rm = TRUE
          ),

        Min_abs_distance_bp =
          min(
            abs(distance),
            na.rm = TRUE
          ),

        Max_abs_distance_bp =
          max(
            abs(distance),
            na.rm = TRUE
          ),

        Median_link_correlation =
          median(
            correlation,
            na.rm = TRUE
          ),

        Positive_link_correlations =
          sum(
            correlation > 0,
            na.rm = TRUE
          ),

        Negative_link_correlations =
          sum(
            correlation < 0,
            na.rm = TRUE
          ),

        Links_FDR_lt_005 =
          sum(
            fdr < 0.05,
            na.rm = TRUE
          ),

        Has_official_ATAC_links =
          TRUE
      )
    }
  ),
  fill = TRUE
)

# Fix Inf values if any gene had only missing distance values
for (
  col in c(
    "Median_abs_distance_bp",
    "Min_abs_distance_bp",
    "Max_abs_distance_bp",
    "Median_link_correlation"
  )
) {

  audit[
    is.infinite(get(col)),
    (col) := NA_real_
  ]
}

# ============================================================
# Add RNA-ranking context
# ============================================================

context_cols <- intersect(
  c(
    "gene",
    "Memory_class",
    "RNA_evidence_index",
    "JAK3_benchmark_domains",
    "RNA_ATAC_priority",
    "Scope_note"
  ),
  names(ranking)
)

audit <- merge(
  audit,
  ranking[
    ,
    ..context_cols
  ],
  by = "gene",
  all.x = TRUE
)

# ============================================================
# Order results
# ============================================================

audit[, Reference :=
  gene == "JAK3"
]

setorder(
  audit,
  -Has_official_ATAC_links,
  -Official_linked_elements,
  -RNA_evidence_index,
  gene
)

# ============================================================
# Save all official candidate links
# ============================================================

fwrite(
  candidate_links,
  file.path(
    outdir,
    "10_all_candidate_DataS7_official_links.csv"
  )
)

# One file per gene
for (g in genes) {

  x <- candidate_links[
    Linked_Gene == g
  ]

  fwrite(
    x,
    file.path(
      outdir,
      paste0(
        "10_",
        g,
        "_DataS7_official_links.csv"
      )
    )
  )
}

# ============================================================
# Save audit
# ============================================================

fwrite(
  audit,
  file.path(
    outdir,
    "10_multikinase_ATAC_link_availability.csv"
  )
)

# ============================================================
# Console report
# ============================================================

cat("\n===== OFFICIAL ATAC LINK AVAILABILITY =====\n")

print(
  audit[
    ,
    .(
      gene,
      Memory_class,
      RNA_evidence_index,
      JAK3_benchmark_domains,
      RNA_ATAC_priority,
      Official_linked_elements,
      Unique_enhancer_IDs,
      Links_FDR_lt_005,
      Median_link_correlation,
      Median_abs_distance_bp,
      Has_official_ATAC_links
    )
  ]
)

cat("\n===== GENES WITHOUT OFFICIAL LINKS =====\n")

print(
  audit[
    Has_official_ATAC_links == FALSE,
    gene
  ]
)

cat("\n===== FJX1 SCOPE NOTE =====\n")

if ("FJX1" %in% audit$gene) {

  print(
    audit[
      gene == "FJX1",
      .(
        gene,
        Memory_class,
        Scope_note,
        Official_linked_elements
      )
    ]
  )
}

cat("\n10 COMPLETE\n")

