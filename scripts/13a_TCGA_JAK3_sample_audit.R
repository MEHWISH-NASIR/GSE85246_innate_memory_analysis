# ============================================================
# 13a — TCGA JAK3 PAN-CANCER RNA-SEQ SAMPLE AUDIT
#
# Purpose:
#   Before downloading JAK3 expression, identify all TCGA
#   projects containing GDC STAR-Counts RNA-seq data and
#   determine which sample types are represented.
#
# IMPORTANT:
#   No expression files are downloaded in this step.
# ============================================================

rm(list = ls())

if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  stop("TCGAbiolinks is not installed.")
}

library(TCGAbiolinks)

outdir <- "results/13_TCGA_JAK3"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n========================================\n")
cat("13a — TCGA JAK3 SAMPLE AUDIT\n")
cat("========================================\n\n")


# ============================================================
# 1. CHECK GDC STATUS
# ============================================================

cat("Checking GDC server...\n")

gdc_info <- tryCatch(
  TCGAbiolinks::getGDCInfo(),
  error = function(e) NULL
)

if (is.null(gdc_info)) {
  stop("Could not connect to the GDC API.")
}

cat("GDC connection: OK\n\n")


# ============================================================
# 2. GET ALL TCGA PROJECTS
# ============================================================

projects_all <- TCGAbiolinks::getGDCprojects()

tcga_projects <- projects_all[
  grepl(
    "^TCGA-",
    projects_all$project_id
  ),
  ,
  drop = FALSE
]

tcga_projects <- tcga_projects[
  order(tcga_projects$project_id),
]

cat(
  "TCGA projects identified:",
  nrow(tcga_projects),
  "\n\n"
)

print(
  tcga_projects$project_id
)


# ============================================================
# 3. QUERY STAR-COUNTS METADATA FOR EACH TCGA PROJECT
# ============================================================

detail_list <- list()
summary_list <- list()

for (i in seq_len(nrow(tcga_projects))) {

  project <- tcga_projects$project_id[i]

  cat("\n========================================\n")
  cat(
    sprintf(
      "[%02d/%02d] %s\n",
      i,
      nrow(tcga_projects),
      project
    )
  )
  cat("========================================\n")

  query <- tryCatch(

    TCGAbiolinks::GDCquery(
      project = project,
      data.category =
        "Transcriptome Profiling",
      data.type =
        "Gene Expression Quantification",
      workflow.type =
        "STAR - Counts"
    ),

    error = function(e) {

      cat(
        "Query failed:",
        conditionMessage(e),
        "\n"
      )

      NULL
    }
  )

  if (is.null(query)) {

    summary_list[[project]] <- data.frame(
      Project = project,
      Total_STAR_files = 0,
      Total_unique_cases = 0,
      Primary_Tumor_cases = 0,
      Solid_Tissue_Normal_cases = 0,
      Primary_Blood_Cancer_cases = 0,
      Recurrent_Tumor_cases = 0,
      Metastatic_cases = 0,
      Has_Primary_Tumor = FALSE,
      Has_Solid_Normal = FALSE,
      Tumor_Normal_available = FALSE,
      stringsAsFactors = FALSE
    )

    next
  }

  res <- TCGAbiolinks::getResults(
    query
  )

  if (
    nrow(res) == 0 ||
    !("sample_type" %in% colnames(res))
  ) {

    cat("No usable STAR-count records.\n")

    next
  }


  # ----------------------------------------------------------
  # Keep only metadata needed for the audit
  # ----------------------------------------------------------

  audit_cols <- intersect(
    c(
      "project",
      "cases",
      "sample_type",
      "file_id",
      "file_name",
      "analysis_workflow_type"
    ),
    colnames(res)
  )

  detail <- res[
    ,
    audit_cols,
    drop = FALSE
  ]

  detail$TCGA_project <- project

  detail_list[[project]] <- detail


  # ----------------------------------------------------------
  # Use unique CASE + SAMPLE TYPE combinations.
  #
  # This avoids counting multiple GDC files for the same
  # patient/sample type as separate biological cases.
  # ----------------------------------------------------------

  case_type <- unique(
    data.frame(
      cases =
        as.character(res$cases),
      sample_type =
        as.character(res$sample_type),
      stringsAsFactors = FALSE
    )
  )

  case_type <- case_type[
    !is.na(case_type$cases) &
    !is.na(case_type$sample_type),
  ]


  cat("\nSample types:\n")

  sample_table <- sort(
    table(case_type$sample_type),
    decreasing = TRUE
  )

  print(sample_table)


  # ----------------------------------------------------------
  # Counts
  # ----------------------------------------------------------

  n_primary <- sum(
    case_type$sample_type ==
      "Primary Tumor"
  )

  n_normal <- sum(
    case_type$sample_type ==
      "Solid Tissue Normal"
  )

  n_blood <- sum(
    grepl(
      "^Primary Blood Derived Cancer",
      case_type$sample_type
    )
  )

  n_recurrent <- sum(
    grepl(
      "Recurrent",
      case_type$sample_type,
      ignore.case = TRUE
    )
  )

  n_metastatic <- sum(
    grepl(
      "Metastatic",
      case_type$sample_type,
      ignore.case = TRUE
    )
  )

  summary_list[[project]] <- data.frame(

    Project = project,

    Total_STAR_files =
      nrow(res),

    Total_unique_cases =
      length(unique(case_type$cases)),

    Primary_Tumor_cases =
      n_primary,

    Solid_Tissue_Normal_cases =
      n_normal,

    Primary_Blood_Cancer_cases =
      n_blood,

    Recurrent_Tumor_cases =
      n_recurrent,

    Metastatic_cases =
      n_metastatic,

    Has_Primary_Tumor =
      n_primary > 0,

    Has_Solid_Normal =
      n_normal > 0,

    Tumor_Normal_available =
      n_primary > 0 &
      n_normal > 0,

    stringsAsFactors = FALSE
  )
}


# ============================================================
# 4. SAVE PAN-CANCER SUMMARY
# ============================================================

summary_table <- do.call(
  rbind,
  summary_list
)

rownames(summary_table) <- NULL

summary_table <- summary_table[
  order(summary_table$Project),
]

write.csv(
  summary_table,
  file.path(
    outdir,
    "13a_TCGA_STARCounts_sample_audit.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 5. SAVE FULL SAMPLE-TYPE DETAILS
# ============================================================

if (length(detail_list) > 0) {

  detail_table <- do.call(
    rbind,
    detail_list
  )

  rownames(detail_table) <- NULL

  write.csv(
    detail_table,
    file.path(
      outdir,
      "13a_TCGA_STARCounts_metadata.csv"
    ),
    row.names = FALSE
  )
}


# ============================================================
# 6. CONSOLE SUMMARY
# ============================================================

cat("\n\n========================================\n")
cat("PAN-CANCER AUDIT SUMMARY\n")
cat("========================================\n\n")

print(
  summary_table,
  row.names = FALSE
)

cat("\nTCGA projects:",
    nrow(summary_table),
    "\n")

cat(
  "Projects with Primary Tumor:",
  sum(summary_table$Has_Primary_Tumor),
  "\n"
)

cat(
  "Projects with Solid Tissue Normal:",
  sum(summary_table$Has_Solid_Normal),
  "\n"
)

cat(
  "Projects with both Primary Tumor + Normal:",
  sum(summary_table$Tumor_Normal_available),
  "\n"
)

cat(
  "Projects containing primary blood-derived cancer:",
  sum(summary_table$Primary_Blood_Cancer_cases > 0),
  "\n"
)

cat("\nSaved:\n")
cat(
  " ",
  file.path(
    outdir,
    "13a_TCGA_STARCounts_sample_audit.csv"
  ),
  "\n"
)

cat(
  " ",
  file.path(
    outdir,
    "13a_TCGA_STARCounts_metadata.csv"
  ),
  "\n"
)

cat("\n========================================\n")
cat("13a COMPLETE\n")
cat("========================================\n")
