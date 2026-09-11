# ============================================================
# 13h — TCGA-KIRC clinical metadata
# Stage / grade / TNM / basic clinical variables
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required.")
}

outdir <- "results/13_TCGA_priority_genes"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y = NA) {
  if (is.null(x) || length(x) == 0) y else x
}

# ------------------------------------------------------------
# Query TCGA-KIRC cases from GDC
# ------------------------------------------------------------

filter_obj <- list(
  op = "in",
  content = list(
    field = "project.project_id",
    value = list("TCGA-KIRC")
  )
)

filters <- jsonlite::toJSON(
  filter_obj,
  auto_unbox = TRUE
)

url <- paste0(
  "https://api.gdc.cancer.gov/cases?",
  "filters=",
  URLencode(filters, reserved = TRUE),
  "&expand=diagnoses,demographic",
  "&size=1000"
)

cat("Downloading TCGA-KIRC clinical metadata...\n")

gdc <- jsonlite::fromJSON(
  url,
  simplifyVector = FALSE
)

hits <- gdc$data$hits

cat("Cases returned:", length(hits), "\n")

# ------------------------------------------------------------
# Extract one primary clinical record per case
# ------------------------------------------------------------

rows <- lapply(hits, function(case) {

  case_id <-
    case$submitter_id %||% NA_character_

  diagnoses <- case$diagnoses

  if (
    is.null(diagnoses) ||
    length(diagnoses) == 0
  ) {

    d <- list()

  } else {

    d <- diagnoses[[1]]
  }

  demo <- case$demographic

  if (is.null(demo)) demo <- list()

  data.frame(
    patient_id =
      case_id,

    primary_diagnosis =
      d$primary_diagnosis %||% NA_character_,

    pathologic_stage =
      d$ajcc_pathologic_stage %||% NA_character_,

    tumor_grade =
      d$tumor_grade %||% NA_character_,

    pathologic_T =
      d$ajcc_pathologic_t %||% NA_character_,

    pathologic_N =
      d$ajcc_pathologic_n %||% NA_character_,

    pathologic_M =
      d$ajcc_pathologic_m %||% NA_character_,

    age_at_diagnosis_days =
      d$age_at_diagnosis %||% NA_real_,

    vital_status =
      demo$vital_status %||% NA_character_,

    days_to_death =
      demo$days_to_death %||% NA_real_,

    gender =
      demo$sex_at_birth %||% NA_character_,

    stringsAsFactors = FALSE
  )
})

clin <- do.call(rbind, rows)

clin$age_at_diagnosis_years <-
  suppressWarnings(
    as.numeric(clin$age_at_diagnosis_days) / 365.25
  )

# ------------------------------------------------------------
# Save raw clinical table
# ------------------------------------------------------------

outfile <- file.path(
  outdir,
  "13h_KIRC_clinical_metadata.csv"
)

write.csv(
  clin,
  outfile,
  row.names = FALSE
)

# ------------------------------------------------------------
# Audit stage / grade availability
# ------------------------------------------------------------

cat("\n========================================\n")
cat("KIRC CLINICAL AUDIT\n")
cat("========================================\n\n")

cat("Patients:", nrow(clin), "\n")
cat(
  "Unique patients:",
  length(unique(clin$patient_id)),
  "\n"
)

cat("\nPATHOLOGIC STAGE\n")
print(
  sort(
    table(
      clin$pathologic_stage,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)

cat("\nTUMOR GRADE\n")
print(
  sort(
    table(
      clin$tumor_grade,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)

cat("\nPATHOLOGIC T\n")
print(
  sort(
    table(
      clin$pathologic_T,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)

cat("\nPATHOLOGIC N\n")
print(
  sort(
    table(
      clin$pathologic_N,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)

cat("\nPATHOLOGIC M\n")
print(
  sort(
    table(
      clin$pathologic_M,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)

cat("\nVITAL STATUS\n")
print(
  table(
    clin$vital_status,
    useNA = "ifany"
  )
)

cat("\nSaved:\n ", outfile, "\n")

cat("\n========================================\n")
cat("13h COMPLETE\n")
cat("========================================\n")
