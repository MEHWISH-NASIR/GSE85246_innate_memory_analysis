# ============================================================
# 14e — TCGA-BLCA survival metadata audit
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required.")
}

outdir <- "results/14_TCGA_BLCA"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y = NA) {
  if (is.null(x) || length(x) == 0) y else x
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_max <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x) & !is.na(x) & x >= 0]

  if (length(x) == 0) {
    return(NA_real_)
  }

  max(x)
}

filter_obj <- list(
  op = "in",
  content = list(
    field = "project.project_id",
    value = list("TCGA-BLCA")
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
  "&expand=diagnoses,demographic,follow_ups",
  "&size=1000"
)

cat("Downloading TCGA-BLCA survival metadata...\n")

gdc <- jsonlite::fromJSON(
  url,
  simplifyVector = FALSE
)

hits <- gdc$data$hits

cat("Cases returned:", length(hits), "\n")

rows <- lapply(hits, function(case) {

  demo <- case$demographic
  if (is.null(demo)) demo <- list()

  diagnoses <- case$diagnoses
  followups <- case$follow_ups

  diag_fu <- NA_real_

  if (!is.null(diagnoses) && length(diagnoses) > 0) {

    diag_fu <- safe_max(
      sapply(
        diagnoses,
        function(d)
          d$days_to_last_follow_up %||% NA
      )
    )
  }

  case_fu <- NA_real_

  if (!is.null(followups) && length(followups) > 0) {

    case_fu <- safe_max(
      sapply(
        followups,
        function(f)
          f$days_to_follow_up %||% NA
      )
    )
  }

  last_followup <- safe_max(
    c(diag_fu, case_fu)
  )

  vital <- tolower(
    as.character(
      demo$vital_status %||% NA_character_
    )
  )

  days_death <- safe_num(
    demo$days_to_death %||% NA_real_
  )

  event <- ifelse(
    vital == "dead",
    1,
    ifelse(
      vital == "alive",
      0,
      NA
    )
  )

  os_days <- ifelse(
    event == 1,
    days_death,
    ifelse(
      event == 0,
      last_followup,
      NA_real_
    )
  )

  data.frame(
    patient_id =
      case$submitter_id %||% NA_character_,

    vital_status =
      vital,

    survival_event =
      event,

    days_to_death =
      days_death,

    diagnosis_last_followup_days =
      diag_fu,

    case_followup_days =
      case_fu,

    last_followup_days =
      last_followup,

    OS_days =
      os_days,

    stringsAsFactors = FALSE
  )
})

surv <- do.call(rbind, rows)

surv$OS_days[
  !is.na(surv$OS_days) &
  surv$OS_days <= 0
] <- NA

surv$OS_years <-
  surv$OS_days / 365.25

surv$survival_usable <-
  !is.na(surv$survival_event) &
  !is.na(surv$OS_days)

outfile <- file.path(
  outdir,
  "14e_BLCA_survival_metadata.csv"
)

write.csv(
  surv,
  outfile,
  row.names = FALSE
)

cat("\n========================================\n")
cat("BLCA SURVIVAL AUDIT\n")
cat("========================================\n\n")

cat("Total patients:", nrow(surv), "\n")

cat("\nVital status:\n")
print(
  table(
    surv$vital_status,
    useNA = "ifany"
  )
)

cat(
  "\nPatients with usable OS data:",
  sum(surv$survival_usable),
  "\n"
)

cat("\nUsable events:\n")
print(
  table(
    surv$survival_event[
      surv$survival_usable
    ]
  )
)

cat("\nOS years summary:\n")
print(
  summary(
    surv$OS_years[
      surv$survival_usable
    ]
  )
)

cat("\nSaved:\n ", outfile, "\n")

cat("\n========================================\n")
cat("14e COMPLETE\n")
cat("========================================\n")
