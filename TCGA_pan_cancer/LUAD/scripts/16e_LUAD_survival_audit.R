# ============================================================
# 15e — TCGA-LUAD survival metadata audit
# Corrected:
# follow_ups are CASE-level in current GDC API
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required.")
}

outdir <- "results/16_TCGA_LUAD"

`%||%` <- function(x, y = NULL) {
  if (is.null(x) || length(x) == 0) y else x
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

valid_time <- function(x) {
  x <- safe_num(x)
  x[is.finite(x) & x >= 0]
}

filter_obj <- list(
  op = "in",
  content = list(
    field = "project.project_id",
    value = list("TCGA-LUAD")
  )
)

filters <- jsonlite::toJSON(
  filter_obj,
  auto_unbox = TRUE
)

# IMPORTANT:
# follow_ups is expanded at CASE level, not diagnoses.follow_ups
url <- paste0(
  "https://api.gdc.cancer.gov/cases?",
  "filters=", URLencode(filters, reserved = TRUE),
  "&expand=diagnoses,demographic,follow_ups",
  "&size=2000"
)

cat("Downloading corrected LUAD survival metadata...\n")

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
  if (is.null(diagnoses)) diagnoses <- list()

  follow_ups <- case$follow_ups
  if (is.null(follow_ups)) follow_ups <- list()

  # ----------------------------------------------------------
  # Diagnosis-level information
  # ----------------------------------------------------------

  age_vals <- numeric(0)
  stage_vals <- character(0)

  diagnosis_follow_vals <- numeric(0)
  diagnosis_death_vals <- numeric(0)

  for (d in diagnoses) {

    age_vals <- c(
      age_vals,
      valid_time(
        d$age_at_diagnosis %||% NA
      )
    )

    st <- d$ajcc_pathologic_stage %||% NA_character_

    if (
      length(st) > 0 &&
      !is.na(st) &&
      st != ""
    ) {
      stage_vals <- c(stage_vals, st)
    }

    diagnosis_follow_vals <- c(
      diagnosis_follow_vals,
      valid_time(
        d$days_to_last_follow_up %||% NA
      )
    )

    diagnosis_death_vals <- c(
      diagnosis_death_vals,
      valid_time(
        d$days_to_death %||% NA
      )
    )
  }

  # ----------------------------------------------------------
  # CASE-level follow-up records
  # ----------------------------------------------------------

  case_follow_vals <- numeric(0)
  followup_death_vals <- numeric(0)

  if (length(follow_ups) > 0) {

    for (f in follow_ups) {

      case_follow_vals <- c(
        case_follow_vals,
        valid_time(
          f$days_to_follow_up %||% NA
        )
      )

      # Retain if present in a given follow-up record
      followup_death_vals <- c(
        followup_death_vals,
        valid_time(
          f$days_to_death %||% NA
        )
      )
    }
  }

  # ----------------------------------------------------------
  # Death time
  # ----------------------------------------------------------

  demographic_death <- valid_time(
    demo$days_to_death %||% NA
  )

  death_candidates <- c(
    demographic_death,
    diagnosis_death_vals,
    followup_death_vals
  )

  death_days <- if (
    length(death_candidates) > 0
  ) {
    max(death_candidates)
  } else {
    NA_real_
  }

  # ----------------------------------------------------------
  # Last known follow-up time
  #
  # Use greatest available value from:
  # 1. diagnoses.days_to_last_follow_up
  # 2. follow_ups.days_to_follow_up
  # ----------------------------------------------------------

  follow_candidates <- c(
    diagnosis_follow_vals,
    case_follow_vals
  )

  follow_days <- if (
    length(follow_candidates) > 0
  ) {
    max(follow_candidates)
  } else {
    NA_real_
  }

  # ----------------------------------------------------------
  # Survival event
  # ----------------------------------------------------------

  vital <- demo$vital_status %||% NA_character_

  event <- ifelse(
    toupper(vital) == "DEAD",
    1L,
    ifelse(
      toupper(vital) == "ALIVE",
      0L,
      NA_integer_
    )
  )

  # Deaths use days_to_death.
  # Living patients are censored at latest available follow-up.
  OS_days <- if (
    !is.na(event) &&
    event == 1
  ) {

    death_days

  } else if (
    !is.na(event) &&
    event == 0
  ) {

    follow_days

  } else {

    NA_real_
  }

  data.frame(
    patient_id =
      case$submitter_id %||% NA_character_,

    vital_status =
      vital,

    survival_event =
      event,

    days_to_death =
      death_days,

    diagnosis_last_followup_days =
      if (length(diagnosis_follow_vals) > 0)
        max(diagnosis_follow_vals)
      else
        NA_real_,

    case_followup_days =
      if (length(case_follow_vals) > 0)
        max(case_follow_vals)
      else
        NA_real_,

    max_followup_days =
      follow_days,

    OS_days =
      OS_days,

    OS_years =
      OS_days / 365.25,

    age_at_diagnosis_years =
      if (length(age_vals) > 0)
        age_vals[1] / 365.25
      else
        NA_real_,

    gender =
      demo$sex_at_birth %||% NA_character_,

    pathologic_stage =
      if (length(stage_vals) > 0)
        stage_vals[1]
      else
        NA_character_,

    stringsAsFactors = FALSE
  )
})

surv <- do.call(
  rbind,
  rows
)

outfile <- file.path(
  outdir,
  "16e_LUAD_survival_metadata.csv"
)

write.csv(
  surv,
  outfile,
  row.names = FALSE
)

# ============================================================
# AUDIT
# ============================================================

cat("\n========================================\n")
cat("CORRECTED LUAD SURVIVAL AUDIT\n")
cat("========================================\n\n")

cat("Total cases:", nrow(surv), "\n")
cat(
  "Unique patients:",
  length(unique(surv$patient_id)),
  "\n"
)

cat("\nVital status:\n")
print(
  table(
    surv$vital_status,
    useNA = "ifany"
  )
)

cat("\nAvailability of follow-up source:\n")

cat(
  "Diagnosis-level follow-up available:",
  sum(!is.na(surv$diagnosis_last_followup_days)),
  "\n"
)

cat(
  "Case-level follow-up available:",
  sum(!is.na(surv$case_followup_days)),
  "\n"
)

cat(
  "Any follow-up available:",
  sum(!is.na(surv$max_followup_days)),
  "\n"
)

cat("\nOS availability by vital status:\n")

print(
  table(
    surv$vital_status,
    !is.na(surv$OS_days),
    useNA = "ifany"
  )
)

usable <- surv[
  !is.na(surv$OS_days) &
  !is.na(surv$survival_event) &
  surv$OS_days >= 0,
]

cat(
  "\nUsable OS patients:",
  nrow(usable),
  "\n"
)

cat(
  "Deaths:",
  sum(usable$survival_event == 1),
  "\n"
)

cat(
  "Censored:",
  sum(usable$survival_event == 0),
  "\n"
)

cat("\nOS years summary:\n")
print(
  summary(
    usable$OS_years
  )
)

cat("\nMissingness among usable OS patients:\n")

cat(
  "Age missing:",
  sum(is.na(usable$age_at_diagnosis_years)),
  "\n"
)

cat(
  "Sex missing:",
  sum(
    is.na(usable$gender) |
    usable$gender == ""
  ),
  "\n"
)

cat(
  "Stage missing:",
  sum(
    is.na(usable$pathologic_stage) |
    usable$pathologic_stage == ""
  ),
  "\n"
)

cat("\nSaved:\n ", outfile, "\n")

cat("\n========================================\n")
cat("16e CORRECTED COMPLETE\n")
cat("========================================\n")
