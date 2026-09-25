library(dplyr)
library(survival)

d <- read.csv(
  "TCGA_JAK3_pan_cancer/results/JAK3_TCGA_CDR_patient_data.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("JAK3 TCGA-CDR SURVIVAL ANALYSIS\n")
cat("========================================\n")

# ------------------------------------------------------------
# Endpoint selection according to TCGA-CDR
# ------------------------------------------------------------

d <- d %>%
  mutate(
    Survival_endpoint = case_when(
      project == "TCGA-LAML" ~ "OS",
      project == "TCGA-PCPG" ~ NA_character_,
      TRUE ~ "PFI"
    ),

    Endpoint_recommendation = case_when(
      project == "TCGA-LAML" ~
        "OS used because PFI unavailable",

      project == "TCGA-PCPG" ~
        "No TCGA-CDR survival endpoint recommended",

      project %in% c(
        "TCGA-DLBC",
        "TCGA-KICH"
      ) ~
        "PFI use with caution",

      TRUE ~
        "PFI"
    )
  )

# ------------------------------------------------------------
# Fit one continuous-expression Cox model per cancer
# ------------------------------------------------------------

fit_one <- function(dat, project_id) {

  endpoint <- unique(dat$Survival_endpoint)

  if (length(endpoint) != 1 || is.na(endpoint)) {

    return(
      data.frame(
        project = project_id,
        Endpoint = NA_character_,
        Recommendation =
          unique(dat$Endpoint_recommendation)[1],
        N = NA_integer_,
        Events = NA_integer_,
        HR = NA_real_,
        CI95_low = NA_real_,
        CI95_high = NA_real_,
        Cox_P = NA_real_,
        PH_P = NA_real_,
        PH_violation = NA,
        Low_event_n = NA,
        Analysis_status = "Not analysed",
        stringsAsFactors = FALSE
      )
    )
  }

  if (endpoint == "PFI") {

    x <- dat %>%
      filter(
        !Redacted %in% TRUE,
        PFI_ready,
        !is.na(log2TPM1)
      ) %>%
      transmute(
        time = PFI.time,
        event = PFI,
        log2TPM1 = log2TPM1
      )

  } else {

    x <- dat %>%
      filter(
        !Redacted %in% TRUE,
        OS_ready,
        !is.na(log2TPM1)
      ) %>%
      transmute(
        time = OS.time,
        event = OS,
        log2TPM1 = log2TPM1
      )
  }

  x <- x %>%
    filter(
      !is.na(time),
      time >= 0,
      event %in% c(0, 1)
    )

  n_patients <- nrow(x)
  n_events <- sum(x$event == 1)

  # Minimum event requirement
  if (n_patients < 10 || n_events < 10) {

    return(
      data.frame(
        project = project_id,
        Endpoint = endpoint,
        Recommendation =
          unique(dat$Endpoint_recommendation)[1],
        N = n_patients,
        Events = n_events,
        HR = NA_real_,
        CI95_low = NA_real_,
        CI95_high = NA_real_,
        Cox_P = NA_real_,
        PH_P = NA_real_,
        PH_violation = NA,
        Low_event_n = n_events < 20,
        Analysis_status =
          "Insufficient events",
        stringsAsFactors = FALSE
      )
    )
  }

  fit <- tryCatch(
    coxph(
      Surv(time, event) ~ log2TPM1,
      data = x,
      ties = "efron"
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {

    return(
      data.frame(
        project = project_id,
        Endpoint = endpoint,
        Recommendation =
          unique(dat$Endpoint_recommendation)[1],
        N = n_patients,
        Events = n_events,
        HR = NA_real_,
        CI95_low = NA_real_,
        CI95_high = NA_real_,
        Cox_P = NA_real_,
        PH_P = NA_real_,
        PH_violation = NA,
        Low_event_n = n_events < 20,
        Analysis_status = "Cox failed",
        stringsAsFactors = FALSE
      )
    )
  }

  s <- summary(fit)

  HR <- s$coefficients[
    "log2TPM1",
    "exp(coef)"
  ]

  Cox_P <- s$coefficients[
    "log2TPM1",
    "Pr(>|z|)"
  ]

  CI_low <- s$conf.int[
    "log2TPM1",
    "lower .95"
  ]

  CI_high <- s$conf.int[
    "log2TPM1",
    "upper .95"
  ]

  ph <- tryCatch(
    cox.zph(fit),
    error = function(e) NULL
  )

  PH_P <- if (is.null(ph)) {
    NA_real_
  } else {
    ph$table["log2TPM1", "p"]
  }

  data.frame(
    project = project_id,
    Endpoint = endpoint,
    Recommendation =
      unique(dat$Endpoint_recommendation)[1],
    N = n_patients,
    Events = n_events,
    HR = HR,
    CI95_low = CI_low,
    CI95_high = CI_high,
    Cox_P = Cox_P,
    PH_P = PH_P,
    PH_violation =
      !is.na(PH_P) && PH_P < 0.05,
    Low_event_n =
      n_events < 20,
    Analysis_status = "Analysed",
    stringsAsFactors = FALSE
  )
}

projects <- sort(unique(d$project))

results <- lapply(
  projects,
  function(p) {
    fit_one(
      d %>% filter(project == p),
      p
    )
  }
) %>%
  bind_rows()

# ------------------------------------------------------------
# FDR only among fitted models
# ------------------------------------------------------------

results$Cox_FDR <- NA_real_

idx <- which(
  results$Analysis_status == "Analysed" &
  !is.na(results$Cox_P)
)

results$Cox_FDR[idx] <- p.adjust(
  results$Cox_P[idx],
  method = "BH"
)

results <- results %>%
  mutate(
    Direction = case_when(
      is.na(HR) ~ NA_character_,
      HR > 1 ~ "Higher JAK3 = higher hazard",
      HR < 1 ~ "Higher JAK3 = lower hazard",
      TRUE ~ "Neutral"
    )
  ) %>%
  arrange(
    is.na(Cox_FDR),
    Cox_FDR
  )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  results,
  "TCGA_JAK3_pan_cancer/results/JAK3_TCGA_CDR_survival_complete.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

cat("\n===== SURVIVAL SUMMARY =====\n")

cat(
  "Projects analysed:",
  sum(results$Analysis_status == "Analysed"),
  "\n"
)

cat(
  "Not analysed:",
  sum(results$Analysis_status != "Analysed"),
  "\n"
)

cat(
  "Cox FDR < 0.05:",
  sum(
    results$Cox_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "PH violations:",
  sum(
    results$PH_violation %in% TRUE,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Low-event analyses (<20 events):",
  sum(
    results$Low_event_n %in% TRUE,
    na.rm = TRUE
  ),
  "\n"
)

cat("\n===== SURVIVAL RESULTS =====\n")

print(
  as.data.frame(results),
  row.names = FALSE
)

cat("\n51c SURVIVAL ANALYSIS COMPLETE\n")
