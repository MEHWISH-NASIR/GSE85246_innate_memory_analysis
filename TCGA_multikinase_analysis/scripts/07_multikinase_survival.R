
library(data.table)
library(survival)

setDTthreads(1)

root <- "TCGA_multikinase_analysis"
outdir <- file.path(root, "results/survival")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("1. Loading multikinase patient-level expression...\n")

# Patient-level primary malignant data generated in Step 06
expr <- fread(
  file.path(
    root,
    "results/stage_grade",
    "06_multikinase_stage_grade_patient_level.csv"
  )
)

cat("Expression rows:", nrow(expr), "\n")
cat("Genes:", uniqueN(expr$gene), "\n")

# ============================================================
# TCGA-CDR reference metadata
# ============================================================

cat("2. Loading TCGA-CDR reference metadata...\n")

cdr <- fread(
  "TCGA_JAK3_pan_cancer/results/JAK3_TCGA_CDR_patient_data.csv"
)

# Keep survival/clinical fields only.
# Remove JAK3-specific expression fields.
cdr_meta <- unique(
  cdr[, .(
    project,
    patient_id = cases.submitter_id,

    type,

    OS,
    OS.time,

    DSS,
    DSS.time,

    DFI,
    DFI.time,

    PFI,
    PFI.time,

    Redacted,
    OS_ready,
    DSS_ready,
    DFI_ready,
    PFI_ready,
    survival_eligible
  )]
)

stopifnot(
  !anyDuplicated(
    cdr_meta[, .(project, patient_id)]
  )
)

cat("TCGA-CDR patients:", nrow(cdr_meta), "\n")

# ============================================================
# Merge survival metadata with all genes
# ============================================================

m <- merge(
  expr,
  cdr_meta,
  by = c("project", "patient_id"),
  all.x = TRUE,
  sort = FALSE
)

cat(
  "Matched expression rows:",
  sum(!is.na(m$type)),
  "/",
  nrow(m),
  "\n"
)

cat(
  "Unmatched rows:",
  sum(is.na(m$type)),
  "\n"
)

fwrite(
  m,
  file.path(
    outdir,
    "07_multikinase_TCGA_CDR_patient_data.csv"
  )
)

# ============================================================
# Endpoint assignment
# Same rules as original JAK3 analysis
# ============================================================

m[, Survival_endpoint :=
  fifelse(
    project == "TCGA-LAML",
    "OS",

    fifelse(
      project == "TCGA-PCPG",
      NA_character_,
      "PFI"
    )
  )
]

m[, Endpoint_recommendation :=
  fifelse(
    project == "TCGA-LAML",
    "OS used because PFI unavailable",

    fifelse(
      project == "TCGA-PCPG",
      "No TCGA-CDR survival endpoint recommended",

      fifelse(
        project %in% c(
          "TCGA-DLBC",
          "TCGA-KICH"
        ),
        "PFI use with caution",
        "PFI"
      )
    )
  )
]

# ============================================================
# Cox fitting function
# ============================================================

fit_one <- function(dat, gene_id, project_id) {

  endpoint <- unique(dat$Survival_endpoint)

  endpoint <- endpoint[!is.na(endpoint)]

  if (length(endpoint) == 0) {

    return(
      data.table(
        gene = gene_id,
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
        Analysis_status = "Not analysed"
      )
    )
  }

  endpoint <- endpoint[1]

  if (endpoint == "PFI") {

    x <- dat[
      !(Redacted %in% TRUE) &
      PFI_ready %in% TRUE &
      !is.na(log2TPM1),
      .(
        time = PFI.time,
        event = PFI,
        log2TPM1
      )
    ]

  } else {

    x <- dat[
      !(Redacted %in% TRUE) &
      OS_ready %in% TRUE &
      !is.na(log2TPM1),
      .(
        time = OS.time,
        event = OS,
        log2TPM1
      )
    ]
  }

  x <- x[
    !is.na(time) &
    time >= 0 &
    event %in% c(0, 1)
  ]

  n_patients <- nrow(x)
  n_events <- sum(x$event == 1)

  # Same minimum requirement as JAK3
  if (
    n_patients < 10 ||
    n_events < 10
  ) {

    return(
      data.table(
        gene = gene_id,
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
        Analysis_status = "Insufficient events"
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
      data.table(
        gene = gene_id,
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
        Analysis_status = "Cox failed"
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
    ph$table[
      "log2TPM1",
      "p"
    ]
  }

  data.table(
    gene = gene_id,
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
      !is.na(PH_P) &&
      PH_P < 0.05,
    Low_event_n =
      n_events < 20,
    Analysis_status = "Analysed"
  )
}

# ============================================================
# Run all gene × cancer models
# ============================================================

groups <- unique(
  m[, .(gene, project)]
)

results_list <- vector(
  "list",
  nrow(groups)
)

cat(
  "3. Running",
  nrow(groups),
  "gene-cancer survival models...\n"
)

for (i in seq_len(nrow(groups))) {

  g <- groups$gene[i]
  cancer <- groups$project[i]

  results_list[[i]] <- fit_one(
    m[
      gene == g &
      project == cancer
    ],
    g,
    cancer
  )

  if (i %% 50 == 0) {
    cat(
      "Processed:",
      i,
      "/",
      nrow(groups),
      "\n"
    )
  }
}

results <- rbindlist(
  results_list,
  fill = TRUE
)

# ============================================================
# Gene-specific BH FDR
# Same logic as original JAK3
# ============================================================

results[, Cox_FDR := NA_real_]

results[
  Analysis_status == "Analysed" &
  !is.na(Cox_P),
  Cox_FDR :=
    p.adjust(
      Cox_P,
      method = "BH"
    ),
  by = gene
]

# ============================================================
# Additional multikinase-wide FDR
# ============================================================

results[, Multikinase_Global_FDR := NA_real_]

idx <- which(
  results$Analysis_status == "Analysed" &
  !is.na(results$Cox_P)
)

results$Multikinase_Global_FDR[idx] <-
  p.adjust(
    results$Cox_P[idx],
    method = "BH"
  )

# ============================================================
# Direction
# ============================================================

results[, Direction :=
  fifelse(
    is.na(HR),
    NA_character_,

    fifelse(
      HR > 1,
      "Higher expression = higher hazard",

      fifelse(
        HR < 1,
        "Higher expression = lower hazard",
        "Neutral"
      )
    )
  )
]

setorder(
  results,
  gene,
  Cox_FDR
)

fwrite(
  results,
  file.path(
    outdir,
    "07_multikinase_TCGA_CDR_survival.csv"
  )
)

# ============================================================
# Gene-level summary
# ============================================================

summary_gene <- results[
  ,
  .(
    Projects_total = .N,

    Projects_analysed =
      sum(
        Analysis_status == "Analysed"
      ),

    Cox_FDR_significant =
      sum(
        Cox_FDR < 0.05,
        na.rm = TRUE
      ),

    Global_FDR_significant =
      sum(
        Multikinase_Global_FDR < 0.05,
        na.rm = TRUE
      ),

    Higher_hazard_significant =
      sum(
        Cox_FDR < 0.05 &
        HR > 1,
        na.rm = TRUE
      ),

    Lower_hazard_significant =
      sum(
        Cox_FDR < 0.05 &
        HR < 1,
        na.rm = TRUE
      ),

    PH_violations =
      sum(
        PH_violation %in% TRUE,
        na.rm = TRUE
      ),

    Low_event_analyses =
      sum(
        Low_event_n %in% TRUE,
        na.rm = TRUE
      )
  ),
  by = gene
]

summary_gene[, Reference :=
  gene == "JAK3"
]

setorder(
  summary_gene,
  -Cox_FDR_significant,
  -Global_FDR_significant
)

fwrite(
  summary_gene,
  file.path(
    outdir,
    "07_survival_gene_comparison_summary.csv"
  )
)

# ============================================================
# Console summary
# ============================================================

cat("\n========================================\n")
cat("MULTIKINASE TCGA-CDR SURVIVAL ANALYSIS\n")
cat("========================================\n")

cat("\n===== GENE-LEVEL SUMMARY =====\n")
print(summary_gene)

cat("\n===== JAK3 SURVIVAL RESULTS =====\n")

print(
  results[
    gene == "JAK3"
  ]
)

cat("\n07 COMPLETE\n")

