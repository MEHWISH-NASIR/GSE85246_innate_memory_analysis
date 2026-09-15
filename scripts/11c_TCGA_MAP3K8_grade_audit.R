# ============================================================
# 11c_TCGA_MAP3K8_grade_audit.R
#
# Audit TCGA clinical grade data before MAP3K8 grade analysis.
# No statistical grade testing yet.
# ============================================================

rm(list = ls())

if (!requireNamespace("UCSCXenaShiny", quietly = TRUE)) {
    stop("UCSCXenaShiny is required.")
}

outdir <- "results/11_TCGA_MAP3K8"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

cat("\n========================================\n")
cat("11c — TCGA GRADE DATA AUDIT\n")
cat("========================================\n\n")

# ------------------------------------------------------------
# 1. Load cleaned TCGA clinical data
# ------------------------------------------------------------

data(
    "tcga_clinical_fine",
    package = "UCSCXenaShiny"
)

clin <- tcga_clinical_fine

cat(
    "Clinical rows:",
    nrow(clin),
    "\n"
)

cat(
    "Clinical columns:",
    ncol(clin),
    "\n\n"
)

cat("Column names:\n\n")

print(
    names(clin)
)

# ------------------------------------------------------------
# 2. Identify grade / stage / sample / cancer columns
# ------------------------------------------------------------

interesting <- grep(
    "sample|cancer|tissue|grade|stage|code",
    names(clin),
    ignore.case = TRUE,
    value = TRUE
)

cat("\n========================================\n")
cat("POTENTIALLY RELEVANT COLUMNS\n")
cat("========================================\n\n")

print(interesting)

if (length(interesting) > 0) {

    cat("\nExample clinical rows:\n\n")

    print(
        head(
            clin[
                ,
                interesting,
                drop = FALSE
            ],
            10
        )
    )
}

# ------------------------------------------------------------
# 3. Grade column
# ------------------------------------------------------------

grade_col <- grep(
    "^Grade$|grade",
    names(clin),
    ignore.case = TRUE,
    value = TRUE
)

cat("\n========================================\n")
cat("GRADE COLUMN CANDIDATES\n")
cat("========================================\n\n")

print(grade_col)

if (length(grade_col) == 0) {
    stop("No grade-like column found.")
}

# Use exact Grade if available
if ("Grade" %in% names(clin)) {

    grade_name <- "Grade"

} else {

    grade_name <- grade_col[1]
}

cat(
    "\nUsing provisional grade column:",
    grade_name,
    "\n"
)

cat("\nOverall grade values:\n\n")

print(
    sort(
        table(
            clin[[grade_name]],
            useNA = "ifany"
        ),
        decreasing = TRUE
    )
)

# ------------------------------------------------------------
# 4. Determine cancer-code column
# ------------------------------------------------------------

candidate_cancer_cols <- c(
    "Cancer",
    "tissue",
    "TCGA"
)

cancer_col <- candidate_cancer_cols[
    candidate_cancer_cols %in%
        names(clin)
]

if (length(cancer_col) == 0) {

    cat(
        "\nCould not automatically identify cancer column.\n"
    )

} else {

    cancer_col <- cancer_col[1]

    cat(
        "\nCancer column:",
        cancer_col,
        "\n"
    )

    # Keep only reported grade values
    g <- clin[
        !is.na(clin[[grade_name]]) &
        clin[[grade_name]] != "",
        ,
        drop = FALSE
    ]

    cat(
        "\nPatients/samples with reported grade:",
        nrow(g),
        "\n\n"
    )

    grade_counts <- as.data.frame(
        table(
            Cancer = g[[cancer_col]],
            Grade = g[[grade_name]]
        ),
        stringsAsFactors = FALSE
    )

    grade_counts <- grade_counts[
        grade_counts$Freq > 0,
        ,
        drop = FALSE
    ]

    grade_counts <- grade_counts[
        order(
            grade_counts$Cancer,
            grade_counts$Grade
        ),
        ,
        drop = FALSE
    ]

    write.csv(
        grade_counts,
        file.path(
            outdir,
            "11c_TCGA_grade_availability.csv"
        ),
        row.names = FALSE
    )

    cat(
        "Grade availability by cancer:\n\n"
    )

    print(
        grade_counts,
        row.names = FALSE
    )
}

cat("\n========================================\n")
cat("11c GRADE AUDIT COMPLETE\n")
cat("========================================\n")
