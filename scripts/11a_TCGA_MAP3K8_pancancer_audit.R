# ============================================================
# 11a_TCGA_MAP3K8_pancancer_audit.R
#
# MAP3K8 PAN-CANCER EXPLORATION
#
# Source:
# UCSC Xena TCGA-TARGET-GTEx Toil RNA-seq
#
# Purpose of this first step:
#   - query MAP3K8 only
#   - inspect cancer / tissue annotations
#   - inspect tumor and normal sample structure
#
# NO overexpression threshold is applied yet.
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Package
# ------------------------------------------------------------

if (!requireNamespace(
    "UCSCXenaShiny",
    quietly = TRUE
)) {

    install.packages(
        "UCSCXenaShiny",
        repos = "https://cloud.r-project.org"
    )
}

# ------------------------------------------------------------
# 2. Output
# ------------------------------------------------------------

outdir <- "results/11_TCGA_MAP3K8"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

cat("\n========================================\n")
cat("11a — MAP3K8 PAN-CANCER AUDIT\n")
cat("========================================\n\n")

# ------------------------------------------------------------
# 3. Query MAP3K8 from uniformly processed
#    TCGA-TARGET-GTEx Toil dataset
# ------------------------------------------------------------

map3k8 <- UCSCXenaShiny::query_toil_value_df(
    identifier = "MAP3K8"
)

map3k8 <- as.data.frame(
    map3k8,
    stringsAsFactors = FALSE
)

cat(
    "Samples returned:",
    nrow(map3k8),
    "\n"
)

cat(
    "Columns returned:",
    ncol(map3k8),
    "\n\n"
)

cat("Column names:\n\n")
print(names(map3k8))

cat("\nFirst rows:\n\n")

print(
    head(
        map3k8,
        10
    )
)

# ------------------------------------------------------------
# 4. Save untouched query result
# ------------------------------------------------------------

write.csv(
    map3k8,
    file.path(
        outdir,
        "11a_MAP3K8_TCGA_TARGET_GTEx_Toil.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 5. Study composition
# ------------------------------------------------------------

if ("study" %in% names(map3k8)) {

    cat("\n========================================\n")
    cat("STUDIES\n")
    cat("========================================\n\n")

    print(
        sort(
            table(map3k8$study),
            decreasing = TRUE
        )
    )
}

# ------------------------------------------------------------
# 6. Sample types
# ------------------------------------------------------------

if ("sample_type" %in% names(map3k8)) {

    cat("\n========================================\n")
    cat("SAMPLE TYPES\n")
    cat("========================================\n\n")

    print(
        sort(
            table(map3k8$sample_type),
            decreasing = TRUE
        )
    )
}

# ------------------------------------------------------------
# 7. Cancer / tissue categories
# ------------------------------------------------------------

if ("detailed_category" %in% names(map3k8)) {

    cat("\n========================================\n")
    cat("DETAILED CATEGORIES\n")
    cat("========================================\n\n")

    print(
        sort(
            table(map3k8$detailed_category),
            decreasing = TRUE
        )
    )
}

# ------------------------------------------------------------
# 8. Primary sites
# ------------------------------------------------------------

if ("primary_site" %in% names(map3k8)) {

    cat("\n========================================\n")
    cat("PRIMARY SITES\n")
    cat("========================================\n\n")

    print(
        sort(
            table(map3k8$primary_site),
            decreasing = TRUE
        )
    )
}

# ------------------------------------------------------------
# 9. Expression audit
# ------------------------------------------------------------

if ("expression" %in% names(map3k8)) {

    cat("\n========================================\n")
    cat("MAP3K8 EXPRESSION\n")
    cat("========================================\n\n")

    cat(
        "Expression range:",
        paste(
            round(
                range(
                    map3k8$expression,
                    na.rm = TRUE
                ),
                3
            ),
            collapse = " to "
        ),
        "\n"
    )

    cat(
        "Median expression:",
        round(
            median(
                map3k8$expression,
                na.rm = TRUE
            ),
            3
        ),
        "\n"
    )
}

cat("\n========================================\n")
cat("11a COMPLETE\n")
cat("========================================\n\n")

cat(
    "Next step:\n",
    "separate TCGA tumors and corresponding normal tissues,\n",
    "then calculate tumor-vs-normal expression and\n",
    "MAP3K8-high frequency by cancer type.\n"
)
