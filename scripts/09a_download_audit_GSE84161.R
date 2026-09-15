# ============================================================
# 09a_download_audit_GSE84161.R
#
# Independent human MAP3K8/TPL2 perturbation validation
#
# PURPOSE:
# Download and audit GSE84161 before any statistical analysis.
#
# We do NOT assume that this dataset agrees with GSE85246.
# We first verify the actual samples, conditions and expression
# matrix supplied by GEO.
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
    BiocManager::install(
        "GEOquery",
        ask = FALSE,
        update = FALSE
    )
}

library(GEOquery)
library(Biobase)

# ------------------------------------------------------------
# 2. Directories
# ------------------------------------------------------------

data_dir <- "data/external/GSE84161"
result_dir <- "results/09_MAP3K8_validation"

dir.create(
    data_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    result_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

cat("\n========================================\n")
cat("GSE84161 DOWNLOAD + AUDIT\n")
cat("========================================\n\n")

# ------------------------------------------------------------
# 3. Download GEO series
# ------------------------------------------------------------

gse_list <- GEOquery::getGEO(
    "GSE84161",
    GSEMatrix = TRUE,
    getGPL = FALSE,
    destdir = data_dir
)

cat(
    "Number of ExpressionSet objects:",
    length(gse_list),
    "\n\n"
)

# ------------------------------------------------------------
# 4. Inspect all returned platforms
# ------------------------------------------------------------

for (i in seq_along(gse_list)) {

    x <- gse_list[[i]]

    cat(
        "ExpressionSet", i, "\n",
        "Platform:", annotation(x), "\n",
        "Features:", nrow(Biobase::exprs(x)), "\n",
        "Samples:", ncol(Biobase::exprs(x)), "\n\n"
    )
}

# ------------------------------------------------------------
# 5. Use primary ExpressionSet
# ------------------------------------------------------------

if (length(gse_list) != 1) {
    cat(
        "NOTE: More than one ExpressionSet was returned.\n",
        "Do not proceed to differential analysis until checked.\n\n"
    )
}

eset <- gse_list[[1]]

expr <- Biobase::exprs(eset)
meta <- Biobase::pData(eset)
features <- Biobase::fData(eset)

# ------------------------------------------------------------
# 6. Save untouched GEO object
# ------------------------------------------------------------

saveRDS(
    eset,
    file.path(
        data_dir,
        "GSE84161_GEO_expressionSet.rds"
    )
)

# ------------------------------------------------------------
# 7. Save complete metadata
# ------------------------------------------------------------

meta$GEO_sample <- rownames(meta)

write.csv(
    meta,
    file.path(
        result_dir,
        "09a_GSE84161_complete_sample_metadata.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 8. Save feature annotation
# ------------------------------------------------------------

if (nrow(features) > 0) {

    features$PROBE_ID <- rownames(features)

    write.csv(
        features,
        file.path(
            result_dir,
            "09a_GSE84161_feature_annotation.csv"
        ),
        row.names = FALSE
    )
}

# ------------------------------------------------------------
# 9. Expression matrix audit
# ------------------------------------------------------------

cat("========================================\n")
cat("EXPRESSION MATRIX\n")
cat("========================================\n\n")

cat(
    "Rows/features:",
    nrow(expr),
    "\n"
)

cat(
    "Samples:",
    ncol(expr),
    "\n"
)

cat(
    "Expression range:",
    paste(
        round(range(expr, na.rm = TRUE), 3),
        collapse = " to "
    ),
    "\n\n"
)

# ------------------------------------------------------------
# 10. Identify informative metadata columns
# ------------------------------------------------------------

possible_cols <- grep(
    paste(
        c(
            "title",
            "geo_accession",
            "source",
            "characteristics",
            "treatment",
            "description"
        ),
        collapse = "|"
    ),
    names(meta),
    ignore.case = TRUE,
    value = TRUE
)

cat("========================================\n")
cat("RELEVANT SAMPLE METADATA\n")
cat("========================================\n\n")

cat(
    "Potentially informative columns:\n",
    paste(possible_cols, collapse = "\n"),
    "\n\n"
)

if (length(possible_cols) > 0) {

    print(
        meta[
            ,
            possible_cols,
            drop = FALSE
        ]
    )
}

# ------------------------------------------------------------
# 11. Save simple sample audit
# ------------------------------------------------------------

audit <- meta[
    ,
    unique(
        c(
            intersect(
                c(
                    "title",
                    "geo_accession",
                    "source_name_ch1"
                ),
                names(meta)
            ),
            possible_cols
        )
    ),
    drop = FALSE
]

write.csv(
    audit,
    file.path(
        result_dir,
        "09a_GSE84161_sample_audit.csv"
    ),
    row.names = FALSE
)

cat("\n========================================\n")
cat("DOWNLOAD/AUDIT COMPLETE\n")
cat("========================================\n\n")

cat(
    "Saved GEO object:\n",
    file.path(
        data_dir,
        "GSE84161_GEO_expressionSet.rds"
    ),
    "\n\n"
)

cat(
    "Saved sample metadata:\n",
    file.path(
        result_dir,
        "09a_GSE84161_complete_sample_metadata.csv"
    ),
    "\n\n"
)

cat(
    "STOP HERE.\n",
    "Sample design must be verified before differential analysis.\n"
)
