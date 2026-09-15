# ============================================================
# 10a_download_audit_GSE116220.R
#
# Independent genetic validation of MAP3K8/TPL2
#
# Dataset:
#   GSE116220
#
# Purpose of this step:
#   - download GEO metadata
#   - download processed RSEM-normalized expression
#   - identify WT and kinase-inactive Map3k8 samples
#   - audit time points and replicates
#
# NO differential analysis is performed here.
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
        "BiocManager",
        repos = "https://cloud.r-project.org"
    )
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

data_dir <-
    "data/external/GSE116220"

result_dir <-
    "results/10_MAP3K8_genetic_validation"

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
cat("10a — GSE116220 DOWNLOAD + AUDIT\n")
cat("========================================\n\n")

# ------------------------------------------------------------
# 3. GEO series metadata
# ------------------------------------------------------------

gse_list <-
    GEOquery::getGEO(
        "GSE116220",
        GSEMatrix = TRUE,
        getGPL = FALSE,
        destdir = data_dir
    )

cat(
    "ExpressionSet objects:",
    length(gse_list),
    "\n"
)

if (length(gse_list) < 1) {
    stop("No GEO ExpressionSet returned.")
}

eset <-
    gse_list[[1]]

meta <-
    Biobase::pData(eset)

cat(
    "GEO samples:",
    nrow(meta),
    "\n\n"
)

# ------------------------------------------------------------
# 4. Save complete GEO metadata
# ------------------------------------------------------------

meta$GEO_sample <-
    rownames(meta)

write.csv(
    meta,
    file.path(
        result_dir,
        "10a_GSE116220_complete_sample_metadata.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 5. Inspect informative metadata
# ------------------------------------------------------------

possible_cols <-
    grep(
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

cat("Informative metadata columns:\n\n")

cat(
    paste(
        possible_cols,
        collapse = "\n"
    ),
    "\n\n"
)

# ------------------------------------------------------------
# 6. Candidate WT / MAP3K8 kinase-inactive samples
# ------------------------------------------------------------

title <-
    as.character(
        meta$title
    )

candidate_idx <-
    grepl(
        "WT|TPL",
        title,
        ignore.case = FALSE
    )

candidate_meta <-
    meta[
        candidate_idx,
        ,
        drop = FALSE
    ]

cat(
    "WT/TPL candidate samples:",
    nrow(candidate_meta),
    "\n\n"
)

if (nrow(candidate_meta) > 0) {

    print(
        candidate_meta[
            ,
            intersect(
                c(
                    "title",
                    "geo_accession",
                    "source_name_ch1",
                    "characteristics_ch1",
                    "characteristics_ch1.1",
                    "characteristics_ch1.2",
                    "description"
                ),
                names(candidate_meta)
            ),
            drop = FALSE
        ]
    )
}

write.csv(
    candidate_meta,
    file.path(
        result_dir,
        "10a_GSE116220_WT_TPL_candidate_metadata.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 7. Download supplementary processed expression data
# ------------------------------------------------------------

cat("\n========================================\n")
cat("DOWNLOADING SUPPLEMENTARY FILES\n")
cat("========================================\n\n")

GEOquery::getGEOSuppFiles(
    "GSE116220",
    makeDirectory = FALSE,
    baseDir = data_dir
)

files <-
    list.files(
        data_dir,
        full.names = TRUE
    )

cat("Files present:\n\n")

print(
    data.frame(
        file = basename(files),
        size_MB =
            round(
                file.info(files)$size /
                1024^2,
                2
            )
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 8. Locate RSEM normalized file
# ------------------------------------------------------------

norm_file <-
    files[
        grepl(
            "RSEM.*normal",
            basename(files),
            ignore.case = TRUE
        )
    ]

if (length(norm_file) != 1) {

    stop(
        paste(
            "Expected exactly one RSEM normalized file.",
            "Found:",
            length(norm_file)
        )
    )
}

cat(
    "\nNormalized expression file:\n",
    norm_file,
    "\n"
)

# ------------------------------------------------------------
# 9. Read processed expression table
# ------------------------------------------------------------

expr_df <-
    read.delim(
        gzfile(norm_file),
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

cat("\n========================================\n")
cat("PROCESSED EXPRESSION TABLE\n")
cat("========================================\n\n")

cat(
    "Rows:",
    nrow(expr_df),
    "\n"
)

cat(
    "Columns:",
    ncol(expr_df),
    "\n\n"
)

cat("First 20 column names:\n\n")

print(
    head(
        names(expr_df),
        20
    )
)

cat("\nFirst few rows/columns:\n\n")

print(
    expr_df[
        seq_len(
            min(
                5,
                nrow(expr_df)
            )
        ),
        seq_len(
            min(
                8,
                ncol(expr_df)
            )
        ),
        drop = FALSE
    ]
)

# ------------------------------------------------------------
# 10. Save column audit
# ------------------------------------------------------------

write.csv(
    data.frame(
        column_number =
            seq_along(
                names(expr_df)
            ),
        column_name =
            names(expr_df)
    ),
    file.path(
        result_dir,
        "10a_GSE116220_expression_columns.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# 11. Basic numeric audit
# ------------------------------------------------------------

numeric_cols <-
    vapply(
        expr_df,
        is.numeric,
        logical(1)
    )

cat(
    "\nNumeric columns:",
    sum(numeric_cols),
    "\n"
)

if (any(numeric_cols)) {

    vals <-
        as.matrix(
            expr_df[
                ,
                numeric_cols,
                drop = FALSE
            ]
        )

    cat(
        "Numeric value range:",
        paste(
            signif(
                range(
                    vals,
                    na.rm = TRUE
                ),
                5
            ),
            collapse = " to "
        ),
        "\n"
    )
}

cat("\n========================================\n")
cat("10a COMPLETE\n")
cat("========================================\n\n")

cat(
    "STOP HERE.\n",
    "The exact WT vs kinase-inactive Map3k8 design\n",
    "must be verified before differential analysis.\n"
)
