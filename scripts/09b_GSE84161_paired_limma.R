# ============================================================
# 09b_GSE84161_paired_limma.R
#
# GSE84161
# Human primary monocytes
#
# Goal:
# Quantify transcriptional effects of:
#   1. acute LPS
#   2. TPL2/MAP3K8 inhibition during LPS
#   3. MEK inhibition during LPS
#
# IMPORTANT:
# This step analyzes GSE84161 independently.
# It does NOT yet force comparison with GSE85246.
# ============================================================

rm(list = ls())

# ============================================================
# 1. PACKAGES
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
        "BiocManager",
        repos = "https://cloud.r-project.org"
    )
}

needed <- c(
    "limma",
    "Biobase",
    "AnnotationDbi",
    "hgu133plus2.db"
)

for (pkg in needed) {

    if (!requireNamespace(pkg, quietly = TRUE)) {

        BiocManager::install(
            pkg,
            ask = FALSE,
            update = FALSE
        )
    }
}

# ============================================================
# 2. PATHS
# ============================================================

input_file <-
    "data/external/GSE84161/GSE84161_GEO_expressionSet.rds"

outdir <-
    "results/09_MAP3K8_validation"

figdir <-
    "figures/09_MAP3K8_validation"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    figdir,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(input_file)) {
    stop("Missing GSE84161 ExpressionSet. Run 09a first.")
}

# ============================================================
# 3. LOAD DATA
# ============================================================

eset <- readRDS(input_file)

expr_raw <- Biobase::exprs(eset)
meta <- Biobase::pData(eset)

meta <- meta[
    colnames(expr_raw),
    ,
    drop = FALSE
]

stopifnot(
    identical(
        rownames(meta),
        colnames(expr_raw)
    )
)

cat("\n========================================\n")
cat("09b — GSE84161 PAIRED LIMMA ANALYSIS\n")
cat("========================================\n\n")

cat("Features:", nrow(expr_raw), "\n")
cat("Samples:", ncol(expr_raw), "\n")
cat(
    "Raw expression range:",
    paste(
        round(range(expr_raw, na.rm = TRUE), 3),
        collapse = " to "
    ),
    "\n\n"
)

# ============================================================
# 4. DEFINE CONDITIONS
# ============================================================

treatment <-
    as.character(
        meta[["treatment name:ch1"]]
    )

condition_map <- c(
    "Untreated" = "Untreated",
    "Tpl-2 SMI" = "TPL2i",
    "MEK SMI" = "MEKi",
    "LPS" = "LPS",
    "Tpl-2 SMI + LPS" = "TPL2i_LPS",
    "MEK SMI + LPS" = "MEKi_LPS"
)

condition <-
    unname(
        condition_map[treatment]
    )

if (anyNA(condition)) {

    stop(
        "Unexpected treatment names detected."
    )
}

condition <- factor(
    condition,
    levels = c(
        "Untreated",
        "TPL2i",
        "MEKi",
        "LPS",
        "TPL2i_LPS",
        "MEKi_LPS"
    )
)

# ============================================================
# 5. DEFINE DONORS
# ============================================================

donor_text <-
    as.character(
        meta[["characteristics_ch1.4"]]
    )

donor <-
    sub(
        "^replicate number:[[:space:]]*",
        "",
        donor_text
    )

donor <- factor(donor)

if (length(levels(donor)) != 5) {
    stop("Expected exactly 5 donors.")
}

sample_info <- data.frame(
    GSM = rownames(meta),
    title = meta$title,
    donor = donor,
    treatment = treatment,
    condition = condition,
    stringsAsFactors = FALSE
)

write.csv(
    sample_info,
    file.path(
        outdir,
        "09b_GSE84161_analysis_manifest.csv"
    ),
    row.names = FALSE
)

cat("Condition counts:\n\n")
print(table(condition))

cat("\nDonor x condition table:\n\n")
print(
    table(
        donor,
        condition
    )
)

if (!all(
    table(donor, condition) == 1
)) {

    stop(
        "Design is not complete: each donor should have one sample per condition."
    )
}

# ============================================================
# 6. LOG2 TRANSFORMATION
# ============================================================

if (
    max(
        expr_raw,
        na.rm = TRUE
    ) > 100
) {

    cat(
        "\nExpression is not log-scale.",
        "Applying log2 transformation.\n"
    )

    expr_log2 <- log2(expr_raw)

} else {

    cat(
        "\nExpression appears already log-scaled.",
        "No transformation applied.\n"
    )

    expr_log2 <- expr_raw
}

cat(
    "Log2 expression range:",
    paste(
        round(
            range(
                expr_log2,
                na.rm = TRUE
            ),
            3
        ),
        collapse = " to "
    ),
    "\n"
)

# ============================================================
# 7. PRE-ANALYSIS QC
# ============================================================

png(
    file.path(
        figdir,
        "09B1_GSE84161_log2_expression_boxplot.png"
    ),
    width = 1800,
    height = 1000,
    res = 160
)

boxplot(
    expr_log2,
    las = 2,
    outline = FALSE,
    main = "GSE84161 log2 expression distributions",
    ylab = "log2 expression"
)

dev.off()

# ============================================================
# 8. AFFYMETRIX PROBE -> HUMAN GENE SYMBOL
# ============================================================

probe_ids <-
    rownames(expr_log2)

symbols <-
    AnnotationDbi::mapIds(
        hgu133plus2.db::hgu133plus2.db,
        keys = probe_ids,
        keytype = "PROBEID",
        column = "SYMBOL",
        multiVals = "first"
    )

symbols <-
    unname(
        symbols[probe_ids]
    )

keep <-
    !is.na(symbols) &
    symbols != ""

cat(
    "\nMapped probes:",
    sum(keep),
    "of",
    length(probe_ids),
    "\n"
)

expr_mapped <-
    expr_log2[
        keep,
        ,
        drop = FALSE
    ]

mapped_symbols <-
    symbols[keep]

# ============================================================
# 9. COLLAPSE MULTIPLE PROBES PER GENE
#
# Condition-blind averaging.
# No significance or direction is used to select probes.
# ============================================================

expr_gene <-
    limma::avereps(
        expr_mapped,
        ID = mapped_symbols
    )

cat(
    "Unique gene symbols after probe collapsing:",
    nrow(expr_gene),
    "\n"
)

write.csv(
    data.frame(
        SYMBOL = rownames(expr_gene),
        expr_gene,
        check.names = FALSE
    ),
    file.path(
        outdir,
        "09b_GSE84161_gene_level_log2_expression.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 10. PCA QC
# ============================================================

pca <-
    prcomp(
        t(expr_gene),
        center = TRUE,
        scale. = FALSE
    )

png(
    file.path(
        figdir,
        "09B2_GSE84161_PCA.png"
    ),
    width = 1400,
    height = 1100,
    res = 160
)

plot(
    pca$x[, 1],
    pca$x[, 2],
    pch = 19,
    xlab = "PC1",
    ylab = "PC2",
    main = "GSE84161 gene-level PCA"
)

text(
    pca$x[, 1],
    pca$x[, 2],
    labels = paste0(
        donor,
        "_",
        condition
    ),
    pos = 3,
    cex = 0.65
)

dev.off()

# ============================================================
# 11. PAIRED LIMMA DESIGN
#
# Every donor contributes every condition.
# Donor is therefore explicitly modeled.
# ============================================================

design <-
    model.matrix(
        ~ 0 + condition + donor
    )

colnames(design) <-
    sub(
        "^condition",
        "",
        colnames(design)
    )

rownames(design) <-
    colnames(expr_gene)

cat("\nDesign matrix rank:")
cat(
    qr(design)$rank,
    "/",
    ncol(design),
    "\n"
)

if (
    qr(design)$rank !=
    ncol(design)
) {

    stop(
        "Design matrix is not full rank."
    )
}

# ============================================================
# 12. CONTRASTS
# ============================================================

contrast_matrix <-
    limma::makeContrasts(

        LPS_vs_Untreated =
            LPS - Untreated,

        TPL2i_LPS_vs_LPS =
            TPL2i_LPS - LPS,

        MEKi_LPS_vs_LPS =
            MEKi_LPS - LPS,

        TPL2i_vs_Untreated =
            TPL2i - Untreated,

        MEKi_vs_Untreated =
            MEKi - Untreated,

        levels = design
    )

# ============================================================
# 13. LIMMA MODEL
# ============================================================

fit <-
    limma::lmFit(
        expr_gene,
        design
    )

fit2 <-
    limma::contrasts.fit(
        fit,
        contrast_matrix
    )

fit2 <-
    limma::eBayes(
        fit2,
        trend = TRUE,
        robust = TRUE
    )

# ============================================================
# 14. SAVE EACH CONTRAST
# ============================================================

contrast_names <-
    colnames(
        contrast_matrix
    )

all_results <- list()

for (cn in contrast_names) {

    tt <-
        limma::topTable(
            fit2,
            coef = cn,
            number = Inf,
            sort.by = "P",
            adjust.method = "BH"
        )

    tt$SYMBOL <-
        rownames(tt)

    names(tt)[
        names(tt) == "logFC"
    ] <- "log2FC"

    names(tt)[
        names(tt) == "adj.P.Val"
    ] <- "FDR"

    tt$Direction <-
        ifelse(
            tt$log2FC > 0,
            "UP",
            ifelse(
                tt$log2FC < 0,
                "DOWN",
                "UNCHANGED"
            )
        )

    tt$FDR_lt_0.05 <-
        tt$FDR < 0.05

    tt$Nominal_P_lt_0.05 <-
        tt$P.Value < 0.05

    tt <-
        tt[
            ,
            c(
                "SYMBOL",
                "log2FC",
                "AveExpr",
                "t",
                "P.Value",
                "FDR",
                "B",
                "Direction",
                "FDR_lt_0.05",
                "Nominal_P_lt_0.05"
            )
        ]

    write.csv(
        tt,
        file.path(
            outdir,
            paste0(
                "09b_",
                cn,
                ".csv"
            )
        ),
        row.names = FALSE
    )

    all_results[[cn]] <- tt
}

# ============================================================
# 15. SUMMARY COUNTS
# ============================================================

summary_table <-
    do.call(
        rbind,
        lapply(
            names(all_results),
            function(cn) {

                x <- all_results[[cn]]

                data.frame(
                    Contrast = cn,

                    Genes_tested =
                        nrow(x),

                    FDR_lt_0.05 =
                        sum(
                            x$FDR < 0.05,
                            na.rm = TRUE
                        ),

                    FDR_UP =
                        sum(
                            x$FDR < 0.05 &
                            x$log2FC > 0,
                            na.rm = TRUE
                        ),

                    FDR_DOWN =
                        sum(
                            x$FDR < 0.05 &
                            x$log2FC < 0,
                            na.rm = TRUE
                        )
                )
            }
        )
    )

write.csv(
    summary_table,
    file.path(
        outdir,
        "09b_contrast_summary.csv"
    ),
    row.names = FALSE
)

cat("\n========================================\n")
cat("DIFFERENTIAL EXPRESSION SUMMARY\n")
cat("========================================\n\n")

print(
    summary_table,
    row.names = FALSE
)

# ============================================================
# 16. MAP3K8 RNA ITSELF
#
# Important:
# TPL2i inhibits MAP3K8 kinase ACTIVITY.
# MAP3K8 mRNA does NOT need to decrease.
# ============================================================

cat("\n========================================\n")
cat("MAP3K8 GENE-LEVEL RESULTS\n")
cat("========================================\n\n")

map3k8_table <- NULL

for (cn in contrast_names) {

    x <-
        all_results[[cn]]

    hit <-
        x[
            x$SYMBOL == "MAP3K8",
            ,
            drop = FALSE
        ]

    if (nrow(hit) > 0) {

        hit$Contrast <- cn

        map3k8_table <-
            rbind(
                map3k8_table,
                hit[
                    ,
                    c(
                        "Contrast",
                        "SYMBOL",
                        "log2FC",
                        "P.Value",
                        "FDR",
                        "Direction"
                    )
                ]
            )
    }
}

print(
    map3k8_table,
    row.names = FALSE
)

write.csv(
    map3k8_table,
    file.path(
        outdir,
        "09b_MAP3K8_gene_results.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 17. TOP TPL2-INHIBITOR EFFECTS DURING LPS
# ============================================================

tpl2_lps <- all_results[["TPL2i_LPS_vs_LPS"]]

cat("\n========================================\n")
cat("TOP TPL2i + LPS vs LPS GENES\n")
cat("========================================\n\n")

print(
    head(
        tpl2_lps[
            ,
            c(
                "SYMBOL",
                "log2FC",
                "P.Value",
                "FDR",
                "Direction"
            )
        ],
        20
    ),
    row.names = FALSE
)

cat("\n========================================\n")
cat("09b COMPLETE\n")
cat("========================================\n\n")

cat(
    "Next step will compare these independent\n",
    "TPL2 perturbation effects with the GSE85246\n",
    "Day-6 persistent memory signature.\n"
)
