# ============================================================
# 10b_GSE116220_MAP3K8_genetic_DE.R
#
# Independent genetic validation of MAP3K8/TPL2
#
# Dataset:
#   GSE116220
#
# Biological comparison:
#   Map3k8[D270A] kinase-inactive BMDMs
#                 versus
#   WT BMDMs
#
# Time points:
#   unstimulated
#   0.5 h LPS
#   1 h LPS
#   2 h LPS
#
# IMPORTANT:
# - These are independent biological replicates, NOT paired.
# - GEO processed matrix is already normalized.
# - Do NOT run DESeq2 on this processed matrix.
# - Do NOT log-transform again.
# ============================================================

rm(list = ls())

# ============================================================
# 1. PACKAGE
# ============================================================

if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required.")
}

# ============================================================
# 2. INPUT / OUTPUT
# ============================================================

expr_file <-
    "data/external/GSE116220/GSE116220_GEO_RSEM_normalised_data.txt.gz"

outdir <-
    "results/10_MAP3K8_genetic_validation"

figdir <-
    "figures/10_MAP3K8_genetic_validation"

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

if (!file.exists(expr_file)) {
    stop("Missing GSE116220 processed expression file.")
}

cat("\n========================================\n")
cat("10b — MAP3K8 GENETIC VALIDATION\n")
cat("========================================\n\n")

# ============================================================
# 3. LOAD PROCESSED EXPRESSION
# ============================================================

dat <- read.delim(
    gzfile(expr_file),
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (!"ID" %in% names(dat)) {
    stop("Expected ID column is absent.")
}

gene_symbol <- as.character(dat$ID)

expr_all <- as.matrix(
    dat[
        ,
        setdiff(
            names(dat),
            "ID"
        ),
        drop = FALSE
    ]
)

storage.mode(expr_all) <- "numeric"

rownames(expr_all) <- gene_symbol

cat(
    "Genes:",
    nrow(expr_all),
    "\n"
)

cat(
    "All samples:",
    ncol(expr_all),
    "\n"
)

cat(
    "Expression range:",
    paste(
        round(
            range(
                expr_all,
                na.rm = TRUE
            ),
            3
        ),
        collapse = " to "
    ),
    "\n\n"
)

# ============================================================
# 4. SELECT ONLY WT vs Map3k8[D270A] EXPERIMENT
# ============================================================

sample_names <- colnames(expr_all)

keep_samples <- grepl(
    "^(unstim|0\\.5hr|1hr|2hr)_(WT|TPL)[0-9]+$",
    sample_names
)

expr <- expr_all[
    ,
    keep_samples,
    drop = FALSE
]

samples <- colnames(expr)

if (ncol(expr) != 46) {
    stop(
        paste(
            "Expected 46 WT/TPL samples; found",
            ncol(expr)
        )
    )
}

# ============================================================
# 5. PARSE SAMPLE DESIGN
# ============================================================

time <- sub(
    "_.*$",
    "",
    samples
)

genotype_code <- ifelse(
    grepl("_WT[0-9]+$", samples),
    "WT",
    "D270A"
)

replicate <- sub(
    "^.*_(WT|TPL)",
    "",
    samples
)

time <- factor(
    time,
    levels = c(
        "unstim",
        "0.5hr",
        "1hr",
        "2hr"
    )
)

genotype <- factor(
    genotype_code,
    levels = c(
        "WT",
        "D270A"
    )
)

manifest <- data.frame(
    sample = samples,
    genotype = genotype,
    time = time,
    replicate = replicate,
    stringsAsFactors = FALSE
)

write.csv(
    manifest,
    file.path(
        outdir,
        "10b_GSE116220_WT_D270A_manifest.csv"
    ),
    row.names = FALSE
)

cat("Genotype x time sample counts:\n\n")

print(
    table(
        genotype,
        time
    )
)

# ============================================================
# 6. CHECK DUPLICATED GENE SYMBOLS
# ============================================================

cat(
    "\nDuplicated gene symbols:",
    sum(
        duplicated(
            rownames(expr)
        )
    ),
    "\n"
)

if (anyDuplicated(rownames(expr))) {

    expr <- limma::avereps(
        expr,
        ID = rownames(expr)
    )
}

cat(
    "Unique genes used:",
    nrow(expr),
    "\n"
)

# ============================================================
# 7. QC — SAMPLE DISTRIBUTIONS
# ============================================================

png(
    file.path(
        figdir,
        "10B1_GSE116220_expression_boxplot.png"
    ),
    width = 2000,
    height = 1100,
    res = 160
)

boxplot(
    expr,
    las = 2,
    outline = FALSE,
    cex.axis = 0.55,
    ylab = "GEO processed normalized expression",
    main = "GSE116220 WT vs Map3k8[D270A]"
)

dev.off()

# ============================================================
# 8. PCA
# ============================================================

gene_var <- apply(
    expr,
    1,
    var,
    na.rm = TRUE
)

top_var <- order(
    gene_var,
    decreasing = TRUE
)

top_var <- top_var[
    seq_len(
        min(
            2000,
            length(top_var)
        )
    )
]

pca <- prcomp(
    t(
        expr[
            top_var,
            ,
            drop = FALSE
        ]
    ),
    center = TRUE,
    scale. = FALSE
)

png(
    file.path(
        figdir,
        "10B2_GSE116220_PCA.png"
    ),
    width = 1500,
    height = 1200,
    res = 160
)

plot(
    pca$x[, 1],
    pca$x[, 2],
    pch = ifelse(
        genotype == "WT",
        16,
        1
    ),
    xlab = "PC1",
    ylab = "PC2",
    main = "GSE116220 WT vs Map3k8[D270A] PCA"
)

text(
    pca$x[, 1],
    pca$x[, 2],
    labels = paste0(
        genotype,
        "_",
        time
    ),
    pos = 3,
    cex = 0.55
)

dev.off()

# ============================================================
# 9. GROUP FACTOR
# ============================================================

group <- interaction(
    genotype,
    time,
    sep = "_"
)

group <- factor(
    group,
    levels = c(
        "WT_unstim",
        "D270A_unstim",
        "WT_0.5hr",
        "D270A_0.5hr",
        "WT_1hr",
        "D270A_1hr",
        "WT_2hr",
        "D270A_2hr"
    )
)

if (anyNA(group)) {
    stop("Failed to construct genotype/time groups.")
}

cat("\nGroup counts:\n\n")
print(table(group))

# ============================================================
# 10. DESIGN
# ============================================================

design <- model.matrix(
    ~ 0 + group
)

colnames(design) <- levels(group)

rownames(design) <- samples

cat(
    "\nDesign rank:",
    qr(design)$rank,
    "/",
    ncol(design),
    "\n"
)

if (qr(design)$rank != ncol(design)) {
    stop("Design matrix is not full rank.")
}

# ============================================================
# 11. CONTRASTS
#
# Positive logFC:
#   higher in kinase-inactive D270A
#
# Negative logFC:
#   lower in kinase-inactive D270A
#
# Therefore:
#
# negative = normally promoted by TPL2 kinase activity
# positive = normally suppressed by TPL2 kinase activity
# ============================================================

contrasts <- limma::makeContrasts(

    D270A_vs_WT_unstim =
        D270A_unstim -
        WT_unstim,

    D270A_vs_WT_0.5hr =
        D270A_0.5hr -
        WT_0.5hr,

    D270A_vs_WT_1hr =
        D270A_1hr -
        WT_1hr,

    D270A_vs_WT_2hr =
        D270A_2hr -
        WT_2hr,

    Interaction_0.5hr =
        (
            D270A_0.5hr -
            D270A_unstim
        ) -
        (
            WT_0.5hr -
            WT_unstim
        ),

    Interaction_1hr =
        (
            D270A_1hr -
            D270A_unstim
        ) -
        (
            WT_1hr -
            WT_unstim
        ),

    Interaction_2hr =
        (
            D270A_2hr -
            D270A_unstim
        ) -
        (
            WT_2hr -
            WT_unstim
        ),

    levels = design
)

# ============================================================
# 12. LIMMA
# ============================================================

fit <- limma::lmFit(
    expr,
    design
)

fit <- limma::contrasts.fit(
    fit,
    contrasts
)

fit <- limma::eBayes(
    fit,
    trend = TRUE,
    robust = TRUE
)

# ============================================================
# 13. SAVE ALL CONTRASTS
# ============================================================

results <- list()

for (cn in colnames(contrasts)) {

    tt <- limma::topTable(
        fit,
        coef = cn,
        number = Inf,
        sort.by = "P",
        adjust.method = "BH"
    )

    tt$SYMBOL <- rownames(tt)

    names(tt)[
        names(tt) == "logFC"
    ] <- "Effect"

    names(tt)[
        names(tt) == "adj.P.Val"
    ] <- "FDR"

    tt$Direction_in_D270A <- ifelse(
        tt$Effect > 0,
        "UP",
        ifelse(
            tt$Effect < 0,
            "DOWN",
            "UNCHANGED"
        )
    )

    tt$FDR_lt_0.05 <-
        tt$FDR < 0.05

    tt$Abs_effect_ge_1 <-
        abs(tt$Effect) >= 1

    tt$FDR05_effect1 <-
        tt$FDR < 0.05 &
        abs(tt$Effect) >= 1

    tt <- tt[
        ,
        c(
            "SYMBOL",
            "Effect",
            "AveExpr",
            "t",
            "P.Value",
            "FDR",
            "B",
            "Direction_in_D270A",
            "FDR_lt_0.05",
            "Abs_effect_ge_1",
            "FDR05_effect1"
        )
    ]

    write.csv(
        tt,
        file.path(
            outdir,
            paste0(
                "10b_",
                cn,
                ".csv"
            )
        ),
        row.names = FALSE
    )

    results[[cn]] <- tt
}

# ============================================================
# 14. SUMMARY
# ============================================================

summary_table <- do.call(
    rbind,
    lapply(
        names(results),
        function(cn) {

            x <- results[[cn]]

            data.frame(
                Contrast = cn,
                Genes = nrow(x),

                FDR_lt_0.05 =
                    sum(
                        x$FDR < 0.05,
                        na.rm = TRUE
                    ),

                FDR_UP_in_D270A =
                    sum(
                        x$FDR < 0.05 &
                        x$Effect > 0,
                        na.rm = TRUE
                    ),

                FDR_DOWN_in_D270A =
                    sum(
                        x$FDR < 0.05 &
                        x$Effect < 0,
                        na.rm = TRUE
                    ),

                FDR05_absEffect_ge1 =
                    sum(
                        x$FDR < 0.05 &
                        abs(x$Effect) >= 1,
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
        "10b_genetic_contrast_summary.csv"
    ),
    row.names = FALSE
)

cat("\n========================================\n")
cat("GENETIC DIFFERENTIAL SUMMARY\n")
cat("========================================\n\n")

print(
    summary_table,
    row.names = FALSE
)

# ============================================================
# 15. PUBLISHED POSITIVE-CONTROL GENES
# ============================================================

control_genes <- c(
    "Fos",
    "Egr1",
    "Egr3",
    "Nr4a1",
    "Ifnb1",
    "Il10"
)

time_contrasts <- c(
    "D270A_vs_WT_unstim",
    "D270A_vs_WT_0.5hr",
    "D270A_vs_WT_1hr",
    "D270A_vs_WT_2hr"
)

control_table <- NULL

for (cn in time_contrasts) {

    x <- results[[cn]]

    hit <- x[
        x$SYMBOL %in%
            control_genes,
        ,
        drop = FALSE
    ]

    if (nrow(hit) > 0) {

        hit$Contrast <- cn

        control_table <- rbind(
            control_table,
            hit[
                ,
                c(
                    "Contrast",
                    "SYMBOL",
                    "Effect",
                    "P.Value",
                    "FDR",
                    "Direction_in_D270A"
                )
            ]
        )
    }
}

write.csv(
    control_table,
    file.path(
        outdir,
        "10b_published_positive_control_genes.csv"
    ),
    row.names = FALSE
)

cat("\n========================================\n")
cat("PUBLISHED POSITIVE-CONTROL GENES\n")
cat("========================================\n\n")

print(
    control_table,
    row.names = FALSE
)

# ============================================================
# 16. MAP3K8 / Map3k8 ITSELF
# ============================================================

map3k8_rows <- NULL

for (cn in time_contrasts) {

    x <- results[[cn]]

    hit <- x[
        tolower(x$SYMBOL) ==
            "map3k8",
        ,
        drop = FALSE
    ]

    if (nrow(hit) > 0) {

        hit$Contrast <- cn

        map3k8_rows <- rbind(
            map3k8_rows,
            hit[
                ,
                c(
                    "Contrast",
                    "SYMBOL",
                    "Effect",
                    "P.Value",
                    "FDR",
                    "Direction_in_D270A"
                )
            ]
        )
    }
}

cat("\n========================================\n")
cat("Map3k8 RNA ITSELF\n")
cat("========================================\n\n")

print(
    map3k8_rows,
    row.names = FALSE
)

write.csv(
    map3k8_rows,
    file.path(
        outdir,
        "10b_Map3k8_gene_results.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 17. TOP GENETIC EFFECTS AT 2 HOURS
# ============================================================

x2 <- results[["D270A_vs_WT_2hr"]]

cat("\n========================================\n")
cat("TOP D270A vs WT GENES — 2 HOURS\n")
cat("========================================\n\n")

print(
    head(
        x2[
            ,
            c(
                "SYMBOL",
                "Effect",
                "P.Value",
                "FDR",
                "Direction_in_D270A"
            )
        ],
        20
    ),
    row.names = FALSE
)

cat("\n========================================\n")
cat("10b COMPLETE\n")
cat("========================================\n\n")

cat(
    "Negative D270A-vs-WT effect:\n",
    "gene normally promoted by TPL2 catalytic activity.\n\n",

    "Positive D270A-vs-WT effect:\n",
    "gene normally suppressed by TPL2 catalytic activity.\n\n",

    "Next step: compare this genetic signature with the\n",
    "independent human TPL2-inhibitor signature from GSE84161.\n"
)
