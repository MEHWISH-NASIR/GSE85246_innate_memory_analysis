# ============================================================
# 09c_memory_vs_TPL2_validation.R
#
# Cross-study validation:
#
# GSE85246:
#   Human Day-6 persistent LPS-memory transcriptome
#
# GSE84161:
#   Human monocytes
#   LPS +/- TPL2/MAP3K8 inhibitor
#   LPS +/- MEK inhibitor
#
# CENTRAL QUESTION:
# Are persistent Day-6 memory genes systematically affected
# when TPL2/MAP3K8 activity is inhibited?
#
# IMPORTANT:
# We do not require the studies to agree.
# Reversal, concordance, or no relationship are all possible.
# ============================================================

rm(list = ls())

# ============================================================
# 1. PACKAGES
# ============================================================

if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required.")
}

# ============================================================
# 2. INPUTS / OUTPUTS
# ============================================================

memory_file <-
    "results/03_memory_kinases/03_Day6_LPS_vs_RPMI_all_genes.csv"

tpl2_file <-
    "results/09_MAP3K8_validation/09b_TPL2i_LPS_vs_LPS.csv"

mek_file <-
    "results/09_MAP3K8_validation/09b_MEKi_LPS_vs_LPS.csv"

acute_lps_file <-
    "results/09_MAP3K8_validation/09b_LPS_vs_Untreated.csv"

expr_file <-
    "results/09_MAP3K8_validation/09b_GSE84161_gene_level_log2_expression.csv"

manifest_file <-
    "results/09_MAP3K8_validation/09b_GSE84161_analysis_manifest.csv"

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

required_files <- c(
    memory_file,
    tpl2_file,
    mek_file,
    acute_lps_file,
    expr_file,
    manifest_file
)

missing_files <-
    required_files[
        !file.exists(required_files)
    ]

if (length(missing_files) > 0) {
    stop(
        paste(
            "Missing required files:",
            paste(missing_files, collapse = "\n")
        )
    )
}

cat("\n========================================\n")
cat("09c — MEMORY vs TPL2 VALIDATION\n")
cat("========================================\n\n")

# ============================================================
# 3. LOAD GSE85246 DAY-6 MEMORY RESULTS
# ============================================================

memory <-
    read.csv(
        memory_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

required_memory_cols <- c(
    "SYMBOL",
    "log2FC_Day6_LPS_vs_RPMI",
    "P.Value",
    "FDR_genome",
    "AveExpr"
)

if (!all(
    required_memory_cols %in%
        names(memory)
)) {
    stop(
        "Unexpected columns in GSE85246 Day-6 result."
    )
}

memory <-
    memory[
        !is.na(memory$SYMBOL) &
        memory$SYMBOL != "",
        ,
        drop = FALSE
    ]

# ============================================================
# 4. COLLAPSE DUPLICATED SYMBOLS
#
# Some ENSEMBL IDs may map to the same gene symbol.
#
# We choose the row with highest average expression.
# This is condition-blind:
# significance or fold-change direction is NOT used to choose
# the representative row.
# ============================================================

memory <-
    memory[
        order(
            memory$SYMBOL,
            -memory$AveExpr
        ),
        ,
        drop = FALSE
    ]

memory <-
    memory[
        !duplicated(memory$SYMBOL),
        ,
        drop = FALSE
    ]

names(memory)[
    names(memory) ==
        "log2FC_Day6_LPS_vs_RPMI"
] <- "Memory_log2FC"

names(memory)[
    names(memory) ==
        "P.Value"
] <- "Memory_P"

names(memory)[
    names(memory) ==
        "FDR_genome"
] <- "Memory_FDR"

cat(
    "Unique GSE85246 genes:",
    nrow(memory),
    "\n"
)

# ============================================================
# 5. DEFINE MEMORY SIGNATURES
#
# PRIMARY:
# genome-wide FDR < 0.05 AND |log2FC| >= 1
#
# SENSITIVITY:
# genome-wide FDR < 0.05 regardless of effect size
# ============================================================

memory$Memory_FDR05 <-
    memory$Memory_FDR < 0.05

memory$Memory_robust <-
    memory$Memory_FDR < 0.05 &
    abs(memory$Memory_log2FC) >= 1

memory$Memory_direction <-
    ifelse(
        memory$Memory_log2FC > 0,
        "UP",
        ifelse(
            memory$Memory_log2FC < 0,
            "DOWN",
            "UNCHANGED"
        )
    )

memory$Robust_memory_class <-
    ifelse(
        memory$Memory_robust &
        memory$Memory_log2FC > 0,
        "MEMORY_UP",
        ifelse(
            memory$Memory_robust &
            memory$Memory_log2FC < 0,
            "MEMORY_DOWN",
            "NOT_ROBUST_MEMORY"
        )
    )

cat(
    "Day-6 genome-wide FDR<0.05 genes:",
    sum(memory$Memory_FDR05, na.rm = TRUE),
    "\n"
)

cat(
    "Robust memory genes (FDR<0.05 & |log2FC|>=1):",
    sum(memory$Memory_robust, na.rm = TRUE),
    "\n"
)

cat(
    "Robust memory UP:",
    sum(
        memory$Robust_memory_class ==
            "MEMORY_UP"
    ),
    "\n"
)

cat(
    "Robust memory DOWN:",
    sum(
        memory$Robust_memory_class ==
            "MEMORY_DOWN"
    ),
    "\n\n"
)

write.csv(
    memory[
        memory$Memory_robust,
        ,
        drop = FALSE
    ],
    file.path(
        outdir,
        "09c_GSE85246_robust_Day6_memory_signature.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 6. LOAD GSE84161 DIFFERENTIAL RESULTS
# ============================================================

load_contrast <- function(
    file,
    prefix
) {

    x <-
        read.csv(
            file,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )

    required <- c(
        "SYMBOL",
        "log2FC",
        "P.Value",
        "FDR"
    )

    if (!all(required %in% names(x))) {
        stop(
            paste(
                "Unexpected columns in",
                file
            )
        )
    }

    x <-
        x[
            ,
            required,
            drop = FALSE
        ]

    names(x)[2:4] <-
        paste0(
            prefix,
            c(
                "_log2FC",
                "_P",
                "_FDR"
            )
        )

    x
}

tpl2 <-
    load_contrast(
        tpl2_file,
        "TPL2i"
    )

mek <-
    load_contrast(
        mek_file,
        "MEKi"
    )

acute_lps <-
    load_contrast(
        acute_lps_file,
        "AcuteLPS"
    )

# ============================================================
# 7. MERGE BOTH STUDIES BY HUMAN GENE SYMBOL
# ============================================================

memory_small <-
    memory[
        ,
        c(
            "SYMBOL",
            "Memory_log2FC",
            "Memory_P",
            "Memory_FDR",
            "Memory_FDR05",
            "Memory_robust",
            "Memory_direction",
            "Robust_memory_class"
        )
    ]

integrated <-
    merge(
        memory_small,
        acute_lps,
        by = "SYMBOL",
        all = FALSE
    )

integrated <-
    merge(
        integrated,
        tpl2,
        by = "SYMBOL",
        all = FALSE
    )

integrated <-
    merge(
        integrated,
        mek,
        by = "SYMBOL",
        all = FALSE
    )

cat(
    "Genes shared between GSE85246 and GSE84161:",
    nrow(integrated),
    "\n\n"
)

# ============================================================
# 8. DIRECTIONAL REVERSAL FLAGS
# ============================================================

integrated$TPL2_opposite <-
    sign(integrated$Memory_log2FC) !=
    sign(integrated$TPL2i_log2FC)

integrated$TPL2_same <-
    sign(integrated$Memory_log2FC) ==
    sign(integrated$TPL2i_log2FC)

integrated$TPL2_significant <-
    integrated$TPL2i_FDR < 0.05

integrated$TPL2_significant_reversal <-
    integrated$Memory_robust &
    integrated$TPL2_significant &
    integrated$TPL2_opposite

integrated$MEK_opposite <-
    sign(integrated$Memory_log2FC) !=
    sign(integrated$MEKi_log2FC)

integrated$MEK_same <-
    sign(integrated$Memory_log2FC) ==
    sign(integrated$MEKi_log2FC)

integrated$MEK_significant <-
    integrated$MEKi_FDR < 0.05

integrated$MEK_significant_reversal <-
    integrated$Memory_robust &
    integrated$MEK_significant &
    integrated$MEK_opposite

write.csv(
    integrated,
    file.path(
        outdir,
        "09c_GSE85246_GSE84161_integrated_gene_table.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 9. DESCRIPTIVE REVERSAL SUMMARY
# ============================================================

make_summary <- function(
    effect,
    fdr,
    name
) {

    robust <-
        integrated$Memory_robust

    up <-
        robust &
        integrated$Memory_log2FC > 0

    down <-
        robust &
        integrated$Memory_log2FC < 0

    data.frame(
        Perturbation = name,

        Robust_memory_genes =
            sum(robust),

        Memory_UP =
            sum(up),

        Memory_DOWN =
            sum(down),

        Opposite_direction_all =
            sum(
                robust &
                sign(
                    integrated$Memory_log2FC
                ) != sign(effect)
            ),

        Same_direction_all =
            sum(
                robust &
                sign(
                    integrated$Memory_log2FC
                ) == sign(effect)
            ),

        Significant_perturbation =
            sum(
                robust &
                fdr < 0.05
            ),

        Significant_reversal =
            sum(
                robust &
                fdr < 0.05 &
                sign(
                    integrated$Memory_log2FC
                ) != sign(effect)
            ),

        Significant_same_direction =
            sum(
                robust &
                fdr < 0.05 &
                sign(
                    integrated$Memory_log2FC
                ) == sign(effect)
            )
    )
}

reversal_summary <-
    rbind(
        make_summary(
            integrated$TPL2i_log2FC,
            integrated$TPL2i_FDR,
            "TPL2_inhibition"
        ),

        make_summary(
            integrated$MEKi_log2FC,
            integrated$MEKi_FDR,
            "MEK_inhibition"
        )
    )

write.csv(
    reversal_summary,
    file.path(
        outdir,
        "09c_descriptive_reversal_summary.csv"
    ),
    row.names = FALSE
)

cat("========================================\n")
cat("DESCRIPTIVE REVERSAL SUMMARY\n")
cat("========================================\n\n")

print(
    reversal_summary,
    row.names = FALSE
)

# ============================================================
# 10. EFFECT-SIZE CORRELATIONS
#
# Negative correlation:
#   compatible with reversal
#
# Positive correlation:
#   compatible with concordance
#
# Near zero:
#   little global relationship
#
# Cross-platform absolute FC magnitudes are NOT directly
# equated; rank correlation is therefore used.
# ============================================================

run_cor <- function(
    x,
    y,
    subset,
    comparison
) {

    ok <-
        subset &
        is.finite(x) &
        is.finite(y)

    test <-
        suppressWarnings(
            cor.test(
                x[ok],
                y[ok],
                method = "spearman",
                exact = FALSE
            )
        )

    data.frame(
        Comparison = comparison,
        N_genes = sum(ok),
        Spearman_rho =
            unname(test$estimate),
        P_value =
            test$p.value
    )
}

cor_results <-
    rbind(
        run_cor(
            integrated$Memory_log2FC,
            integrated$TPL2i_log2FC,
            rep(TRUE, nrow(integrated)),
            "All genes: Memory vs TPL2 inhibition"
        ),

        run_cor(
            integrated$Memory_log2FC,
            integrated$TPL2i_log2FC,
            integrated$Memory_robust,
            "Robust memory genes: Memory vs TPL2 inhibition"
        ),

        run_cor(
            integrated$Memory_log2FC,
            integrated$MEKi_log2FC,
            rep(TRUE, nrow(integrated)),
            "All genes: Memory vs MEK inhibition"
        ),

        run_cor(
            integrated$Memory_log2FC,
            integrated$MEKi_log2FC,
            integrated$Memory_robust,
            "Robust memory genes: Memory vs MEK inhibition"
        ),

        run_cor(
            integrated$TPL2i_log2FC,
            integrated$MEKi_log2FC,
            rep(TRUE, nrow(integrated)),
            "All genes: TPL2 inhibition vs MEK inhibition"
        ),

        run_cor(
            integrated$TPL2i_log2FC,
            integrated$MEKi_log2FC,
            integrated$Memory_robust,
            "Robust memory genes: TPL2 inhibition vs MEK inhibition"
        )
    )

cor_results$FDR <-
    p.adjust(
        cor_results$P_value,
        method = "BH"
    )

write.csv(
    cor_results,
    file.path(
        outdir,
        "09c_effect_correlation_results.csv"
    ),
    row.names = FALSE
)

cat("\n========================================\n")
cat("EFFECT CORRELATIONS\n")
cat("========================================\n\n")

print(
    cor_results,
    row.names = FALSE
)

# ============================================================
# 11. LOAD GSE84161 EXPRESSION FOR CAMERA
# ============================================================

expr_df <-
    read.csv(
        expr_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

if (!"SYMBOL" %in% names(expr_df)) {
    stop(
        "SYMBOL column absent from GSE84161 gene expression file."
    )
}

expr <-
    as.matrix(
        expr_df[
            ,
            setdiff(
                names(expr_df),
                "SYMBOL"
            ),
            drop = FALSE
        ]
    )

rownames(expr) <-
    expr_df$SYMBOL

storage.mode(expr) <- "numeric"

manifest <-
    read.csv(
        manifest_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

idx <-
    match(
        manifest$GSM,
        colnames(expr)
    )

if (anyNA(idx)) {
    stop(
        "Could not match GSE84161 expression columns to manifest."
    )
}

expr <-
    expr[
        ,
        idx,
        drop = FALSE
    ]

stopifnot(
    identical(
        colnames(expr),
        manifest$GSM
    )
)

# ============================================================
# 12. REBUILD PAIRED GSE84161 DESIGN
# ============================================================

condition <-
    factor(
        manifest$condition,
        levels = c(
            "Untreated",
            "TPL2i",
            "MEKi",
            "LPS",
            "TPL2i_LPS",
            "MEKi_LPS"
        )
    )

donor <-
    factor(
        manifest$donor
    )

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
    manifest$GSM

if (
    qr(design)$rank !=
    ncol(design)
) {
    stop(
        "CAMERA design matrix is not full rank."
    )
}

contrasts <-
    limma::makeContrasts(

        TPL2i_LPS_vs_LPS =
            TPL2i_LPS - LPS,

        MEKi_LPS_vs_LPS =
            MEKi_LPS - LPS,

        levels = design
    )

# ============================================================
# 13. BUILD PRE-DEFINED MEMORY GENE SETS
# ============================================================

gene_sets <- list(

    Robust_Memory_UP =
        memory$SYMBOL[
            memory$Memory_FDR < 0.05 &
            memory$Memory_log2FC >= 1
        ],

    Robust_Memory_DOWN =
        memory$SYMBOL[
            memory$Memory_FDR < 0.05 &
            memory$Memory_log2FC <= -1
        ],

    FDR05_Memory_UP =
        memory$SYMBOL[
            memory$Memory_FDR < 0.05 &
            memory$Memory_log2FC > 0
        ],

    FDR05_Memory_DOWN =
        memory$SYMBOL[
            memory$Memory_FDR < 0.05 &
            memory$Memory_log2FC < 0
        ]
)

gene_sets <-
    lapply(
        gene_sets,
        function(x) {
            intersect(
                unique(x),
                rownames(expr)
            )
        }
    )

cat("\n========================================\n")
cat("MEMORY GENE SET SIZES IN GSE84161\n")
cat("========================================\n\n")

print(
    sapply(
        gene_sets,
        length
    )
)

indices <-
    lapply(
        gene_sets,
        function(gs) {
            match(
                gs,
                rownames(expr)
            )
        }
    )

indices <-
    lapply(
        indices,
        function(x) {
            x[!is.na(x)]
        }
    )

# ============================================================
# 14. CAMERA COMPETITIVE GENE-SET TEST
#
# For TPL2i_LPS_vs_LPS:
#
# Supportive of reversal if:
#
# Robust_Memory_UP   -> Direction DOWN
# Robust_Memory_DOWN -> Direction UP
#
# Same logic for MEK.
# ============================================================

run_camera <- function(
    contrast_name
) {

    result <-
        limma::camera(
            expr,
            index = indices,
            design = design,
            contrast =
                contrasts[
                    ,
                    contrast_name
                ],
            sort = FALSE
        )

    data.frame(
        Contrast =
            contrast_name,
        GeneSet =
            rownames(result),
        result,
        row.names = NULL,
        check.names = FALSE
    )
}

camera_tpl2 <-
    run_camera(
        "TPL2i_LPS_vs_LPS"
    )

camera_mek <-
    run_camera(
        "MEKi_LPS_vs_LPS"
    )

camera_results <-
    rbind(
        camera_tpl2,
        camera_mek
    )

write.csv(
    camera_results,
    file.path(
        outdir,
        "09c_CAMERA_memory_signature_results.csv"
    ),
    row.names = FALSE
)

cat("\n========================================\n")
cat("CAMERA MEMORY-SIGNATURE TEST\n")
cat("========================================\n\n")

print(
    camera_results,
    row.names = FALSE
)

# ============================================================
# 15. MAP3K8-SPECIFIC CROSS-STUDY ROW
# ============================================================

map3k8 <-
    integrated[
        integrated$SYMBOL ==
            "MAP3K8",
        ,
        drop = FALSE
    ]

cat("\n========================================\n")
cat("MAP3K8 CROSS-STUDY RESULT\n")
cat("========================================\n\n")

print(
    map3k8[
        ,
        c(
            "SYMBOL",
            "Memory_log2FC",
            "Memory_FDR",
            "AcuteLPS_log2FC",
            "AcuteLPS_FDR",
            "TPL2i_log2FC",
            "TPL2i_FDR",
            "MEKi_log2FC",
            "MEKi_FDR"
        )
    ],
    row.names = FALSE
)

write.csv(
    map3k8,
    file.path(
        outdir,
        "09c_MAP3K8_cross_study_result.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 16. SCATTERPLOT — MEMORY vs TPL2 INHIBITION
# ============================================================

png(
    file.path(
        figdir,
        "09C1_Day6_memory_vs_TPL2_inhibition.png"
    ),
    width = 1400,
    height = 1200,
    res = 160
)

plot(
    integrated$Memory_log2FC,
    integrated$TPL2i_log2FC,
    pch = 16,
    cex = 0.45,
    xlab =
        "GSE85246 Day-6 LPS vs RPMI log2FC",
    ylab =
        "GSE84161 TPL2i+LPS vs LPS log2FC",
    main =
        "Persistent memory effect vs TPL2 inhibition"
)

abline(
    h = 0,
    lty = 2
)

abline(
    v = 0,
    lty = 2
)

rob <-
    integrated$Memory_robust

points(
    integrated$Memory_log2FC[rob],
    integrated$TPL2i_log2FC[rob],
    pch = 1,
    cex = 0.8
)

if (any(
    integrated$SYMBOL ==
        "MAP3K8"
)) {

    m <-
        integrated[
            integrated$SYMBOL ==
                "MAP3K8",
            ,
            drop = FALSE
        ]

    text(
        m$Memory_log2FC,
        m$TPL2i_log2FC,
        labels = "MAP3K8",
        pos = 3
    )
}

dev.off()

# ============================================================
# 17. SCATTERPLOT — MEMORY vs MEK INHIBITION
# ============================================================

png(
    file.path(
        figdir,
        "09C2_Day6_memory_vs_MEK_inhibition.png"
    ),
    width = 1400,
    height = 1200,
    res = 160
)

plot(
    integrated$Memory_log2FC,
    integrated$MEKi_log2FC,
    pch = 16,
    cex = 0.45,
    xlab =
        "GSE85246 Day-6 LPS vs RPMI log2FC",
    ylab =
        "GSE84161 MEKi+LPS vs LPS log2FC",
    main =
        "Persistent memory effect vs MEK inhibition"
)

abline(
    h = 0,
    lty = 2
)

abline(
    v = 0,
    lty = 2
)

points(
    integrated$Memory_log2FC[rob],
    integrated$MEKi_log2FC[rob],
    pch = 1,
    cex = 0.8
)

dev.off()

# ============================================================
# 18. FINAL MACHINE-READABLE SUMMARY
# ============================================================

sink(
    file.path(
        outdir,
        "09c_validation_summary.txt"
    )
)

cat(
    "GSE85246 vs GSE84161 MAP3K8/TPL2 validation\n\n"
)

cat(
    "Shared genes:",
    nrow(integrated),
    "\n"
)

cat(
    "Robust Day-6 memory genes:",
    sum(
        integrated$Memory_robust
    ),
    "\n\n"
)

cat(
    "Descriptive reversal summary:\n"
)

print(
    reversal_summary,
    row.names = FALSE
)

cat(
    "\nEffect correlations:\n"
)

print(
    cor_results,
    row.names = FALSE
)

cat(
    "\nCAMERA results:\n"
)

print(
    camera_results,
    row.names = FALSE
)

cat(
    "\nMAP3K8:\n"
)

print(
    map3k8[
        ,
        c(
            "SYMBOL",
            "Memory_log2FC",
            "Memory_FDR",
            "AcuteLPS_log2FC",
            "AcuteLPS_FDR",
            "TPL2i_log2FC",
            "TPL2i_FDR",
            "MEKi_log2FC",
            "MEKi_FDR"
        )
    ],
    row.names = FALSE
)

sink()

cat("\n========================================\n")
cat("09c COMPLETE\n")
cat("========================================\n\n")

cat(
    "IMPORTANT:\n",
    "Do not interpret this as proof of memory maintenance.\n",
    "GSE84161 tests acute TPL2 dependence.\n",
    "We are testing whether the persistent GSE85246\n",
    "memory signature is functionally connected to\n",
    "TPL2/MAP3K8-regulated transcription.\n"
)
