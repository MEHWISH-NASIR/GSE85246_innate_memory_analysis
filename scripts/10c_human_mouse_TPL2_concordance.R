# ============================================================
# 10c_human_mouse_TPL2_concordance.R
#
# CROSS-SPECIES FUNCTIONAL VALIDATION OF MAP3K8/TPL2
#
# HUMAN:
#   GSE84161 primary monocytes
#   pharmacological TPL2 inhibition
#
# MOUSE:
#   GSE116220 BMDMs
#   kinase-inactive Map3k8[D270A]
#
# PRIMARY COMPARISON:
#
# Human:
# (TPL2i+LPS - TPL2i) - (LPS - Untreated)
#
# Mouse:
# (D270A_LPS - D270A_unstim)
#      -
# (WT_LPS - WT_unstim)
#
# Only high-confidence one-to-one Ensembl orthologues are used.
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
    "biomaRt"
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

human_expr_file <-
    "results/09_MAP3K8_validation/09b_GSE84161_gene_level_log2_expression.csv"

human_manifest_file <-
    "results/09_MAP3K8_validation/09b_GSE84161_analysis_manifest.csv"

mouse_files <- c(
    "0.5hr" =
        "results/10_MAP3K8_genetic_validation/10b_Interaction_0.5hr.csv",

    "1hr" =
        "results/10_MAP3K8_genetic_validation/10b_Interaction_1hr.csv",

    "2hr" =
        "results/10_MAP3K8_genetic_validation/10b_Interaction_2hr.csv"
)

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

required_files <- c(
    human_expr_file,
    human_manifest_file,
    mouse_files
)

if (!all(file.exists(required_files))) {

    cat("Missing files:\n")

    print(
        required_files[
            !file.exists(required_files)
        ]
    )

    stop("Required file missing.")
}

cat("\n========================================\n")
cat("10c — HUMAN/MOUSE TPL2 CONCORDANCE\n")
cat("========================================\n\n")

# ============================================================
# 3. HUMAN GSE84161 EXPRESSION
# ============================================================

human_df <- read.csv(
    human_expr_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (!"SYMBOL" %in% names(human_df)) {
    stop("Human expression file lacks SYMBOL.")
}

human_expr <- as.matrix(
    human_df[
        ,
        setdiff(
            names(human_df),
            "SYMBOL"
        ),
        drop = FALSE
    ]
)

storage.mode(human_expr) <- "numeric"

rownames(human_expr) <-
    human_df$SYMBOL

human_manifest <- read.csv(
    human_manifest_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_manifest_cols <- c(
    "GSM",
    "donor",
    "condition"
)

if (!all(
    required_manifest_cols %in%
        names(human_manifest)
)) {
    stop("Unexpected human manifest columns.")
}

idx <- match(
    human_manifest$GSM,
    colnames(human_expr)
)

if (anyNA(idx)) {
    stop("Human expression/manifest mismatch.")
}

human_expr <-
    human_expr[
        ,
        idx,
        drop = FALSE
    ]

stopifnot(
    identical(
        colnames(human_expr),
        human_manifest$GSM
    )
)

# ============================================================
# 4. HUMAN PAIRED DESIGN
# ============================================================

condition <- factor(
    human_manifest$condition,
    levels = c(
        "Untreated",
        "TPL2i",
        "MEKi",
        "LPS",
        "TPL2i_LPS",
        "MEKi_LPS"
    )
)

donor <- factor(
    human_manifest$donor
)

design_h <- model.matrix(
    ~ 0 + condition + donor
)

colnames(design_h) <- sub(
    "^condition",
    "",
    colnames(design_h)
)

rownames(design_h) <-
    human_manifest$GSM

if (
    qr(design_h)$rank !=
    ncol(design_h)
) {
    stop("Human design is not full rank.")
}

cat(
    "Human design rank:",
    qr(design_h)$rank,
    "/",
    ncol(design_h),
    "\n"
)

# ============================================================
# 5. HUMAN TPL2 x LPS INTERACTION
#
# (TPL2i_LPS - TPL2i)
#           -
# (LPS - Untreated)
#
# Negative effect:
# TPL2 inhibition SUPPRESSES the LPS response.
# Gene is normally positively supported by TPL2 activity.
#
# Positive effect:
# TPL2 inhibition ENHANCES the LPS response.
# Gene is normally restrained by TPL2 activity.
# ============================================================

contrast_h <- limma::makeContrasts(

    TPL2_LPS_interaction =
        (
            TPL2i_LPS -
            TPL2i
        ) -
        (
            LPS -
            Untreated
        ),

    levels = design_h
)

fit_h <- limma::lmFit(
    human_expr,
    design_h
)

fit_h <- limma::contrasts.fit(
    fit_h,
    contrast_h
)

fit_h <- limma::eBayes(
    fit_h,
    trend = TRUE,
    robust = TRUE
)

human <- limma::topTable(
    fit_h,
    coef = "TPL2_LPS_interaction",
    number = Inf,
    sort.by = "P",
    adjust.method = "BH"
)

human$SYMBOL <-
    rownames(human)

names(human)[
    names(human) == "logFC"
] <- "Human_effect"

names(human)[
    names(human) == "adj.P.Val"
] <- "Human_FDR"

names(human)[
    names(human) == "t"
] <- "Human_t"

names(human)[
    names(human) == "P.Value"
] <- "Human_P"

human$Human_direction <- ifelse(
    human$Human_effect > 0,
    "UP",
    ifelse(
        human$Human_effect < 0,
        "DOWN",
        "UNCHANGED"
    )
)

human <- human[
    ,
    c(
        "SYMBOL",
        "Human_effect",
        "Human_t",
        "Human_P",
        "Human_FDR",
        "Human_direction"
    )
]

write.csv(
    human,
    file.path(
        outdir,
        "10c_GSE84161_TPL2_LPS_interaction.csv"
    ),
    row.names = FALSE
)

cat(
    "\nHuman interaction FDR<0.05:",
    sum(
        human$Human_FDR < 0.05,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Human interaction DOWN:",
    sum(
        human$Human_FDR < 0.05 &
        human$Human_effect < 0,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Human interaction UP:",
    sum(
        human$Human_FDR < 0.05 &
        human$Human_effect > 0,
        na.rm = TRUE
    ),
    "\n\n"
)

# ============================================================
# 6. LOAD MOUSE INTERACTION RESULTS
# ============================================================

load_mouse <- function(
    file,
    time_label
) {

    x <- read.csv(
        file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    required <- c(
        "SYMBOL",
        "Effect",
        "t",
        "P.Value",
        "FDR"
    )

    if (!all(required %in% names(x))) {
        stop(
            paste(
                "Unexpected mouse columns:",
                file
            )
        )
    }

    out <- x[
        ,
        required,
        drop = FALSE
    ]

    names(out) <- c(
        "Mouse_SYMBOL",
        "Mouse_effect",
        "Mouse_t",
        "Mouse_P",
        "Mouse_FDR"
    )

    out$Time <- time_label

    out
}

mouse <- lapply(
    names(mouse_files),
    function(tm) {

        load_mouse(
            mouse_files[[tm]],
            tm
        )
    }
)

names(mouse) <-
    names(mouse_files)

mouse_symbols <- unique(
    unlist(
        lapply(
            mouse,
            function(x) {
                x$Mouse_SYMBOL
            }
        )
    )
)

cat(
    "Mouse symbols requiring orthology mapping:",
    length(mouse_symbols),
    "\n"
)

# ============================================================
# 7. MGI / ALLIANCE ONE-TO-ONE HUMAN-MOUSE ORTHOLOGY
#
# BioMart archive access returned HTTP 405.
#
# We therefore use the official MGI static report:
# HOM_ProteinCoding.rpt
#
# The report contains protein-coding mouse genes with
# one-to-one human orthology from the Alliance of Genome
# Resources.
# ============================================================

cat("\n========================================\n")
cat("LOADING MGI ONE-TO-ONE ORTHOLOGUES\n")
cat("========================================\n\n")

orth_file <-
    "data/external/GSE116220/HOM_ProteinCoding.rpt"

if (!file.exists(orth_file)) {
    stop(
        paste(
            "Missing MGI orthology file:",
            orth_file
        )
    )
}

mgi_orth <- read.delim(
    orth_file,
    header = FALSE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    fill = TRUE
)

if (ncol(mgi_orth) != 6) {
    stop(
        paste(
            "Expected 6 columns in HOM_ProteinCoding.rpt;",
            "found",
            ncol(mgi_orth)
        )
    )
}

names(mgi_orth) <- c(
    "MGI_ID",
    "Mouse_SYMBOL",
    "Mouse_Entrez",
    "HGNC_ID",
    "Human_SYMBOL",
    "Human_Entrez"
)

mgi_orth <- mgi_orth[
    mgi_orth$Mouse_SYMBOL != "" &
    mgi_orth$Human_SYMBOL != "",
    ,
    drop = FALSE
]

# ------------------------------------------------------------
# Resolve symbol-level duplicates conservatively
#
# HOM_ProteinCoding.rpt represents one-to-one orthology at
# the gene-record level. Because our expression matrices are
# keyed by gene SYMBOL, we additionally require an unambiguous
# symbol-to-symbol relationship.
#
# We do NOT arbitrarily choose the first ambiguous mapping.
# ------------------------------------------------------------

raw_orth_n <- nrow(mgi_orth)

# First remove completely duplicated records/pairs
mgi_orth <- unique(
    mgi_orth
)

after_exact_dedup_n <- nrow(mgi_orth)

# Number of distinct human symbols associated with each
# mouse symbol
mouse_human_n <- ave(
    mgi_orth$Human_SYMBOL,
    mgi_orth$Mouse_SYMBOL,
    FUN = function(z) {
        length(unique(z))
    }
)

# Number of distinct mouse symbols associated with each
# human symbol
human_mouse_n <- ave(
    mgi_orth$Mouse_SYMBOL,
    mgi_orth$Human_SYMBOL,
    FUN = function(z) {
        length(unique(z))
    }
)

ambiguous_mouse_symbols <- unique(
    mgi_orth$Mouse_SYMBOL[
        mouse_human_n > 1
    ]
)

ambiguous_human_symbols <- unique(
    mgi_orth$Human_SYMBOL[
        human_mouse_n > 1
    ]
)

cat(
    "Raw MGI rows:",
    raw_orth_n,
    "\n"
)

cat(
    "Exact duplicate rows removed:",
    raw_orth_n - after_exact_dedup_n,
    "\n"
)

cat(
    "Mouse symbols mapping to >1 human symbol:",
    length(ambiguous_mouse_symbols),
    "\n"
)

cat(
    "Human symbols mapping to >1 mouse symbol:",
    length(ambiguous_human_symbols),
    "\n"
)

if (length(ambiguous_mouse_symbols) > 0) {

    cat(
        "\nFirst ambiguous mouse symbols:\n"
    )

    print(
        head(
            ambiguous_mouse_symbols,
            20
        )
    )
}

if (length(ambiguous_human_symbols) > 0) {

    cat(
        "\nFirst ambiguous human symbols:\n"
    )

    print(
        head(
            ambiguous_human_symbols,
            20
        )
    )
}

# Keep only unambiguous symbol-level 1:1 mappings
mgi_orth <- mgi_orth[
    mouse_human_n == 1 &
    human_mouse_n == 1,
    ,
    drop = FALSE
]

# The same mouse-human symbol pair may still have duplicated
# identifier rows. Collapse those to exactly one symbol pair.
pair_key <- paste(
    mgi_orth$Mouse_SYMBOL,
    mgi_orth$Human_SYMBOL,
    sep = "|||"
)

mgi_orth <- mgi_orth[
    !duplicated(pair_key),
    ,
    drop = FALSE
]

# Final safety checks
if (anyDuplicated(mgi_orth$Mouse_SYMBOL)) {
    stop(
        "Mouse symbols remain duplicated after ambiguity filtering."
    )
}

if (anyDuplicated(mgi_orth$Human_SYMBOL)) {
    stop(
        "Human symbols remain duplicated after ambiguity filtering."
    )
}

cat(
    "Final unambiguous symbol-level one-to-one pairs:",
    nrow(mgi_orth),
    "\n\n"
)

cat(
    "MGI one-to-one orthologue pairs:",
    nrow(mgi_orth),
    "\n"
)

map3k8_check <- mgi_orth[
    mgi_orth$Mouse_SYMBOL == "Map3k8",
    ,
    drop = FALSE
]

cat("\nMAP3K8 orthology check:\n\n")

print(
    map3k8_check,
    row.names = FALSE
)

if (
    nrow(map3k8_check) != 1 ||
    map3k8_check$Human_SYMBOL != "MAP3K8"
) {
    stop(
        "MAP3K8 one-to-one orthology check failed."
    )
}

# Keep names expected by downstream integration code.
orth1 <- data.frame(
    Mouse_SYMBOL =
        mgi_orth$Mouse_SYMBOL,

    Human_SYMBOL =
        mgi_orth$Human_SYMBOL,

    Orthology_type =
        "MGI_Alliance_one2one",

    Orthology_confidence =
        1,

    Mouse_percent_identity =
        NA_real_,

    stringsAsFactors = FALSE
)

write.csv(
    mgi_orth,
    file.path(
        outdir,
        "10c_MGI_Alliance_one2one_mouse_human_orthologues.csv"
    ),
    row.names = FALSE
)

cat(
    "\nMouse expression genes:",
    length(mouse_symbols),
    "\n"
)

cat(
    "Available MGI one-to-one mappings:",
    nrow(orth1),
    "\n\n"
)

# ============================================================
# 9. BUILD CROSS-SPECIES TABLES
# ============================================================

integrate_time <- function(
    mouse_result,
    time_label
) {

    temp <- merge(
        orth1,
        mouse_result,
        by = "Mouse_SYMBOL",
        all = FALSE
    )

    temp <- merge(
        temp,
        human,
        by.x = "Human_SYMBOL",
        by.y = "SYMBOL",
        all = FALSE
    )

    temp$Same_direction <-
        sign(
            temp$Human_effect
        ) ==
        sign(
            temp$Mouse_effect
        )

    temp$Human_significant <-
        temp$Human_FDR < 0.05

    temp$Mouse_significant <-
        temp$Mouse_FDR < 0.05

    temp$Both_significant <-
        temp$Human_significant &
        temp$Mouse_significant

    temp$Time <- time_label

    temp
}

integrated <- lapply(
    names(mouse),
    function(tm) {

        integrate_time(
            mouse[[tm]],
            tm
        )
    }
)

names(integrated) <-
    names(mouse)

for (tm in names(integrated)) {

    write.csv(
        integrated[[tm]],
        file.path(
            outdir,
            paste0(
                "10c_human_mouse_integrated_",
                tm,
                ".csv"
            )
        ),
        row.names = FALSE
    )
}

# ============================================================
# 10. STATISTICAL HELPERS
# ============================================================

safe_cor <- function(
    x,
    y
) {

    ok <-
        is.finite(x) &
        is.finite(y)

    if (sum(ok) < 3) {

        return(
            c(
                N = sum(ok),
                rho = NA,
                P = NA
            )
        )
    }

    z <- suppressWarnings(
        cor.test(
            x[ok],
            y[ok],
            method = "spearman",
            exact = FALSE
        )
    )

    c(
        N = sum(ok),
        rho = unname(
            z$estimate
        ),
        P = z$p.value
    )
}

safe_concordance <- function(
    same
) {

    same <- same[
        !is.na(same)
    ]

    n <- length(same)

    if (n == 0) {

        return(
            c(
                N = 0,
                Same = 0,
                Fraction = NA,
                P = NA
            )
        )
    }

    k <- sum(same)

    bt <- binom.test(
        k,
        n,
        p = 0.5,
        alternative = "greater"
    )

    c(
        N = n,
        Same = k,
        Fraction = k / n,
        P = bt$p.value
    )
}

# ============================================================
# 11. CORRELATION + SIGN CONCORDANCE
# ============================================================

summary_rows <- list()

for (tm in names(integrated)) {

    x <- integrated[[tm]]

    all_cor <- safe_cor(
        x$Human_effect,
        x$Mouse_effect
    )

    human_sig <- x$Human_FDR < 0.05

    hs_cor <- safe_cor(
        x$Human_effect[human_sig],
        x$Mouse_effect[human_sig]
    )

    hs_sign <- safe_concordance(
        x$Same_direction[
            human_sig
        ]
    )

    both_sig <-
        x$Human_FDR < 0.05 &
        x$Mouse_FDR < 0.05

    both_sign <- safe_concordance(
        x$Same_direction[
            both_sig
        ]
    )

    summary_rows[[tm]] <- data.frame(
        Mouse_time = tm,

        Shared_one2one_genes =
            as.integer(
                all_cor["N"]
            ),

        All_genes_Spearman_rho =
            as.numeric(
                all_cor["rho"]
            ),

        All_genes_Spearman_P =
            as.numeric(
                all_cor["P"]
            ),

        Human_FDR05_genes =
            sum(
                human_sig,
                na.rm = TRUE
            ),

        Human_FDR05_Spearman_rho =
            as.numeric(
                hs_cor["rho"]
            ),

        Human_FDR05_Spearman_P =
            as.numeric(
                hs_cor["P"]
            ),

        Human_FDR05_same_direction =
            as.integer(
                hs_sign["Same"]
            ),

        Human_FDR05_sign_fraction =
            as.numeric(
                hs_sign["Fraction"]
            ),

        Human_FDR05_sign_P =
            as.numeric(
                hs_sign["P"]
            ),

        Both_FDR05_genes =
            sum(
                both_sig,
                na.rm = TRUE
            ),

        Both_FDR05_same_direction =
            as.integer(
                both_sign["Same"]
            ),

        Both_FDR05_sign_fraction =
            as.numeric(
                both_sign["Fraction"]
            ),

        Both_FDR05_sign_P =
            as.numeric(
                both_sign["P"]
            ),

        stringsAsFactors = FALSE
    )
}

summary_table <-
    do.call(
        rbind,
        summary_rows
    )

summary_table$All_genes_Spearman_FDR <-
    p.adjust(
        summary_table$All_genes_Spearman_P,
        method = "BH"
    )

summary_table$Human_FDR05_Spearman_FDR <-
    p.adjust(
        summary_table$Human_FDR05_Spearman_P,
        method = "BH"
    )

summary_table$Human_FDR05_sign_FDR <-
    p.adjust(
        summary_table$Human_FDR05_sign_P,
        method = "BH"
    )

summary_table$Both_FDR05_sign_FDR <-
    p.adjust(
        summary_table$Both_FDR05_sign_P,
        method = "BH"
    )

write.csv(
    summary_table,
    file.path(
        outdir,
        "10c_human_mouse_concordance_summary.csv"
    ),
    row.names = FALSE
)

cat("========================================\n")
cat("CROSS-SPECIES CONCORDANCE\n")
cat("========================================\n\n")

print(
    summary_table,
    row.names = FALSE
)

# ============================================================
# 12. PRE-RANKED GENE-SET TEST
#
# Human genes significantly affected by TPL2 inhibition are
# split by direction.
#
# Human DOWN:
# normally promoted by TPL2.
# Expected mouse direction = DOWN in D270A.
#
# Human UP:
# normally suppressed by TPL2.
# Expected mouse direction = UP in D270A.
# ============================================================

camera_results <- list()

for (tm in names(integrated)) {

    x <- integrated[[tm]]

    statistic <-
        x$Mouse_t

    names(statistic) <-
        x$Human_SYMBOL

    index <- list(

        Human_TPL2_promoted = which(
            x$Human_FDR < 0.05 &
            x$Human_effect < 0
        ),

        Human_TPL2_suppressed = which(
            x$Human_FDR < 0.05 &
            x$Human_effect > 0
        )
    )

    index <- index[
        vapply(
            index,
            length,
            integer(1)
        ) >= 2
    ]

    if (length(index) == 0) {
        next
    }

    cam <- limma::cameraPR(
        statistic = statistic,
        index = index,
        sort = FALSE
    )

    cam <- data.frame(
        Mouse_time = tm,
        GeneSet = rownames(cam),
        cam,
        row.names = NULL,
        check.names = FALSE
    )

    camera_results[[tm]] <-
        cam
}

camera_table <-
    do.call(
        rbind,
        camera_results
    )

if (!is.null(camera_table)) {

    camera_table$Global_FDR <-
        p.adjust(
            camera_table$PValue,
            method = "BH"
        )

    write.csv(
        camera_table,
        file.path(
            outdir,
            "10c_cameraPR_human_TPL2_sets_in_mouse.csv"
        ),
        row.names = FALSE
    )

    cat("\n========================================\n")
    cat("GENE-SET CONCORDANCE\n")
    cat("========================================\n\n")

    print(
        camera_table,
        row.names = FALSE
    )
}

# ============================================================
# 13. SCATTERPLOTS
# ============================================================

for (tm in names(integrated)) {

    x <- integrated[[tm]]

    png(
        file.path(
            figdir,
            paste0(
                "10C_human_TPL2_vs_mouse_D270A_",
                tm,
                ".png"
            )
        ),
        width = 1400,
        height = 1200,
        res = 160
    )

    plot(
        x$Human_effect,
        x$Mouse_effect,
        pch = 16,
        cex = 0.4,

        xlab =
            "Human TPL2i x LPS interaction effect",

        ylab =
            paste0(
                "Mouse D270A x LPS interaction effect (",
                tm,
                ")"
            ),

        main =
            paste0(
                "TPL2/MAP3K8 perturbation concordance: ",
                tm
            )
    )

    abline(
        h = 0,
        lty = 2
    )

    abline(
        v = 0,
        lty = 2
    )

    sig <- x$Human_FDR < 0.05

    points(
        x$Human_effect[sig],
        x$Mouse_effect[sig],
        pch = 1,
        cex = 0.75
    )

    dev.off()
}

# ============================================================
# 14. SELECT TOP CONCORDANT GENES
# ============================================================

top_tables <- list()

for (tm in names(integrated)) {

    x <- integrated[[tm]]

    hits <- x[
        x$Human_FDR < 0.05 &
        x$Mouse_FDR < 0.05 &
        x$Same_direction,
        ,
        drop = FALSE
    ]

    hits$combined_rank_score <-
        -log10(
            pmax(
                hits$Human_FDR,
                1e-300
            )
        ) +
        -log10(
            pmax(
                hits$Mouse_FDR,
                1e-300
            )
        )

    hits <- hits[
        order(
            -hits$combined_rank_score
        ),
        ,
        drop = FALSE
    ]

    top_tables[[tm]] <-
        hits

    write.csv(
        hits,
        file.path(
            outdir,
            paste0(
                "10c_concordant_both_FDR05_",
                tm,
                ".csv"
            )
        ),
        row.names = FALSE
    )
}

# ============================================================
# 15. FINAL SUMMARY
# ============================================================

cat("\n========================================\n")
cat("10c COMPLETE\n")
cat("========================================\n\n")

cat(
    "Interpretation:\n\n",
    "Positive human-mouse correlation and significant\n",
    "same-direction enrichment support conservation of a\n",
    "TPL2/MAP3K8-dependent LPS transcriptional program.\n\n",

    "Human TPL2-promoted genes should trend DOWN after\n",
    "mouse Map3k8 catalytic inactivation.\n\n",

    "Human TPL2-suppressed genes should trend UP after\n",
    "mouse Map3k8 catalytic inactivation.\n\n",

    "This validates TPL2/MAP3K8 pathway function;\n",
    "it does NOT by itself prove Day-6 memory maintenance.\n"
)
