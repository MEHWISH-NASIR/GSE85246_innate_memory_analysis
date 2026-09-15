# ============================================================
# 11b_TCGA_MAP3K8_overexpression_frequency.R
#
# MAP3K8 PAN-CANCER OVEREXPRESSION
#
# Expression:
# UCSC Xena TCGA-TARGET-GTEx Toil
# log2(TPM + 0.001)
#
# Normal-tissue matching:
# UCSCXenaShiny curated tcga_gtex mapping
#
# QUESTIONS:
# 1. Which TCGA cancers overexpress MAP3K8?
# 2. How large is the difference?
# 3. What fraction of tumors are MAP3K8-high?
#
# Primary frequency definition:
# tumor expression > 95th percentile of corresponding normals
# ============================================================

rm(list = ls())

# ============================================================
# 1. PACKAGE
# ============================================================

if (!requireNamespace("UCSCXenaShiny", quietly = TRUE)) {
    stop("UCSCXenaShiny is required. Run 11a first.")
}

# ============================================================
# 2. PATHS
# ============================================================

input_file <-
    "results/11_TCGA_MAP3K8/11a_MAP3K8_TCGA_TARGET_GTEx_Toil.csv"

outdir <-
    "results/11_TCGA_MAP3K8"

figdir <-
    "figures/11_TCGA_MAP3K8"

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
    stop("Missing 11a MAP3K8 Toil result.")
}

cat("\n========================================\n")
cat("11b — MAP3K8 PAN-CANCER OVEREXPRESSION\n")
cat("========================================\n\n")

# ============================================================
# 3. LOAD MAP3K8 EXPRESSION
# ============================================================

x <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_cols <- c(
    "sample",
    "expression",
    "study"
)

if (!all(required_cols %in% names(x))) {
    stop("Unexpected 11a columns.")
}

x$expression <-
    as.numeric(
        x$expression
    )

cat(
    "11a samples:",
    nrow(x),
    "\n"
)

# ============================================================
# 4. LOAD CURATED TCGA-GTEx MATCHING
# ============================================================

data(
    "tcga_gtex",
    package = "UCSCXenaShiny"
)

pheno <- tcga_gtex

required_pheno <- c(
    "sample",
    "tissue",
    "type",
    "type2"
)

if (!all(
    required_pheno %in%
        names(pheno)
)) {
    stop("Unexpected tcga_gtex structure.")
}

cat(
    "Samples in curated TCGA-GTEx map:",
    nrow(pheno),
    "\n"
)

# ============================================================
# 5. TCGA CANCER NAMES
# ============================================================

data(
    "TCGA.organ",
    package = "UCSCXenaShiny"
)

cancer_names <- TCGA.organ[
    ,
    c(
        "TCGA",
        "Detail"
    )
]

names(cancer_names) <-
    c(
        "tissue",
        "Cancer_name"
    )

# ============================================================
# 6. JOIN EXPRESSION TO CURATED PHENOTYPE
# ============================================================

dat <- merge(
    pheno,
    x[
        ,
        c(
            "sample",
            "expression"
        )
    ],
    by = "sample",
    all = FALSE
)

dat <- dat[
    dat$type2 %in%
        c(
            "tumor",
            "normal"
        ) &
    is.finite(
        dat$expression
    ),
    ,
    drop = FALSE
]

cat(
    "Samples after TCGA-GTEx matching:",
    nrow(dat),
    "\n"
)

cat("\nMapped sample classes:\n\n")

print(
    table(
        dat$type2
    )
)

cat("\nTCGA cancer codes represented:\n")

print(
    sort(
        unique(
            dat$tissue
        )
    )
)

# ============================================================
# 7. CHECK DUPLICATION WITHIN TCGA CANCER MAPPING
#
# IMPORTANT:
# The same GTEx normal sample can legitimately be mapped to
# more than one TCGA cancer type when that normal tissue is
# used as reference for multiple cancers.
#
# Therefore sample IDs do NOT need to be globally unique.
#
# What must be unique is:
#     sample + TCGA cancer code (tissue)
#
# This follows the UCSCXenaShiny implementation, which
# de-duplicates samples within each tissue/cancer group.
# ============================================================

global_dup_n <- sum(
    duplicated(
        dat$sample
    )
)

sample_tissue_key <- paste(
    dat$sample,
    dat$tissue,
    sep = "|||"
)

within_tissue_dup_n <- sum(
    duplicated(
        sample_tissue_key
    )
)

cat(
    "\nGlobally repeated sample IDs:",
    global_dup_n,
    "\n"
)

cat(
    "Repeated sample+tissue combinations:",
    within_tissue_dup_n,
    "\n"
)

# Remove only duplicated entries WITHIN the same cancer mapping.
# A GTEx sample mapped to two different cancers is retained
# once in each relevant cancer comparison.

if (within_tissue_dup_n > 0) {

    dat <- dat[
        !duplicated(
            sample_tissue_key
        ),
        ,
        drop = FALSE
    ]

    cat(
        "Removed duplicate rows within cancer mappings:",
        within_tissue_dup_n,
        "\n"
    )
}

# Final safety check

final_key <- paste(
    dat$sample,
    dat$tissue,
    sep = "|||"
)

if (anyDuplicated(final_key)) {
    stop(
        "Duplicate sample+tissue combinations remain."
    )
}

cat(
    "Final mapped sample-cancer rows:",
    nrow(dat),
    "\n"
)

# ============================================================
# 8. CANCER-BY-CANCER ANALYSIS
# ============================================================

tissues <- sort(
    unique(
        dat$tissue
    )
)

results <- list()

for (cc in tissues) {

    d <- dat[
        dat$tissue == cc,
        ,
        drop = FALSE
    ]

    tumor <- d$expression[
        d$type2 == "tumor"
    ]

    normal <- d$expression[
        d$type2 == "normal"
    ]

    n_tumor <- length(tumor)
    n_normal <- length(normal)

    if (
        n_tumor == 0 ||
        n_normal == 0
    ) {
        next
    }

    median_tumor <-
        median(
            tumor,
            na.rm = TRUE
        )

    median_normal <-
        median(
            normal,
            na.rm = TRUE
        )

    delta <-
        median_tumor -
        median_normal

    # --------------------------------------------
    # Tumor vs normal test
    # --------------------------------------------

    wt <- suppressWarnings(
        wilcox.test(
            tumor,
            normal,
            alternative = "two.sided",
            exact = FALSE
        )
    )

    # --------------------------------------------
    # High-expression threshold
    #
    # Pre-defined before seeing results:
    # > 95th percentile of corresponding normals
    # --------------------------------------------

    q95_normal <-
        as.numeric(
            quantile(
                normal,
                probs = 0.95,
                na.rm = TRUE,
                names = FALSE,
                type = 7
            )
        )

    n_high <-
        sum(
            tumor > q95_normal,
            na.rm = TRUE
        )

    pct_high <-
        100 *
        n_high /
        n_tumor

    # --------------------------------------------
    # Normal-source composition
    # --------------------------------------------

    normal_types <- d$type[
        d$type2 == "normal"
    ]

    n_normal_tcga <-
        sum(
            grepl(
                "normal_TCGA",
                normal_types,
                fixed = TRUE
            )
        )

    n_normal_gtex <-
        sum(
            grepl(
                "normal_GTEx",
                normal_types,
                fixed = TRUE
            )
        )

    normal_source <- if (
        n_normal_tcga > 0 &&
        n_normal_gtex > 0
    ) {
        "TCGA+GTEx"
    } else if (
        n_normal_tcga > 0
    ) {
        "TCGA"
    } else if (
        n_normal_gtex > 0
    ) {
        "GTEx"
    } else {
        "Other"
    }

    results[[cc]] <- data.frame(

        Cancer = cc,

        N_tumor = n_tumor,

        N_normal = n_normal,

        N_normal_TCGA =
            n_normal_tcga,

        N_normal_GTEx =
            n_normal_gtex,

        Normal_source =
            normal_source,

        Median_tumor_log2TPM =
            median_tumor,

        Median_normal_log2TPM =
            median_normal,

        Median_difference =
            delta,

        P_value =
            wt$p.value,

        Normal_95th_percentile =
            q95_normal,

        N_MAP3K8_high =
            n_high,

        Percent_MAP3K8_high =
            pct_high,

        stringsAsFactors = FALSE
    )
}

summary_table <- do.call(
    rbind,
    results
)

rownames(summary_table) <- NULL

# ============================================================
# 9. MULTIPLE TESTING
# ============================================================

summary_table$FDR <-
    p.adjust(
        summary_table$P_value,
        method = "BH"
    )

summary_table$Expression_status <-
    ifelse(
        summary_table$FDR < 0.05 &
        summary_table$Median_difference > 0,
        "OVEREXPRESSED",
        ifelse(
            summary_table$FDR < 0.05 &
            summary_table$Median_difference < 0,
            "UNDEREXPRESSED",
            "NOT_SIGNIFICANT"
        )
    )

# Reliable primary analysis:
# at least 20 tumors and 10 normals
summary_table$Reliable_comparison <-
    summary_table$N_tumor >= 20 &
    summary_table$N_normal >= 10

# Approximate pseudocount-adjusted ratio:
# because expression is log2(TPM + 0.001)
summary_table$Median_expression_ratio <-
    2^(
        summary_table$Median_difference
    )

# ============================================================
# 10. ADD FULL CANCER NAMES
# ============================================================

summary_table <- merge(
    summary_table,
    cancer_names,
    by.x = "Cancer",
    by.y = "tissue",
    all.x = TRUE
)

# Put useful columns first
summary_table <- summary_table[
    ,
    c(
        "Cancer",
        "Cancer_name",
        "N_tumor",
        "N_normal",
        "N_normal_TCGA",
        "N_normal_GTEx",
        "Normal_source",
        "Median_tumor_log2TPM",
        "Median_normal_log2TPM",
        "Median_difference",
        "Median_expression_ratio",
        "P_value",
        "FDR",
        "Expression_status",
        "Normal_95th_percentile",
        "N_MAP3K8_high",
        "Percent_MAP3K8_high",
        "Reliable_comparison"
    )
]

# Sort strongest reliable overexpression first
summary_table <- summary_table[
    order(
        !summary_table$Reliable_comparison,
        summary_table$Expression_status !=
            "OVEREXPRESSED",
        -summary_table$Percent_MAP3K8_high,
        summary_table$FDR
    ),
    ,
    drop = FALSE
]

# ============================================================
# 11. SAVE COMPLETE RESULTS
# ============================================================

write.csv(
    summary_table,
    file.path(
        outdir,
        "11b_MAP3K8_pancancer_overexpression_frequency_ALL.csv"
    ),
    row.names = FALSE
)

reliable <- summary_table[
    summary_table$Reliable_comparison,
    ,
    drop = FALSE
]

write.csv(
    reliable,
    file.path(
        outdir,
        "11b_MAP3K8_pancancer_overexpression_frequency_RELIABLE.csv"
    ),
    row.names = FALSE
)

overexpressed <- reliable[
    reliable$Expression_status ==
        "OVEREXPRESSED",
    ,
    drop = FALSE
]

write.csv(
    overexpressed,
    file.path(
        outdir,
        "11b_MAP3K8_significantly_overexpressed_cancers.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 12. CLASSIFY EVERY TCGA TUMOR AS MAP3K8-HIGH/NOT HIGH
#
# This table will be used in Step 11c for grade.
# ============================================================

thresholds <- setNames(
    summary_table$Normal_95th_percentile,
    summary_table$Cancer
)

threshold_reliable <- setNames(
    summary_table$Reliable_comparison,
    summary_table$Cancer
)

tumor_samples <- dat[
    dat$type2 == "tumor",
    ,
    drop = FALSE
]

tumor_samples$Normal_95th_percentile <-
    unname(
        thresholds[
            as.character(
                tumor_samples$tissue
            )
        ]
    )

tumor_samples$Threshold_reliable <-
    unname(
        threshold_reliable[
            as.character(
                tumor_samples$tissue
            )
        ]
    )

tumor_samples$MAP3K8_high <-
    tumor_samples$expression >
    tumor_samples$Normal_95th_percentile

tumor_samples <- merge(
    tumor_samples,
    cancer_names,
    by = "tissue",
    all.x = TRUE
)

write.csv(
    tumor_samples,
    file.path(
        outdir,
        "11b_TCGA_MAP3K8_tumor_sample_classification.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 13. CONSOLE RESULTS
# ============================================================

cat("\n========================================\n")
cat("PAN-CANCER RESULTS\n")
cat("========================================\n\n")

cat(
    "Cancer types with tumor + normal data:",
    nrow(summary_table),
    "\n"
)

cat(
    "Reliable comparisons (>=20 tumors, >=10 normals):",
    nrow(reliable),
    "\n"
)

cat(
    "Significantly MAP3K8-overexpressed cancers:",
    nrow(overexpressed),
    "\n"
)

cat(
    "Significantly MAP3K8-underexpressed cancers:",
    sum(
        reliable$Expression_status ==
            "UNDEREXPRESSED"
    ),
    "\n"
)

cat("\n========================================\n")
cat("MAP3K8-OVEREXPRESSED CANCERS\n")
cat("========================================\n\n")

if (nrow(overexpressed) == 0) {

    cat(
        "No reliable cancer type met the overexpression criterion.\n"
    )

} else {

    print(
        overexpressed[
            ,
            c(
                "Cancer",
                "Cancer_name",
                "N_tumor",
                "N_normal",
                "Median_difference",
                "Median_expression_ratio",
                "FDR",
                "N_MAP3K8_high",
                "Percent_MAP3K8_high",
                "Normal_source"
            )
        ],
        row.names = FALSE
    )
}

# ============================================================
# 14. TOP CANCERS BY HIGH-EXPRESSION FREQUENCY
# ============================================================

top_freq <- reliable[
    order(
        -reliable$Percent_MAP3K8_high
    ),
    ,
    drop = FALSE
]

cat("\n========================================\n")
cat("TOP MAP3K8-HIGH FREQUENCIES\n")
cat("========================================\n\n")

print(
    head(
        top_freq[
            ,
            c(
                "Cancer",
                "Cancer_name",
                "N_tumor",
                "N_normal",
                "Expression_status",
                "Median_difference",
                "FDR",
                "Percent_MAP3K8_high"
            )
        ],
        20
    ),
    row.names = FALSE
)

# ============================================================
# 15. SIMPLE FREQUENCY FIGURE
# ============================================================

plot_dat <- head(
    top_freq,
    20
)

if (nrow(plot_dat) > 0) {

    plot_dat <- plot_dat[
        order(
            plot_dat$Percent_MAP3K8_high
        ),
        ,
        drop = FALSE
    ]

    png(
        file.path(
            figdir,
            "11B_MAP3K8_high_frequency_top20.png"
        ),
        width = 1600,
        height = 1200,
        res = 160
    )

    par(
        mar = c(
            5,
            12,
            4,
            2
        )
    )

    barplot(
        plot_dat$Percent_MAP3K8_high,
        names.arg =
            plot_dat$Cancer,
        horiz = TRUE,
        las = 1,
        xlab =
            "% tumors above normal 95th percentile",
        main =
            "MAP3K8-high frequency across TCGA cancers"
    )

    dev.off()
}

# ============================================================
# 16. COMPLETE
# ============================================================

cat("\n========================================\n")
cat("11b COMPLETE\n")
cat("========================================\n\n")

cat(
    "Definition used:\n",
    "MAP3K8-high = tumor expression above the\n",
    "95th percentile of matched TCGA/GTEx normal tissue.\n\n",
    "Next step: integrate TCGA clinical grade and determine\n",
    "whether MAP3K8 expression/high-frequency increases\n",
    "with tumor grade.\n"
)
