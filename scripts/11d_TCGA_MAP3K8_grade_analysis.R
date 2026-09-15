# ============================================================
# 11d_TCGA_MAP3K8_grade_analysis.R
#
# MAP3K8 vs HISTOLOGIC GRADE
#
# Primary analysis:
# cancers already identified as MAP3K8-overexpressed
# AND with usable TCGA grade data:
#
#   PAAD
#   KIRC
#   LGG
#   STAD
#
# Only TCGA primary tumor samples (Code = TP) are used.
#
# Outcomes:
#   - continuous MAP3K8 expression by grade
#   - MAP3K8-high frequency by grade
#   - ordinal expression trend
# ============================================================

rm(list = ls())

# ============================================================
# 1. PACKAGE
# ============================================================

if (!requireNamespace("UCSCXenaShiny", quietly = TRUE)) {
    stop("UCSCXenaShiny is required.")
}

# ============================================================
# 2. INPUT / OUTPUT
# ============================================================

tumor_file <-
    "results/11_TCGA_MAP3K8/11b_TCGA_MAP3K8_tumor_sample_classification.csv"

over_file <-
    "results/11_TCGA_MAP3K8/11b_MAP3K8_significantly_overexpressed_cancers.csv"

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

if (!file.exists(tumor_file)) {
    stop("Missing 11b tumor classification.")
}

if (!file.exists(over_file)) {
    stop("Missing 11b overexpressed-cancer table.")
}

cat("\n========================================\n")
cat("11d — MAP3K8 vs TUMOR GRADE\n")
cat("========================================\n\n")

# ============================================================
# 3. LOAD EXPRESSION / HIGH CLASSIFICATION
# ============================================================

tumor <- read.csv(
    tumor_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_tumor <- c(
    "sample",
    "tissue",
    "expression",
    "MAP3K8_high"
)

if (!all(required_tumor %in% names(tumor))) {
    stop("Unexpected 11b tumor-classification columns.")
}

tumor$expression <-
    as.numeric(
        tumor$expression
    )

tumor$MAP3K8_high <-
    as.logical(
        tumor$MAP3K8_high
    )

tumor <- tumor[
    is.finite(tumor$expression) &
    !is.na(tumor$MAP3K8_high),
    ,
    drop = FALSE
]

# ============================================================
# 4. LOAD TCGA CLINICAL DATA
# ============================================================

data(
    "tcga_clinical_fine",
    package = "UCSCXenaShiny"
)

clin <- as.data.frame(
    tcga_clinical_fine,
    stringsAsFactors = FALSE
)

required_clin <- c(
    "Sample",
    "Cancer",
    "Code",
    "Grade"
)

if (!all(required_clin %in% names(clin))) {
    stop("Unexpected TCGA clinical columns.")
}

# Primary solid tumors only
clin <- clin[
    clin$Code == "TP" &
    clin$Grade %in%
        c(
            "G1",
            "G2",
            "G3",
            "G4"
        ),
    ,
    drop = FALSE
]

cat(
    "Clinical primary-tumor grade rows:",
    nrow(clin),
    "\n"
)

# ============================================================
# 5. REMOVE EXACT DUPLICATES / CONFLICTING GRADES
# ============================================================

clin <- unique(
    clin[
        ,
        c(
            "Sample",
            "Cancer",
            "Grade"
        )
    ]
)

key <- paste(
    clin$Sample,
    clin$Cancer,
    sep = "|||"
)

grade_n <- ave(
    clin$Grade,
    key,
    FUN = function(z) {
        length(
            unique(z)
        )
    }
)

conflicting_keys <- unique(
    key[
        grade_n > 1
    ]
)

cat(
    "Sample-cancer records with conflicting grades:",
    length(conflicting_keys),
    "\n"
)

if (length(conflicting_keys) > 0) {

    clin <- clin[
        !key %in%
            conflicting_keys,
        ,
        drop = FALSE
    ]
}

key <- paste(
    clin$Sample,
    clin$Cancer,
    sep = "|||"
)

clin <- clin[
    !duplicated(key),
    ,
    drop = FALSE
]

# ============================================================
# 6. MERGE EXPRESSION WITH CLINICAL GRADE
# ============================================================

dat <- merge(
    tumor,
    clin,
    by.x = c(
        "sample",
        "tissue"
    ),
    by.y = c(
        "Sample",
        "Cancer"
    ),
    all = FALSE
)

dat$Grade_num <-
    as.integer(
        sub(
            "^G",
            "",
            dat$Grade
        )
    )

cat(
    "Expression samples with grade:",
    nrow(dat),
    "\n\n"
)

# ============================================================
# 7. IDENTIFY OVEREXPRESSED CANCERS
# ============================================================

over <- read.csv(
    over_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

over_cancers <- unique(
    over$Cancer
)

focus <- intersect(
    over_cancers,
    unique(
        dat$tissue
    )
)

cat(
    "Overexpressed cancers with grade data:",
    paste(
        focus,
        collapse = ", "
    ),
    "\n\n"
)

# Expected:
# PAAD KIRC LGG STAD

focus_dat <- dat[
    dat$tissue %in%
        focus,
    ,
    drop = FALSE
]

write.csv(
    focus_dat,
    file.path(
        outdir,
        "11d_MAP3K8_grade_matched_samples.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 8. GRADE-LEVEL SUMMARY
# ============================================================

grade_summary_list <- list()

counter <- 1

for (cc in focus) {

    d <- focus_dat[
        focus_dat$tissue == cc,
        ,
        drop = FALSE
    ]

    grades <- sort(
        unique(
            d$Grade_num
        )
    )

    for (g in grades) {

        z <- d[
            d$Grade_num == g,
            ,
            drop = FALSE
        ]

        grade_summary_list[[counter]] <-
            data.frame(

                Cancer = cc,

                Grade =
                    paste0(
                        "G",
                        g
                    ),

                Grade_num = g,

                N =
                    nrow(z),

                Median_MAP3K8 =
                    median(
                        z$expression,
                        na.rm = TRUE
                    ),

                Mean_MAP3K8 =
                    mean(
                        z$expression,
                        na.rm = TRUE
                    ),

                IQR_MAP3K8 =
                    IQR(
                        z$expression,
                        na.rm = TRUE
                    ),

                N_MAP3K8_high =
                    sum(
                        z$MAP3K8_high,
                        na.rm = TRUE
                    ),

                Percent_MAP3K8_high =
                    100 *
                    mean(
                        z$MAP3K8_high,
                        na.rm = TRUE
                    ),

                stringsAsFactors = FALSE
            )

        counter <- counter + 1
    }
}

grade_summary <- do.call(
    rbind,
    grade_summary_list
)

write.csv(
    grade_summary,
    file.path(
        outdir,
        "11d_MAP3K8_expression_frequency_by_grade.csv"
    ),
    row.names = FALSE
)

cat("========================================\n")
cat("MAP3K8 BY GRADE\n")
cat("========================================\n\n")

print(
    grade_summary,
    row.names = FALSE
)

# ============================================================
# 9. WITHIN-CANCER STATISTICAL TESTS
# ============================================================

test_results <- list()

counter <- 1

for (cc in focus) {

    d <- focus_dat[
        focus_dat$tissue == cc,
        ,
        drop = FALSE
    ]

    grade_counts <- table(
        d$Grade_num
    )

    # Require at least two usable grade groups
    usable_grades <- as.integer(
        names(
            grade_counts[
                grade_counts >= 10
            ]
        )
    )

    test_d <- d[
        d$Grade_num %in%
            usable_grades,
        ,
        drop = FALSE
    ]

    if (
        length(
            unique(
                test_d$Grade_num
            )
        ) < 2
    ) {
        next
    }

    # --------------------------------------------
    # Spearman ordinal trend
    # --------------------------------------------

    sp <- suppressWarnings(
        cor.test(
            test_d$Grade_num,
            test_d$expression,
            method = "spearman",
            exact = FALSE
        )
    )

    # --------------------------------------------
    # Kruskal-Wallis across grade groups
    # --------------------------------------------

    kw <- kruskal.test(
        expression ~
            factor(
                Grade_num
            ),
        data = test_d
    )

    # --------------------------------------------
    # Logistic trend in MAP3K8-high frequency
    #
    # Only calculate if enough high/low cases.
    # --------------------------------------------

    n_high <- sum(
        test_d$MAP3K8_high
    )

    n_low <- sum(
        !test_d$MAP3K8_high
    )

    logistic_or <- NA_real_
    logistic_low <- NA_real_
    logistic_high <- NA_real_
    logistic_p <- NA_real_

    if (
        n_high >= 10 &&
        n_low >= 10
    ) {

        fit <- suppressWarnings(
            glm(
                MAP3K8_high ~
                    Grade_num,
                data = test_d,
                family = binomial()
            )
        )

        co <- summary(
            fit
        )$coefficients

        beta <-
            co[
                "Grade_num",
                "Estimate"
            ]

        se <-
            co[
                "Grade_num",
                "Std. Error"
            ]

        logistic_or <-
            exp(beta)

        logistic_low <-
            exp(
                beta -
                1.96 * se
            )

        logistic_high <-
            exp(
                beta +
                1.96 * se
            )

        logistic_p <-
            co[
                "Grade_num",
                "Pr(>|z|)"
            ]
    }

    test_results[[counter]] <-
        data.frame(

            Cancer = cc,

            N =
                nrow(test_d),

            Grades_tested =
                paste(
                    paste0(
                        "G",
                        sort(
                            unique(
                                test_d$Grade_num
                            )
                        )
                    ),
                    collapse = ","
                ),

            Spearman_rho =
                unname(
                    sp$estimate
                ),

            Spearman_P =
                sp$p.value,

            Kruskal_P =
                kw$p.value,

            High_frequency_OR_per_grade =
                logistic_or,

            High_frequency_OR_95CI_low =
                logistic_low,

            High_frequency_OR_95CI_high =
                logistic_high,

            High_frequency_trend_P =
                logistic_p,

            stringsAsFactors = FALSE
        )

    counter <- counter + 1
}

tests <- do.call(
    rbind,
    test_results
)

tests$Spearman_FDR <-
    p.adjust(
        tests$Spearman_P,
        method = "BH"
    )

tests$Kruskal_FDR <-
    p.adjust(
        tests$Kruskal_P,
        method = "BH"
    )

tests$High_frequency_trend_FDR <-
    p.adjust(
        tests$High_frequency_trend_P,
        method = "BH"
    )

write.csv(
    tests,
    file.path(
        outdir,
        "11d_MAP3K8_grade_statistical_tests.csv"
    ),
    row.names = FALSE
)

cat("\n========================================\n")
cat("GRADE ASSOCIATION TESTS\n")
cat("========================================\n\n")

print(
    tests,
    row.names = FALSE
)

# ============================================================
# 10. FIGURES FOR EACH CANCER
# ============================================================

for (cc in focus) {

    d <- focus_dat[
        focus_dat$tissue == cc,
        ,
        drop = FALSE
    ]

    d$Grade <- factor(
        d$Grade,
        levels = c(
            "G1",
            "G2",
            "G3",
            "G4"
        )
    )

    png(
        file.path(
            figdir,
            paste0(
                "11D_MAP3K8_grade_",
                cc,
                ".png"
            )
        ),
        width = 1200,
        height = 1000,
        res = 160
    )

    boxplot(
        expression ~ Grade,
        data = d,
        outline = FALSE,
        xlab = "Tumor grade",
        ylab = "MAP3K8 log2(TPM + 0.001)",
        main =
            paste0(
                cc,
                ": MAP3K8 expression by grade"
            )
    )

    stripchart(
        expression ~ Grade,
        data = d,
        vertical = TRUE,
        method = "jitter",
        pch = 16,
        cex = 0.45,
        add = TRUE
    )

    dev.off()
}

# ============================================================
# 11. SEPARATE HISTORICAL GLIOMA GRADE ANALYSIS
#
# TCGA LGG:
#   WHO grade II and III
#
# TCGA GBM:
#   historically WHO grade IV
#
# This is reported separately and must NOT be interpreted
# using modern molecular WHO classification.
# ============================================================

glioma <- tumor[
    tumor$tissue %in%
        c(
            "LGG",
            "GBM"
        ),
    ,
    drop = FALSE
]

# LGG exact clinical G2/G3
lgg_clin <- clin[
    clin$Cancer == "LGG" &
    clin$Grade %in%
        c(
            "G2",
            "G3"
        ),
    ,
    drop = FALSE
]

lgg <- merge(
    glioma[
        glioma$tissue ==
            "LGG",
        ,
        drop = FALSE
    ],
    lgg_clin,
    by.x = c(
        "sample",
        "tissue"
    ),
    by.y = c(
        "Sample",
        "Cancer"
    ),
    all = FALSE
)

lgg$Historical_grade_num <-
    as.integer(
        sub(
            "^G",
            "",
            lgg$Grade
        )
    )

# GBM treated as historical WHO grade IV
gbm <- glioma[
    glioma$tissue ==
        "GBM",
    ,
    drop = FALSE
]

gbm$Historical_grade_num <-
    4L

glioma_grade <- rbind(
    data.frame(
        sample =
            lgg$sample,

        tissue =
            lgg$tissue,

        expression =
            lgg$expression,

        MAP3K8_high =
            lgg$MAP3K8_high,

        Historical_grade_num =
            lgg$Historical_grade_num
    ),

    data.frame(
        sample =
            gbm$sample,

        tissue =
            gbm$tissue,

        expression =
            gbm$expression,

        MAP3K8_high =
            gbm$MAP3K8_high,

        Historical_grade_num =
            gbm$Historical_grade_num
    )
)

glioma_summary <- do.call(
    rbind,
    lapply(
        sort(
            unique(
                glioma_grade$Historical_grade_num
            )
        ),
        function(g) {

            z <- glioma_grade[
                glioma_grade$Historical_grade_num ==
                    g,
                ,
                drop = FALSE
            ]

            data.frame(
                Historical_grade =
                    paste0(
                        "G",
                        g
                    ),

                N =
                    nrow(z),

                Median_MAP3K8 =
                    median(
                        z$expression,
                        na.rm = TRUE
                    ),

                N_MAP3K8_high =
                    sum(
                        z$MAP3K8_high,
                        na.rm = TRUE
                    ),

                Percent_MAP3K8_high =
                    100 *
                    mean(
                        z$MAP3K8_high,
                        na.rm = TRUE
                    )
            )
        }
    )
)

glioma_sp <- cor.test(
    glioma_grade$Historical_grade_num,
    glioma_grade$expression,
    method = "spearman",
    exact = FALSE
)

glioma_kw <- kruskal.test(
    expression ~
        factor(
            Historical_grade_num
        ),
    data = glioma_grade
)

write.csv(
    glioma_summary,
    file.path(
        outdir,
        "11d_MAP3K8_historical_glioma_G2_G3_G4.csv"
    ),
    row.names = FALSE
)

cat("\n========================================\n")
cat("HISTORICAL GLIOMA GRADES II-IV\n")
cat("========================================\n\n")

print(
    glioma_summary,
    row.names = FALSE
)

cat(
    "\nGlioma grade-expression Spearman rho:",
    round(
        unname(
            glioma_sp$estimate
        ),
        4
    ),
    "\n"
)

cat(
    "Spearman P:",
    signif(
        glioma_sp$p.value,
        5
    ),
    "\n"
)

cat(
    "Kruskal-Wallis P:",
    signif(
        glioma_kw$p.value,
        5
    ),
    "\n"
)

# ============================================================
# 12. COMPLETE
# ============================================================

cat("\n========================================\n")
cat("11d COMPLETE\n")
cat("========================================\n\n")

cat(
    "Interpretation rules:\n\n",

    "Positive Spearman rho with FDR<0.05:\n",
    "MAP3K8 expression increases with grade.\n\n",

    "High-frequency OR >1 with FDR<0.05:\n",
    "odds of MAP3K8-high status increase with grade.\n\n",

    "Grade comparisons are cancer-specific.\n",
    "Do not pool unrelated cancer grading systems.\n"
)
