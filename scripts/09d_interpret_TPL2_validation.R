# ============================================================
# 09d_interpret_TPL2_validation.R
#
# PURPOSE
#
# 1. Identify the robust persistent-UP genes contributing
#    to the TPL2 CAMERA signal.
#
# 2. Identify exact persistent genes significantly reversed
#    by TPL2 inhibition.
#
# 3. Positive-control comparison:
#    GSE85246 Day-1 acute LPS response
#             vs
#    GSE84161 acute LPS response
#
# This step interprets 09c rather than adding another dataset.
# ============================================================

rm(list = ls())

# ============================================================
# 1. INPUTS / OUTPUTS
# ============================================================

integrated_file <-
    "results/09_MAP3K8_validation/09c_GSE85246_GSE84161_integrated_gene_table.csv"

day1_file <-
    "results/02_initial_LPS/02_Day1_LPS_vs_RPMI_all_genes.csv"

acute84161_file <-
    "results/09_MAP3K8_validation/09b_LPS_vs_Untreated.csv"

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
    integrated_file,
    day1_file,
    acute84161_file
)

missing_files <-
    required_files[
        !file.exists(required_files)
    ]

if (length(missing_files) > 0) {

    stop(
        paste(
            "Missing files:",
            paste(
                missing_files,
                collapse = "\n"
            )
        )
    )
}

cat("\n========================================\n")
cat("09d — INTERPRETING TPL2 VALIDATION\n")
cat("========================================\n\n")

# ============================================================
# 2. LOAD 09c INTEGRATED TABLE
# ============================================================

integrated <-
    read.csv(
        integrated_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

# ============================================================
# 3. THE ROBUST MEMORY-UP GENE SET
# ============================================================

memory_up <-
    integrated[
        integrated$Memory_robust &
        integrated$Memory_log2FC > 0,
        ,
        drop = FALSE
    ]

memory_up <-
    memory_up[
        order(
            memory_up$TPL2i_FDR,
            -abs(memory_up$TPL2i_log2FC)
        ),
        ,
        drop = FALSE
    ]

cat(
    "Robust persistent-UP genes shared with GSE84161:",
    nrow(memory_up),
    "\n\n"
)

print(
    memory_up[
        ,
        c(
            "SYMBOL",
            "Memory_log2FC",
            "Memory_FDR",
            "AcuteLPS_log2FC",
            "TPL2i_log2FC",
            "TPL2i_FDR",
            "MEKi_log2FC",
            "MEKi_FDR"
        )
    ],
    row.names = FALSE
)

write.csv(
    memory_up[
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
    file.path(
        outdir,
        "09d_robust_memory_UP_13_genes.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 4. SIGNIFICANT TPL2 REVERSALS
# ============================================================

tpl2_reversal <-
    integrated[
        integrated$Memory_robust &
        integrated$TPL2i_FDR < 0.05 &
        sign(
            integrated$Memory_log2FC
        ) !=
        sign(
            integrated$TPL2i_log2FC
        ),
        ,
        drop = FALSE
    ]

tpl2_reversal <-
    tpl2_reversal[
        order(
            tpl2_reversal$TPL2i_FDR
        ),
        ,
        drop = FALSE
    ]

cat("\n========================================\n")
cat("SIGNIFICANT TPL2 REVERSALS\n")
cat("========================================\n\n")

cat(
    "Number of significant reversals:",
    nrow(tpl2_reversal),
    "\n\n"
)

print(
    tpl2_reversal[
        ,
        c(
            "SYMBOL",
            "Memory_log2FC",
            "Memory_FDR",
            "TPL2i_log2FC",
            "TPL2i_FDR",
            "MEKi_log2FC",
            "MEKi_FDR"
        )
    ],
    row.names = FALSE
)

write.csv(
    tpl2_reversal,
    file.path(
        outdir,
        "09d_significant_TPL2_reversal_genes.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 5. SIGNIFICANT SAME-DIRECTION TPL2 EFFECTS
# ============================================================

tpl2_same <-
    integrated[
        integrated$Memory_robust &
        integrated$TPL2i_FDR < 0.05 &
        sign(
            integrated$Memory_log2FC
        ) ==
        sign(
            integrated$TPL2i_log2FC
        ),
        ,
        drop = FALSE
    ]

write.csv(
    tpl2_same,
    file.path(
        outdir,
        "09d_significant_TPL2_same_direction_genes.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 6. LOAD GSE85246 DAY-1 ACUTE LPS RESULT
# ============================================================

day1 <-
    read.csv(
        day1_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

required_day1 <- c(
    "SYMBOL",
    "log2FC_Day1_LPS_vs_RPMI",
    "AveExpr",
    "P.Value",
    "FDR_genome"
)

if (!all(
    required_day1 %in%
        names(day1)
)) {

    stop(
        "Unexpected GSE85246 Day-1 columns."
    )
}

day1 <-
    day1[
        !is.na(day1$SYMBOL) &
        day1$SYMBOL != "",
        ,
        drop = FALSE
    ]

# ============================================================
# 7. COLLAPSE DUPLICATED HUMAN SYMBOLS
#
# Choose highest average-expression row.
# This does not use P value or fold-change direction.
# ============================================================

day1 <-
    day1[
        order(
            day1$SYMBOL,
            -day1$AveExpr
        ),
        ,
        drop = FALSE
    ]

day1 <-
    day1[
        !duplicated(day1$SYMBOL),
        ,
        drop = FALSE
    ]

names(day1)[
    names(day1) ==
        "log2FC_Day1_LPS_vs_RPMI"
] <- "GSE85246_Day1_log2FC"

names(day1)[
    names(day1) ==
        "P.Value"
] <- "GSE85246_Day1_P"

names(day1)[
    names(day1) ==
        "FDR_genome"
] <- "GSE85246_Day1_FDR"

# ============================================================
# 8. LOAD GSE84161 ACUTE LPS RESPONSE
# ============================================================

acute84161 <-
    read.csv(
        acute84161_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

acute84161 <-
    acute84161[
        ,
        c(
            "SYMBOL",
            "log2FC",
            "P.Value",
            "FDR"
        ),
        drop = FALSE
    ]

names(acute84161)[2:4] <-
    c(
        "GSE84161_LPS_log2FC",
        "GSE84161_LPS_P",
        "GSE84161_LPS_FDR"
    )

# ============================================================
# 9. MERGE ACUTE LPS RESULTS
# ============================================================

acute_compare <-
    merge(
        day1[
            ,
            c(
                "SYMBOL",
                "GSE85246_Day1_log2FC",
                "GSE85246_Day1_P",
                "GSE85246_Day1_FDR"
            )
        ],
        acute84161,
        by = "SYMBOL",
        all = FALSE
    )

cat("\n========================================\n")
cat("ACUTE LPS CROSS-STUDY COMPARISON\n")
cat("========================================\n\n")

cat(
    "Shared genes:",
    nrow(acute_compare),
    "\n"
)

# ============================================================
# 10. GLOBAL SPEARMAN CORRELATION
# ============================================================

ok <-
    is.finite(
        acute_compare$GSE85246_Day1_log2FC
    ) &
    is.finite(
        acute_compare$GSE84161_LPS_log2FC
    )

cor_all <-
    cor.test(
        acute_compare$GSE85246_Day1_log2FC[ok],
        acute_compare$GSE84161_LPS_log2FC[ok],
        method = "spearman",
        exact = FALSE
    )

# ============================================================
# 11. STRONG GSE84161 ACUTE-LPS REFERENCE GENES
#
# GSE84161 has 5 paired donors and much greater power.
# We therefore define the positive-control LPS signature there.
# ============================================================

acute_compare$GSE84161_strong_LPS <-
    acute_compare$GSE84161_LPS_FDR < 0.05 &
    abs(
        acute_compare$GSE84161_LPS_log2FC
    ) >= 1

strong <-
    acute_compare[
        acute_compare$GSE84161_strong_LPS,
        ,
        drop = FALSE
    ]

strong$Direction_concordant <-
    sign(
        strong$GSE85246_Day1_log2FC
    ) ==
    sign(
        strong$GSE84161_LPS_log2FC
    )

n_concordant <-
    sum(
        strong$Direction_concordant,
        na.rm = TRUE
    )

n_strong <-
    nrow(strong)

concordance_fraction <-
    n_concordant /
    n_strong

binom_result <-
    binom.test(
        n_concordant,
        n_strong,
        p = 0.5,
        alternative = "greater"
    )

# ============================================================
# 12. STRONG-SIGNATURE SPEARMAN CORRELATION
# ============================================================

cor_strong <-
    cor.test(
        strong$GSE85246_Day1_log2FC,
        strong$GSE84161_LPS_log2FC,
        method = "spearman",
        exact = FALSE
    )

acute_summary <-
    data.frame(
        Test = c(
            "All_shared_genes_Spearman",
            "Strong_GSE84161_LPS_genes_Spearman",
            "Strong_GSE84161_LPS_direction_concordance"
        ),

        N = c(
            sum(ok),
            n_strong,
            n_strong
        ),

        Statistic = c(
            unname(
                cor_all$estimate
            ),
            unname(
                cor_strong$estimate
            ),
            concordance_fraction
        ),

        P_value = c(
            cor_all$p.value,
            cor_strong$p.value,
            binom_result$p.value
        )
    )

cat(
    "\nGlobal Day1-LPS Spearman rho:",
    round(
        unname(
            cor_all$estimate
        ),
        4
    ),
    "\n"
)

cat(
    "Global correlation P:",
    signif(
        cor_all$p.value,
        4
    ),
    "\n\n"
)

cat(
    "Strong GSE84161 acute-LPS genes:",
    n_strong,
    "\n"
)

cat(
    "Direction concordant:",
    n_concordant,
    "/",
    n_strong,
    "=",
    round(
        100 *
        concordance_fraction,
        1
    ),
    "%\n"
)

cat(
    "Concordance binomial P:",
    signif(
        binom_result$p.value,
        4
    ),
    "\n"
)

cat(
    "Strong-gene Spearman rho:",
    round(
        unname(
            cor_strong$estimate
        ),
        4
    ),
    "\n"
)

cat(
    "Strong-gene correlation P:",
    signif(
        cor_strong$p.value,
        4
    ),
    "\n"
)

write.csv(
    acute_summary,
    file.path(
        outdir,
        "09d_acute_LPS_crossstudy_summary.csv"
    ),
    row.names = FALSE
)

write.csv(
    strong[
        order(
            -abs(
                strong$GSE84161_LPS_log2FC
            )
        ),
        ,
        drop = FALSE
    ],
    file.path(
        outdir,
        "09d_strong_acute_LPS_crossstudy_genes.csv"
    ),
    row.names = FALSE
)

# ============================================================
# 13. ACUTE LPS SCATTERPLOT
# ============================================================

png(
    file.path(
        figdir,
        "09D1_GSE85246_vs_GSE84161_acute_LPS.png"
    ),
    width = 1400,
    height = 1200,
    res = 160
)

plot(
    acute_compare$GSE85246_Day1_log2FC,
    acute_compare$GSE84161_LPS_log2FC,
    pch = 16,
    cex = 0.4,
    xlab =
        "GSE85246 Day-1 LPS vs RPMI log2FC",
    ylab =
        "GSE84161 LPS vs untreated log2FC",
    main =
        "Independent human acute-LPS response comparison"
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
    strong$GSE85246_Day1_log2FC,
    strong$GSE84161_LPS_log2FC,
    pch = 1,
    cex = 0.8
)

dev.off()

# ============================================================
# 14. FINAL SUMMARY FILE
# ============================================================

sink(
    file.path(
        outdir,
        "09d_interpretation_summary.txt"
    )
)

cat(
    "09d MAP3K8/TPL2 validation interpretation\n\n"
)

cat(
    "Robust memory-UP genes:",
    nrow(memory_up),
    "\n"
)

cat(
    "Significant TPL2 reversal genes:",
    nrow(tpl2_reversal),
    "\n\n"
)

cat(
    "TPL2 reversal genes:\n"
)

print(
    tpl2_reversal[
        ,
        c(
            "SYMBOL",
            "Memory_log2FC",
            "TPL2i_log2FC",
            "TPL2i_FDR"
        )
    ],
    row.names = FALSE
)

cat(
    "\nAcute LPS cross-study validation:\n"
)

print(
    acute_summary,
    row.names = FALSE
)

sink()

cat("\n========================================\n")
cat("09d COMPLETE\n")
cat("========================================\n\n")

cat(
    "Next decision depends on:\n",
    "1. identity of TPL2-reversed memory genes\n",
    "2. acute LPS concordance between studies\n"
)
