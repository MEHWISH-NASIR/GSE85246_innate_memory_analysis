# ============================================================
# 10d_conserved_TPL2_vs_memory.R
#
# FINAL TRANSCRIPTOMIC INTEGRATION
#
# Question:
# Is the GSE85246 Day-6 persistent-memory program connected
# to a TPL2/MAP3K8-dependent transcriptional program that is
# supported independently by:
#
#   1. human TPL2 pharmacological inhibition (GSE84161)
#   2. mouse Map3k8 kinase-inactive genetics (GSE116220)
#
# Primary hypothesis:
# Persistent MEMORY-UP genes should be suppressed when
# TPL2/MAP3K8 function is lost.
# ============================================================

rm(list = ls())

# ============================================================
# 1. INPUTS
# ============================================================

memory_file <-
    "results/03_memory_kinases/03_Day6_LPS_vs_RPMI_all_genes.csv"

human_file <-
    "results/10_MAP3K8_genetic_validation/10c_GSE84161_TPL2_LPS_interaction.csv"

mouse_files <- c(
    "0.5hr" =
        "results/10_MAP3K8_genetic_validation/10c_human_mouse_integrated_0.5hr.csv",

    "1hr" =
        "results/10_MAP3K8_genetic_validation/10c_human_mouse_integrated_1hr.csv",

    "2hr" =
        "results/10_MAP3K8_genetic_validation/10c_human_mouse_integrated_2hr.csv"
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

required <- c(
    memory_file,
    human_file,
    mouse_files
)

if (!all(file.exists(required))) {

    print(
        required[
            !file.exists(required)
        ]
    )

    stop("Required input file missing.")
}

cat("\n========================================\n")
cat("10d — CONSERVED TPL2 PROGRAM vs MEMORY\n")
cat("========================================\n\n")

# ============================================================
# 2. LOAD GSE85246 MEMORY
# ============================================================

memory <- read.csv(
    memory_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_memory <- c(
    "SYMBOL",
    "log2FC_Day6_LPS_vs_RPMI",
    "FDR_genome",
    "AveExpr"
)

if (!all(
    required_memory %in%
        names(memory)
)) {
    stop("Unexpected memory-result columns.")
}

memory <- memory[
    !is.na(memory$SYMBOL) &
    memory$SYMBOL != "",
    ,
    drop = FALSE
]

# Condition-blind duplicate collapse:
# retain highest average-expression row
memory <- memory[
    order(
        memory$SYMBOL,
        -memory$AveExpr
    ),
    ,
    drop = FALSE
]

memory <- memory[
    !duplicated(memory$SYMBOL),
    ,
    drop = FALSE
]

names(memory)[
    names(memory) ==
        "log2FC_Day6_LPS_vs_RPMI"
] <- "Memory_effect"

names(memory)[
    names(memory) ==
        "FDR_genome"
] <- "Memory_FDR"

memory$Memory_robust <-
    memory$Memory_FDR < 0.05 &
    abs(memory$Memory_effect) >= 1

memory$Memory_UP <-
    memory$Memory_robust &
    memory$Memory_effect > 0

memory$Memory_DOWN <-
    memory$Memory_robust &
    memory$Memory_effect < 0

cat(
    "Robust memory genes:",
    sum(memory$Memory_robust),
    "\n"
)

cat(
    "Memory UP:",
    sum(memory$Memory_UP),
    "\n"
)

cat(
    "Memory DOWN:",
    sum(memory$Memory_DOWN),
    "\n\n"
)

# ============================================================
# 3. LOAD HUMAN TPL2 INTERACTION
# ============================================================

human <- read.csv(
    human_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_human <- c(
    "SYMBOL",
    "Human_effect",
    "Human_FDR"
)

if (!all(
    required_human %in%
        names(human)
)) {
    stop("Unexpected human TPL2 columns.")
}

# ============================================================
# 4. HUMAN MEMORY INTEGRATION
# ============================================================

mh <- merge(
    memory[
        ,
        c(
            "SYMBOL",
            "Memory_effect",
            "Memory_FDR",
            "Memory_robust",
            "Memory_UP",
            "Memory_DOWN"
        )
    ],
    human[
        ,
        c(
            "SYMBOL",
            "Human_effect",
            "Human_FDR"
        )
    ],
    by = "SYMBOL",
    all = FALSE
)

mh$Human_reverses_memory <-
    sign(mh$Memory_effect) !=
    sign(mh$Human_effect)

mh$Human_TPL2_significant <-
    mh$Human_FDR < 0.05

# ============================================================
# 5. LOAD EACH MOUSE TIMEPOINT
# ============================================================

results <- list()
summary_rows <- list()

for (tm in names(mouse_files)) {

    mouse <- read.csv(
        mouse_files[[tm]],
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    required_mouse <- c(
        "Human_SYMBOL",
        "Mouse_SYMBOL",
        "Mouse_effect",
        "Mouse_FDR"
    )

    if (!all(
        required_mouse %in%
            names(mouse)
    )) {
        stop(
            paste(
                "Unexpected mouse-integrated columns at",
                tm
            )
        )
    }

    mm <- merge(
        mh,
        mouse[
            ,
            c(
                "Human_SYMBOL",
                "Mouse_SYMBOL",
                "Mouse_effect",
                "Mouse_FDR"
            )
        ],
        by.x = "SYMBOL",
        by.y = "Human_SYMBOL",
        all = FALSE
    )

    mm$Mouse_reverses_memory <-
        sign(mm$Memory_effect) !=
        sign(mm$Mouse_effect)

    mm$Mouse_significant <-
        mm$Mouse_FDR < 0.05

    mm$Both_perturbations_reverse <-
        mm$Human_reverses_memory &
        mm$Mouse_reverses_memory

    mm$Both_perturbations_significant <-
        mm$Human_TPL2_significant &
        mm$Mouse_significant

    mm$Conserved_significant_reversal <-
        mm$Both_perturbations_reverse &
        mm$Both_perturbations_significant

    mm$Mouse_time <- tm

    results[[tm]] <- mm

    # --------------------------------------------
    # Robust memory-UP branch
    # --------------------------------------------

    up <- mm$Memory_UP

    up_human_down <-
        up &
        mm$Human_effect < 0

    up_mouse_down <-
        up &
        mm$Mouse_effect < 0

    up_both_down <-
        up &
        mm$Human_effect < 0 &
        mm$Mouse_effect < 0

    up_both_sig_down <-
        up &
        mm$Human_effect < 0 &
        mm$Mouse_effect < 0 &
        mm$Human_FDR < 0.05 &
        mm$Mouse_FDR < 0.05

    # --------------------------------------------
    # Robust memory-DOWN branch
    # --------------------------------------------

    down <- mm$Memory_DOWN

    down_human_up <-
        down &
        mm$Human_effect > 0

    down_mouse_up <-
        down &
        mm$Mouse_effect > 0

    down_both_up <-
        down &
        mm$Human_effect > 0 &
        mm$Mouse_effect > 0

    down_both_sig_up <-
        down &
        mm$Human_effect > 0 &
        mm$Mouse_effect > 0 &
        mm$Human_FDR < 0.05 &
        mm$Mouse_FDR < 0.05

    summary_rows[[tm]] <- data.frame(
        Mouse_time = tm,

        Shared_genes =
            nrow(mm),

        Robust_memory_genes =
            sum(mm$Memory_robust),

        Memory_UP =
            sum(up),

        Memory_UP_human_down =
            sum(up_human_down),

        Memory_UP_mouse_down =
            sum(up_mouse_down),

        Memory_UP_both_down =
            sum(up_both_down),

        Memory_UP_both_significant_down =
            sum(up_both_sig_down),

        Memory_DOWN =
            sum(down),

        Memory_DOWN_human_up =
            sum(down_human_up),

        Memory_DOWN_mouse_up =
            sum(down_mouse_up),

        Memory_DOWN_both_up =
            sum(down_both_up),

        Memory_DOWN_both_significant_up =
            sum(down_both_sig_up),

        stringsAsFactors = FALSE
    )

    write.csv(
        mm,
        file.path(
            outdir,
            paste0(
                "10d_memory_TPL2_integrated_",
                tm,
                ".csv"
            )
        ),
        row.names = FALSE
    )
}

summary_table <- do.call(
    rbind,
    summary_rows
)

write.csv(
    summary_table,
    file.path(
        outdir,
        "10d_memory_conserved_TPL2_summary.csv"
    ),
    row.names = FALSE
)

cat("========================================\n")
cat("DESCRIPTIVE INTEGRATION\n")
cat("========================================\n\n")

print(
    summary_table,
    row.names = FALSE
)

# ============================================================
# 6. RANK-BASED MEMORY GENE-SET TESTS
#
# We test the pre-defined GSE85246 robust memory sets against:
#
#   Human TPL2 interaction t-statistic
#   Mouse D270A interaction t-statistic
#
# Supportive:
#
# Memory_UP   -> DOWN under loss of TPL2
# Memory_DOWN -> UP under loss of TPL2
# ============================================================

memory_up_symbols <- memory$SYMBOL[
    memory$Memory_UP
]

memory_down_symbols <- memory$SYMBOL[
    memory$Memory_DOWN
]

genesets <- list(
    Memory_UP =
        memory_up_symbols,

    Memory_DOWN =
        memory_down_symbols
)

# ------------------------------------------------------------
# Human cameraPR
# ------------------------------------------------------------

human_stat <- human$Human_effect
names(human_stat) <- human$SYMBOL

human_index <- lapply(
    genesets,
    function(gs) {

        which(
            names(human_stat) %in%
                gs
        )
    }
)

human_index <- human_index[
    vapply(
        human_index,
        length,
        integer(1)
    ) >= 2
]

human_cam <- limma::cameraPR(
    statistic = human_stat,
    index = human_index,
    sort = FALSE
)

human_cam <- data.frame(
    Dataset = "Human_GSE84161",
    Time = "6h",
    GeneSet = rownames(human_cam),
    human_cam,
    row.names = NULL,
    check.names = FALSE
)

# ------------------------------------------------------------
# Mouse cameraPR by time
# ------------------------------------------------------------

mouse_camera <- list()

for (tm in names(results)) {

    x <- results[[tm]]

    mouse_stat <-
        x$Mouse_effect

    names(mouse_stat) <-
        x$SYMBOL

    index <- lapply(
        genesets,
        function(gs) {

            which(
                names(mouse_stat) %in%
                    gs
            )
        }
    )

    index <- index[
        vapply(
            index,
            length,
            integer(1)
        ) >= 2
    ]

    cam <- limma::cameraPR(
        statistic = mouse_stat,
        index = index,
        sort = FALSE
    )

    mouse_camera[[tm]] <- data.frame(
        Dataset = "Mouse_GSE116220",
        Time = tm,
        GeneSet = rownames(cam),
        cam,
        row.names = NULL,
        check.names = FALSE
    )
}

camera_table <- do.call(
    rbind,
    c(
        list(human_cam),
        mouse_camera
    )
)

camera_table$Global_FDR <-
    p.adjust(
        camera_table$PValue,
        method = "BH"
    )

write.csv(
    camera_table,
    file.path(
        outdir,
        "10d_memory_gene_set_tests.csv"
    ),
    row.names = FALSE
)

cat("\n========================================\n")
cat("MEMORY GENE-SET TESTS\n")
cat("========================================\n\n")

print(
    camera_table,
    row.names = FALSE
)

# ============================================================
# 7. IDENTIFY CONSERVED MEMORY-UP CANDIDATES
#
# Primary:
# memory UP
# human TPL2 effect DOWN
# mouse D270A effect DOWN
#
# Rank across mouse time points.
# ============================================================

candidate_rows <- list()

for (tm in names(results)) {

    x <- results[[tm]]

    hit <- x[
        x$Memory_UP &
        x$Human_effect < 0 &
        x$Mouse_effect < 0,
        ,
        drop = FALSE
    ]

    if (nrow(hit) > 0) {

        candidate_rows[[tm]] <- hit[
            ,
            c(
                "SYMBOL",
                "Mouse_SYMBOL",
                "Memory_effect",
                "Memory_FDR",
                "Human_effect",
                "Human_FDR",
                "Mouse_effect",
                "Mouse_FDR",
                "Mouse_time"
            )
        ]
    }
}

candidates <- do.call(
    rbind,
    candidate_rows
)

if (!is.null(candidates)) {

    candidates$Human_significant <-
        candidates$Human_FDR < 0.05

    candidates$Mouse_significant <-
        candidates$Mouse_FDR < 0.05

    candidates$Both_significant <-
        candidates$Human_significant &
        candidates$Mouse_significant

    candidate_count <- aggregate(
        Mouse_time ~ SYMBOL,
        data = candidates,
        FUN = length
    )

    names(candidate_count)[2] <-
        "N_mouse_timepoints_concordant"

    candidates <- merge(
        candidates,
        candidate_count,
        by = "SYMBOL",
        all.x = TRUE
    )

    candidates <- candidates[
        order(
            -candidates$Both_significant,
            -candidates$N_mouse_timepoints_concordant,
            candidates$Human_FDR,
            candidates$Mouse_FDR
        ),
        ,
        drop = FALSE
    ]

    write.csv(
        candidates,
        file.path(
            outdir,
            "10d_conserved_TPL2_supported_memory_UP_genes.csv"
        ),
        row.names = FALSE
    )

    cat("\n========================================\n")
    cat("CONSERVED MEMORY-UP CANDIDATES\n")
    cat("========================================\n\n")

    print(
        candidates,
        row.names = FALSE
    )
}

# ============================================================
# 8. MAP3K8 ITSELF
# ============================================================

cat("\n========================================\n")
cat("MAP3K8 STATUS\n")
cat("========================================\n\n")

m <- memory[
    memory$SYMBOL == "MAP3K8",
    c(
        "SYMBOL",
        "Memory_effect",
        "Memory_FDR"
    ),
    drop = FALSE
]

print(
    m,
    row.names = FALSE
)

cat(
    "\nRemember: MAP3K8 itself did not meet genome-wide\n",
    "Day-6 FDR<0.05; it was prioritized through the\n",
    "kinase-family analysis and multi-layer evidence.\n"
)

# ============================================================
# 9. COMPLETE
# ============================================================

cat("\n========================================\n")
cat("10d COMPLETE\n")
cat("========================================\n\n")

cat(
    "Primary interpretation criterion:\n\n",
    "If the GSE85246 Memory_UP set trends DOWN after\n",
    "TPL2 inhibition in human cells AND after Map3k8\n",
    "catalytic inactivation in mouse macrophages,\n",
    "this supports a conserved MAP3K8/TPL2-dependent\n",
    "component of the persistent memory-UP program.\n\n",
    "This still does NOT prove that MAP3K8 maintains\n",
    "the Day-6 state after washout; that requires a\n",
    "post-washout perturbation experiment.\n"
)
