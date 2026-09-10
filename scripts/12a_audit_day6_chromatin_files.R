# ============================================================
# 12a — AUDIT REQUIRED DAY-6 CHROMATIN FILES
#
# Goal:
# Identify exact supplementary BigWig files for:
#
#   RPMI / LPS / BG
#   rep1 / rep2
#
# Marks:
#   H3K27ac
#   H3K4me1
#   H3K4me3
#
# No files are downloaded yet.
# ============================================================

rm(list = ls())

if (!requireNamespace("GEOquery", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages(
            "BiocManager",
            repos = "https://cloud.r-project.org"
        )
    }

    BiocManager::install(
        "GEOquery",
        ask = FALSE,
        update = FALSE
    )
}

library(GEOquery)

cat("\n========================================\n")
cat("12a — DAY-6 CHROMATIN FILE AUDIT\n")
cat("========================================\n\n")

samples <- data.frame(

    GSM = c(

        # H3K27ac
        "GSM2263027", # RPMI d6 rep1
        "GSM2262963", # RPMI d6 rep2
        "GSM2262956", # LPS d6 rep1
        "GSM2263053", # LPS d6 rep2
        "GSM2262932", # BG d6 rep1
        "GSM2263034", # BG d6 rep2

        # H3K4me1
        "GSM2263033", # RPMI d6 rep1
        "GSM2263007", # RPMI d6 rep2
        "GSM2262996", # LPS d6 rep1
        "GSM2262953", # LPS d6 rep2
        "GSM2263030", # BG d6 rep1
        "GSM2262961", # BG d6 rep2

        # H3K4me3
        "GSM2262940", # RPMI d6 rep1
        "GSM2263015", # RPMI d6 rep2
        "GSM2263031", # LPS d6 rep1
        "GSM2263046", # LPS d6 rep2
        "GSM2262980", # BG d6 rep1
        "GSM2263054"  # BG d6 rep2
    ),

    Condition = rep(
        c(
            "RPMI_rep1",
            "RPMI_rep2",
            "LPS_rep1",
            "LPS_rep2",
            "BG_rep1",
            "BG_rep2"
        ),
        3
    ),

    Mark = rep(
        c(
            "H3K27ac",
            "H3K4me1",
            "H3K4me3"
        ),
        each = 6
    ),

    stringsAsFactors = FALSE
)

audit <- list()

for (i in seq_len(nrow(samples))) {

    gsm <- samples$GSM[i]

    cat(
        sprintf(
            "[%02d/%02d] %s  %s  %s\n",
            i,
            nrow(samples),
            gsm,
            samples$Mark[i],
            samples$Condition[i]
        )
    )

    x <- GEOquery::getGEOSuppFiles(
        gsm,
        fetch_files = FALSE
    )

    if (nrow(x) == 0) {

        audit[[i]] <- data.frame(
            GSM = gsm,
            Mark = samples$Mark[i],
            Condition = samples$Condition[i],
            Filename = NA_character_,
            URL = NA_character_,
            stringsAsFactors = FALSE
        )

    } else {

        audit[[i]] <- data.frame(
            GSM = gsm,
            Mark = samples$Mark[i],
            Condition = samples$Condition[i],
            Filename = rownames(x),
            URL = x$url,
            stringsAsFactors = FALSE
        )
    }
}

audit <- do.call(
    rbind,
    audit
)

dir.create(
    "results/12_kinase_reevaluation",
    recursive = TRUE,
    showWarnings = FALSE
)

write.csv(
    audit,
    "results/12_kinase_reevaluation/12a_day6_chromatin_file_manifest.csv",
    row.names = FALSE
)

cat("\n========================================\n")
cat("AUDIT RESULT\n")
cat("========================================\n\n")

print(
    audit,
    row.names = FALSE
)

cat("\nFiles identified:", nrow(audit), "\n")

cat(
    "Missing URLs:",
    sum(is.na(audit$URL)),
    "\n"
)

cat("\n========================================\n")
cat("12a COMPLETE\n")
cat("========================================\n")
