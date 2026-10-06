
library(data.table)

setDTthreads(1)

root <- "TCGA_multikinase_analysis"

indir <- file.path(
  root,
  "results/final_integration"
)

outdir <- file.path(
  root,
  "results/final_summary"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n========================================\n")
cat("STEP 15: FINAL RESULTS INTERPRETATION\n")
cat("========================================\n")

# ============================================================
# 1. Load finalized data
# ============================================================

compact <- fread(
  file.path(
    indir,
    "13_supervisor_final_multikinase_comparison.csv"
  ),
  nThread = 1
)

atac <- fread(
  file.path(
    indir,
    "13_significant_ATAC_RNA_associations.csv"
  ),
  nThread = 1
)

if (nrow(compact) != 10) {
  stop("Expected 10 candidate genes.")
}

cat("Candidate genes:", nrow(compact), "\n")
cat("Significant ATAC-RNA associations:", nrow(atac), "\n")

# ============================================================
# 2. Explicit biological follow-up groups
#
# These are interpretation groups, not formal statistical tiers.
# ============================================================

compact[, Biological_followup_group :=
  fifelse(
    gene == "JAK3",
    "Reference",

    fifelse(
      gene %in% c("EPHB2", "HCK", "LYN", "IRAK2"),
      "Primary follow-up",

      fifelse(
        gene %in% c("DDR1", "RPS6KA2"),
        "Secondary follow-up",

        fifelse(
          gene %in% c("MET", "EPHB1"),
          "Contextual / lower priority",

          fifelse(
            gene == "FJX1",
            "Separate scope check",
            "Other"
          )
        )
      )
    )
  )
]

# ============================================================
# 3. Add concise interpretation per gene
# ============================================================

compact[, Interpretation := ""]

compact[
  gene == "JAK3",
  Interpretation :=
    paste0(
      "Reference persistent kinase with balanced RNA, clinical and ",
      "ATAC-RNA evidence; retained as benchmark rather than a new target."
    )
]

compact[
  gene == "EPHB2",
  Interpretation :=
    paste0(
      "Washout-emergent candidate with the highest final multi-omic index, ",
      "strong early-stage and grade-associated RNA evidence, survival support, ",
      "and broad cancer-specific ATAC-RNA coupling."
    )
]

compact[
  gene == "HCK",
  Interpretation :=
    paste0(
      "Strict persistent candidate with a multi-omic profile closest to JAK3; ",
      "strong tumor-normal, early-stage, immune-adjusted and regulatory support."
    )
]

compact[
  gene == "LYN",
  Interpretation :=
    paste0(
      "Strict persistent candidate with particularly strong normalized ATAC-RNA ",
      "support despite having only one official linked regulatory element."
    )
]

compact[
  gene == "IRAK2",
  Interpretation :=
    paste0(
      "Strict persistent candidate with balanced RNA, immune-adjusted and ",
      "ATAC-RNA evidence across multiple cancer types."
    )
]

compact[
  gene == "DDR1",
  Interpretation :=
    paste0(
      "Strict persistent candidate with strong individual ATAC-RNA correlations ",
      "but less broad regulatory support than the leading candidates."
    )
]

compact[
  gene == "RPS6KA2",
  Interpretation :=
    paste0(
      "Strict persistent candidate with broad ATAC-RNA support but weaker ",
      "overall RNA and clinical evidence."
    )
]

compact[
  gene == "MET",
  Interpretation :=
    paste0(
      "Strict persistent RNA candidate with substantial TCGA RNA and ATAC-RNA ",
      "associations, but weaker integrated progression and immune-adjusted evidence; ",
      "original macrophage chromatin support did not survive corrected normalization."
    )
]

compact[
  gene == "EPHB1",
  Interpretation :=
    paste0(
      "Strict persistent candidate with several ATAC-RNA associations but weaker ",
      "overall integrated RNA and clinical support."
    )
]

compact[
  gene == "FJX1",
  Interpretation :=
    paste0(
      "Strong RNA-level signal but limited normalized ATAC-RNA support; ",
      "kinase-family classification should be verified before mechanistic ",
      "kinase-level interpretation."
    )
]

# ============================================================
# 4. Supervisor-ready compact table
# ============================================================

priority_table <- compact[
  ,
  .(
    Gene = gene,
    Discovery_support,
    Biological_followup_group,

    Day6_log2FC,
    Day6_FDR,

    RNA_evidence_index,
    ATAC_evidence_support,
    Final_multiomic_index,

    Official_ATAC_elements,
    ATAC_FDR_sig,
    ATAC_sig_cancers,

    Median_abs_ATAC_RNA_rho,

    JAK3_domains_final,

    Interpretation
  )
]

group_order <- c(
  "Reference",
  "Primary follow-up",
  "Secondary follow-up",
  "Contextual / lower priority",
  "Separate scope check",
  "Other"
)

priority_table[, group_order :=
  match(
    Biological_followup_group,
    group_order
  )
]

setorder(
  priority_table,
  group_order,
  -Final_multiomic_index
)

priority_table[, group_order := NULL]

fwrite(
  priority_table,
  file.path(
    outdir,
    "15_final_candidate_priority_table.csv"
  )
)

# ============================================================
# 5. Pull values used in narrative
# ============================================================

getv <- function(g, col) {
  compact[
    gene == g,
    get(col)
  ][1]
}

fmt3 <- function(x) {
  sprintf("%.3f", x)
}

jak3_index  <- getv("JAK3", "Final_multiomic_index")
ephb2_index <- getv("EPHB2", "Final_multiomic_index")
hck_index   <- getv("HCK", "Final_multiomic_index")
lyn_index   <- getv("LYN", "Final_multiomic_index")
irak2_index <- getv("IRAK2", "Final_multiomic_index")

jak3_atac   <- getv("JAK3", "ATAC_evidence_support")
ephb2_atac  <- getv("EPHB2", "ATAC_evidence_support")
hck_atac    <- getv("HCK", "ATAC_evidence_support")
lyn_atac    <- getv("LYN", "ATAC_evidence_support")
irak2_atac  <- getv("IRAK2", "ATAC_evidence_support")

# ============================================================
# 6. Final results narrative
# ============================================================

results_text <- c(

"# Final multikinase RNA + ATAC results",

"",

"## Analysis framework",

paste0(
  "The final analysis integrated six evidence domains: primary tumor RNA, ",
  "early-stage RNA, stage/grade progression, survival, immune-adjusted ",
  "progression, and patient-matched ATAC-RNA regulatory coupling. ",
  "JAK3 was retained as the completed reference benchmark rather than ",
  "re-analyzed as a new candidate."
),

"",

"## Main findings",

paste0(
  "EPHB2 showed the highest final multi-omic evidence index (",
  fmt3(ephb2_index),
  "), exceeding the JAK3 reference index (",
  fmt3(jak3_index),
  "). EPHB2 therefore represents the strongest overall integrated signal, ",
  "although its evidence profile differs from JAK3 and it is classified as ",
  "washout-emergent rather than strict persistent in the original macrophage analysis."
),

"",

paste0(
  "HCK produced one of the most JAK3-like integrated profiles, with a final ",
  "multi-omic index of ",
  fmt3(hck_index),
  " and normalized ATAC-RNA support of ",
  fmt3(hck_atac),
  " compared with ",
  fmt3(jak3_atac),
  " for JAK3."
),

"",

paste0(
  "LYN showed the strongest normalized ATAC-RNA evidence among the panel ",
  "(ATAC support ",
  fmt3(lyn_atac),
  "), despite being represented by only one official linked regulatory element. ",
  "This indicates strong and recurrent coupling between accessibility at the ",
  "LYN-linked element and LYN expression across multiple cancer types."
),

"",

paste0(
  "IRAK2 showed a balanced profile across RNA, immune-adjusted and chromatin ",
  "evidence, with a final multi-omic index of ",
  fmt3(irak2_index),
  " and ATAC-RNA support of ",
  fmt3(irak2_atac),
  "."
),

"",

paste0(
  "DDR1 and RPS6KA2 retained meaningful regulatory evidence but showed less ",
  "balanced multi-domain support. MET showed extensive TCGA RNA and ATAC-RNA ",
  "associations, but its integrated progression and immune-adjusted evidence ",
  "was weaker; importantly, corrected macrophage chromatin analysis had already ",
  "reclassified MET as RNA-only in the original innate-memory dataset."
),

"",

paste0(
  "FJX1 showed a strong RNA-level profile but very limited normalized ATAC-RNA ",
  "support. It should therefore not be prioritized for mechanistic follow-up ",
  "without first confirming its kinase-family classification."
),

"",

"## Interpretation",

paste0(
  "Considering the original innate-memory classification together with TCGA ",
  "RNA, clinical and regulatory evidence, the most compelling follow-up ",
  "candidates are EPHB2, HCK, LYN and IRAK2. These candidates are not ",
  "equivalent: EPHB2 is strongest at the integrated RNA/clinical level, ",
  "HCK most closely resembles the balanced JAK3 profile, LYN is particularly ",
  "strong at the regulatory-chromatin level, and IRAK2 provides balanced ",
  "multi-layer support."
),

"",

paste0(
  "These categories are intended as biological follow-up priorities rather ",
  "than formal probability estimates or definitive mechanistic rankings."
  )
)

writeLines(
  results_text,
  file.path(
    outdir,
    "15_final_results_summary.md"
  )
)

# ============================================================
# 7. Figure legends
# ============================================================

legends <- c(

"# Final figure legends",

"",

"## Figure 1. Six-domain multi-omic support across candidate genes",

paste0(
  "Heatmap showing normalized evidence support across six domains: primary ",
  "tumor RNA, early-stage RNA, stage/grade progression, survival, ",
  "immune-adjusted progression, and patient-matched ATAC-RNA coupling. ",
  "Values represent domain-specific support metrics used in the integrated ",
  "comparison. JAK3 serves as the completed reference benchmark."
),

"",

"## Figure 2. Integrated multi-omic evidence across candidate genes",

paste0(
  "Final multi-omic evidence index for each candidate, calculated as the mean ",
  "support across the six evidence domains. Colors indicate the original ",
  "discovery classification: JAK3 reference, strict persistent candidates, ",
  "or washout-emergent candidate. The index is a descriptive integration ",
  "metric and is not a formal probability or statistical ranking."
),

"",

"## Figure 3. Multi-domain evidence profiles relative to JAK3",

paste0(
  "Comparison of normalized evidence profiles for JAK3 and the principal ",
  "follow-up candidates HCK, EPHB2, LYN and IRAK2 across six RNA, clinical ",
  "and chromatin domains. The figure highlights differences in the type of ",
  "evidence supporting each candidate rather than implying mechanistic equivalence."
),

"",

"## Figure 4. Cancer-specific ATAC-RNA coupling of candidate regulatory elements",

paste0(
  "Heatmap of significant patient-matched ATAC-RNA associations across TCGA ",
  "cancer types. Color represents the strongest significant Spearman ",
  "correlation observed for each gene-cancer pair, while numbers indicate ",
  "the number of linked regulatory elements significant at the gene-level ",
  "BH-adjusted FDR threshold. Blank cells indicate no FDR-significant ",
  "association. Positive correlations indicate higher accessibility associated ",
  "with higher expression; negative correlations indicate inverse coupling."
),

"",

"## Figure 5. RNA versus ATAC evidence landscape",

paste0(
  "Candidate genes positioned according to integrated RNA evidence and ",
  "normalized ATAC-RNA regulatory support. Dashed lines indicate the JAK3 ",
  "reference values on each axis. Point size represents the number of evidence ",
  "domains in which the candidate equals or exceeds the JAK3 benchmark, while ",
  "color indicates the original discovery classification. This representation ",
  "illustrates distinct evidence profiles rather than a single linear ranking."
)
)

writeLines(
  legends,
  file.path(
    outdir,
    "15_final_figure_legends.md"
  )
)

# ============================================================
# 8. Short supervisor summary
# ============================================================

supervisor_summary <- c(

"# Supervisor summary",

"",

paste0(
  "The expanded TCGA analysis confirms that the persistent/washout-associated ",
  "candidate set contains several genes with reproducible cancer-level RNA and ",
  "regulatory evidence beyond JAK3."
),

"",

paste0(
  "The strongest follow-up group is EPHB2, HCK, LYN and IRAK2. ",
  "EPHB2 has the highest overall integrated evidence index; HCK most closely ",
  "approximates the balanced JAK3 profile; LYN shows particularly strong ",
  "normalized ATAC-RNA coupling; and IRAK2 shows balanced RNA, immune-adjusted ",
  "and chromatin evidence."
),

"",

paste0(
  "DDR1 and RPS6KA2 remain secondary candidates. MET is retained mainly as an ",
  "RNA-associated candidate because its original macrophage chromatin support ",
  "did not survive corrected normalization. FJX1 should be kept separate until ",
  "its kinase-family classification is verified."
)
)

writeLines(
  supervisor_summary,
  file.path(
    outdir,
    "15_supervisor_summary.md"
  )
)

# ============================================================
# 9. Console summary
# ============================================================

cat("\n===== FINAL BIOLOGICAL FOLLOW-UP TABLE =====\n")

print(
  priority_table[
    ,
    .(
      Gene,
      Discovery_support,
      Biological_followup_group,
      RNA_evidence_index,
      ATAC_evidence_support,
      Final_multiomic_index,
      ATAC_FDR_sig,
      ATAC_sig_cancers,
      JAK3_domains_final
    )
  ]
)

cat("\n===== PRIMARY FOLLOW-UP CANDIDATES =====\n")

print(
  priority_table[
    Biological_followup_group ==
    "Primary follow-up",
    .(
      Gene,
      Final_multiomic_index,
      ATAC_evidence_support,
      ATAC_FDR_sig,
      ATAC_sig_cancers
    )
  ]
)

cat("\nFiles written to:\n")
cat(outdir, "\n")

print(
  list.files(
    outdir,
    full.names = FALSE
  )
)

cat("\n15 COMPLETE\n")

