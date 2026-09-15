# ============================================================
# 12n — Update analysis_report.Rmd after corrected chromatin audit
# ============================================================

rm(list = ls())

file <- "analysis_report.Rmd"

if (!file.exists(file)) {
  stop("analysis_report.Rmd not found")
}

x <- readLines(file, warn = FALSE, encoding = "UTF-8")

# Backup only if it does not already exist
backup <- "analysis_report_before_Tobias_correction.Rmd"

if (!file.exists(backup)) {
  file.copy(file, backup)
}

txt <- paste(x, collapse = "\n")

replace_section <- function(text, start_heading, end_heading, replacement) {

  start <- regexpr(
    paste0("(?m)^", start_heading, "$"),
    text,
    perl = TRUE
  )

  end <- regexpr(
    paste0("(?m)^", end_heading, "$"),
    text,
    perl = TRUE
  )

  if (start[1] == -1 || end[1] == -1 || end[1] <= start[1]) {
    stop(
      paste(
        "Could not locate section:",
        start_heading,
        "->",
        end_heading
      )
    )
  }

  paste0(
    substr(text, 1, start[1] - 1),
    replacement,
    "\n\n",
    substr(text, end[1], nchar(text))
  )
}

# ------------------------------------------------------------
# 1. Replace old Step-08 final-table loader
# ------------------------------------------------------------

old_loader <- paste(
  'top <- safe_read(',
  '  "results/08_final_chromatin/08_TOP_integrated_candidates.csv"',
  ')',
  sep = "\n"
)

new_loader <- paste(
  '# Historical Step-08 table retained for provenance',
  'legacy_top <- safe_read(',
  '  "results/08_final_chromatin/08_TOP_integrated_candidates.csv"',
  ')',
  '',
  '# Corrected formal chromatin audit outputs',
  'final_status <- safe_read(',
  '  "results/12_kinase_reevaluation/12m_FINAL_candidate_status.csv"',
  ')',
  '',
  'local_promoter <- safe_read(',
  '  "results/12_kinase_reevaluation/12m_local_formal_promoter_results.csv"',
  ')',
  '',
  'tobias_benchmark <- safe_read(',
  '  "results/12_kinase_reevaluation/12m_Tobias_independent_benchmark.csv"',
  ')',
  '',
  'audit_repro <- safe_read(',
  '  "results/12_kinase_reevaluation/12m_reproducibility_summary.csv"',
  ')',
  '',
  'priority_candidates <- final_status$SYMBOL[',
  '  final_status$Downstream_status == "PRIORITIZE"',
  ']',
  sep = "\n"
)

if (!grepl(old_loader, txt, fixed = TRUE)) {
  stop("Could not find old Step-08 loader block")
}

txt <- sub(
  old_loader,
  new_loader,
  txt,
  fixed = TRUE
)

# ------------------------------------------------------------
# 2. Executive Summary
# ------------------------------------------------------------

executive <- paste(
  "# Executive Summary",
  "",
  "This report presents a downstream computational re-analysis of the human monocyte innate-memory dataset **GSE85246/GSE85245**, integrating RNA-seq state, β-glucan rescue/restimulation behavior, and chromatin information.",
  "",
  "The RNA analysis identified `r n_memory` Day-6 kinase-family candidates and assigned integrated trajectory priorities using acute response, persistent Day-6 state, β-glucan rescue, and restimulation information.",
  "",
  "The original Step-07/07b/08 chromatin workflow used directional H3K27ac/H3K4me1 evidence and produced MAP3K8, BMPR1A, and JAK3 as the strongest exploratory integrated candidates. A subsequent formal promoter audit showed that this earlier ranking should **not** be treated as the final promoter-level statistical result.",
  "",
  "The corrected audit uses TSS ±2 kb promoter windows, mean BigWig signal, log transformation, an unpaired limma design, and kinase-restricted statistical testing. Three RPMI Day-6 replicate-2 tracks labelled `NotNormalized` were explicitly audited and rescaled locally.",
  "",
  "The local reconstruction did **not** reproduce every numerical FDR from the independent corrected analysis. Therefore, local statistics and independently reported corrected promoter results are kept separate.",
  "",
  "The current downstream priority produced by the integration workflow is:",
  "",
  "```{r corrected-priority-summary, results='asis'}",
  'cat(paste0("**", paste(priority_candidates, collapse = ", "), "**"))',
  "```",
  "",
  "**JAK3** is supported by the project's RNA analysis and by independently corrected H3K27ac/H3K4me3 promoter evidence. **EPHB2** is a washout-emergent RNA candidate with independently corrected H3K27ac promoter support.",
  "",
  "MET, MAP3K8, and BMPR1A remain RNA-supported candidates but are not currently classified as formally promoter-chromatin-supported.",
  "",
  "> The exact corrected Tobias FDR values are treated as independent benchmark results and are not claimed as locally reproduced statistics.",
  sep = "\n"
)

txt <- replace_section(
  txt,
  "# Executive Summary",
  "# 1\\. Biological Question",
  executive
)

# ------------------------------------------------------------
# 3. Analytical workflow
# ------------------------------------------------------------

workflow <- paste(
  "# 3. Analytical Workflow",
  "",
  "The project now contains two linked workflows.",
  "",
  "The original RNA and exploratory chromatin workflow is:",
  "",
  "```text",
  "01_setup_and_data.R",
  "        ↓",
  "02_initial_LPS_response.R",
  "        ↓",
  "03_persistent_kinase_memory.R",
  "        ↓",
  "04_BG_rescue_and_restimulation.R",
  "        ↓",
  "05_integrated_kinase_trajectory.R",
  "        ↓",
  "06_make_final_kinase_figure.R",
  "        ↓",
  "07_epigenetic_persistence.R",
  "        ↓",
  "07b_unbiased_distal_chromatin.R",
  "        ↓",
  "08_make_final_chromatin_figures.R",
  "```",
  "",
  "This historical workflow can be rerun with:",
  "",
  "```bash",
  "Rscript --vanilla run_all.R",
  "```",
  "",
  "The corrected Day-6 promoter audit is run separately because it requires local GEO BigWig files that are intentionally excluded from Git:",
  "",
  "```text",
  "12a file audit",
  "    ↓",
  "12c BigWig global QC",
  "    ↓",
  "12d TSS ±2 kb promoter quantification",
  "    ↓",
  "12e RPMI_d6_rep2 normalization correction",
  "    ↓",
  "12j local kinase-restricted unpaired limma",
  "    ↓",
  "12m transparent final integration",
  "```",
  "",
  "Run the corrected audit with:",
  "",
  "```bash",
  "Rscript --vanilla run_corrected_chromatin_audit.R",
  "```",
  "",
  "The final integration explicitly distinguishes locally reproduced chromatin statistics from the independently corrected promoter results.",
  sep = "\n"
)

txt <- replace_section(
  txt,
  "# 3\\. Analytical Workflow",
  "# 4\\. Day-1 Initial LPS Response",
  workflow
)

# ------------------------------------------------------------
# 4. Dataset chromatin description
# ------------------------------------------------------------

txt <- gsub(
  "- \\*\\*Chromatin:\\*\\* H3K27ac and H3K4me1 at predefined promoters and condition-blind selected non-promoter loci\\.",
  paste(
    "- **Chromatin:** H3K27ac, H3K4me1, and H3K4me3.",
    "The earlier exploratory workflow evaluated H3K27ac/H3K4me1 promoter and non-promoter patterns;",
    "the corrected formal promoter audit prioritizes H3K27ac and H3K4me3 at TSS ±2 kb.",
    sep = " "
  ),
  txt
)

# ------------------------------------------------------------
# 5. Replace Sections 8-15
# ------------------------------------------------------------

corrected_chromatin <- paste(
  "# 8. Legacy Exploratory Chromatin Analysis",
  "",
  "Steps 07, 07b, and 08 are retained as analysis provenance.",
  "",
  "The earlier chromatin workflow evaluated directional H3K27ac and H3K4me1 patterns at promoter and condition-blind selected non-promoter regions. Step 07b corrected the original post-hoc distal-window selection problem by selecting one non-promoter region per gene without using condition-specific direction.",
  "",
  "That exploratory workflow produced MAP3K8, BMPR1A, and JAK3 as the highest integrated candidates. These outputs remain useful for historical and locus-level exploration, but they no longer define the final formal promoter-supported candidate hierarchy.",
  "",
  "```{r legacy-chromatin-figure}",
  'safe_include("figures/08_final_chromatin/08A_corrected_chromatin_support_all24.png")',
  "```",
  "",
  "# 9. Corrected Day-6 Promoter Audit",
  "",
  "A subsequent audit was performed after identifying three RPMI Day-6 replicate-2 BigWig tracks labelled `NotNormalized`:",
  "",
  "- GSM2262963 — H3K27ac",
  "- GSM2263007 — H3K4me1",
  "- GSM2263015 — H3K4me3",
  "",
  "The corrected local audit uses:",
  "",
  "- TSS ±2 kb promoter windows;",
  "- mean BigWig promoter signal;",
  "- log-transformed signal;",
  "- an unpaired limma model;",
  "- statistical testing restricted to kinase genes.",
  "",
  "The locally implemented transformation is `log2(signal + 1)`. The exact pseudocount used in the independent analysis was not specified.",
  "",
  "The three unnormalized RPMI Day-6 replicate-2 tracks were rescaled locally using normalized-peer global-mean matching. This is an approximation of the correction described during independent review; the exact independent scaling factors were not supplied.",
  "",
  "```{r audit-repro-table}",
  "knitr::kable(",
  "  audit_repro,",
  '  caption = "Corrected chromatin audit reproducibility status"',
  ")",
  "```",
  "",
  "# 10. Local Formal Promoter Results",
  "",
  "The local kinase-restricted unpaired limma reconstruction is reported separately from the independent benchmark.",
  "",
  "```{r local-promoter-table}",
  "local_show <- local_promoter[",
  '  local_promoter$SYMBOL %in% c("JAK3", "EPHB2", "MET", "MAP3K8", "BMPR1A"),',
  "  ,",
  "  drop = FALSE",
  "]",
  "",
  "knitr::kable(",
  "  local_show,",
  '  caption = "Local formal Day-6 promoter reproduction"',
  ")",
  "```",
  "",
  "The local model shows strong nominal promoter signal for JAK3 and MET, but the kinase-family adjusted FDRs do not reproduce the independent corrected values. Most importantly, the independently reported EPHB2 H3K27ac result is not reproduced locally.",
  "",
  "Therefore the local analysis is retained as a transparent reproduction audit rather than presented as an exact replication.",
  "",
  "# 11. Independent Corrected Promoter Benchmark",
  "",
  "The following corrected values were reported independently during methodological review and are stored separately from local statistics:",
  "",
  "```{r tobias-benchmark-table}",
  "knitr::kable(",
  "  tobias_benchmark,",
  '  caption = "Independent corrected promoter benchmark"',
  ")",
  "```",
  "",
  "The independently corrected interpretation is:",
  "",
  "- **JAK3:** H3K27ac FDR 0.0065 and H3K4me3 FDR 0.0141;",
  "- **EPHB2:** significant H3K27ac support, reported FDR 0.030;",
  "- **MET:** H3K27ac changed from FDR 0.032 before correction to 0.051 after correction;",
  "- **MAP3K8:** no significant corrected formal promoter support;",
  "- **BMPR1A:** no significant corrected formal promoter support.",
  "",
  "For EPHB2 and MET, the quoted correspondence did not explicitly identify the relevant Day-6 contrast in the same sentence as the reported FDR. The benchmark table therefore avoids assigning those values to multiple contrasts.",
  "",
  "# 12. Current RNA + Chromatin Integration",
  "",
  "The final Step-12m integration combines the project's RNA evidence with the independently corrected promoter review while explicitly recording whether the chromatin result was locally reproduced.",
  "",
  "```{r final-status-table}",
  "knitr::kable(",
  "  final_status,",
  '  caption = "Current integrated candidate status"',
  ")",
  "```",
  "",
  "The generated downstream priority is:",
  "",
  "```{r final-priority, results='asis'}",
  'cat(paste0("**", paste(priority_candidates, collapse = ", "), "**"))',
  "```",
  "",
  "This means that **JAK3 and EPHB2 are prioritized for downstream external validation**.",
  "",
  "It does **not** mean that the local chromatin model independently discovered both genes.",
  "",
  "# 13. Candidate Interpretation",
  "",
  "## JAK3",
  "",
  "JAK3 remains the strongest cross-layer candidate. It has persistent RNA support in the project and independently corrected promoter support for both H3K27ac and H3K4me3.",
  "",
  "The local reconstruction also shows strong positive nominal promoter effects, although the exact independent FDR values were not reproduced.",
  "",
  "## EPHB2",
  "",
  "EPHB2 is not part of the original 24 Day-6 kinase-FDR set. It emerged from the broader washout-associated RNA review and gained independently corrected H3K27ac promoter support.",
  "",
  "The local promoter reconstruction does not reproduce the reported significant EPHB2 result. This discrepancy is retained explicitly as a reproducibility limitation.",
  "",
  "## MET",
  "",
  "MET remains RNA-supported. Its initially significant H3K27ac result became borderline/non-significant after correction of the RPMI Day-6 normalization issue (reported FDR 0.051), so it is not currently classified as a confirmed RNA + promoter-chromatin hit.",
  "",
  "## MAP3K8 and BMPR1A",
  "",
  "MAP3K8 and BMPR1A retain RNA-level biological interest but do not have significant corrected formal promoter support in the independent review. Their earlier H3K27ac/H3K4me1 directional patterns remain exploratory provenance.",
  "",
  "# 14. Corrected Final Integration Figure",
  "",
  "```{r corrected-final-figure}",
  'safe_include("figures/12_kinase_reevaluation/12m_FINAL_candidate_summary.png")',
  "```",
  "",
  "This figure represents the **downstream integration decision**, not an assertion that all chromatin statistics were reproduced locally.",
  "",
  "# 15. IGV and Locus-Level Interpretation",
  "",
  "The previous IGV inspections of MAP3K8, BMPR1A, and JAK3 remain useful as historical qualitative checks, but IGV visualization is not a statistical test and does not override the corrected promoter analysis.",
  "",
  "Future locus-level visualization should prioritize **JAK3 and EPHB2**, particularly the Day-6 H3K27ac and H3K4me3 promoter tracks in hg19.",
  "",
  "The earlier non-promoter regions remain computationally defined regions rather than experimentally validated enhancers.",
  sep = "\n"
)

txt <- replace_section(
  txt,
  "# 8\\. Chromatin Persistence Analysis",
  "# 16\\. Limitations",
  corrected_chromatin
)

# ------------------------------------------------------------
# Write updated report
# ------------------------------------------------------------

writeLines(
  strsplit(txt, "\n", fixed = TRUE)[[1]],
  file,
  useBytes = TRUE
)

cat("\nanalysis_report.Rmd updated successfully.\n")
