# GSE85246 Innate Immune Memory Analysis

Integrated RNA-seq and ChIP-seq re-analysis of **GSE85246 / GSE85245**, focused on persistent kinase-associated transcriptional states following LPS exposure, β-glucan rescue/restimulation, and corrected Day-6 promoter chromatin analysis using H3K27ac, H3K4me1, and H3K4me3.

## Live Analysis Report

[Open the complete interactive HTML report](https://mehwish-nasir.github.io/GSE85246_innate_memory_analysis/)


## Main result

A formal Day-6 promoter audit revised the interpretation of the earlier exploratory chromatin analysis.

The current integration is:

```text
RNA-seq kinase analysis
        ↓
persistent / washout-associated candidates
        ↓
corrected Day-6 promoter audit
        ↓
local reproduction + independent corrected benchmark
        ↓
downstream priority
        ↓
JAK3 + EPHB2
```

### Current downstream priority

- **JAK3** — persistent RNA-supported candidate with independent corrected promoter H3K27ac and H3K4me3 support.
- **EPHB2** — washout-associated candidate carried forward because of independent corrected H3K27ac promoter support.

Independent corrected promoter benchmarks:

```text
JAK3  H3K27ac  FDR = 0.0065
JAK3  H3K4me3  FDR = 0.0141
EPHB2 H3K27ac  FDR = 0.030
MET   H3K27ac  FDR = 0.051 after correction
```

**MET, MAP3K8, and BMPR1A remain RNA-supported candidates but are not currently classified as formally promoter-chromatin supported.**

> Important: the exact independent corrected FDR values were not fully reproduced by the local limma reconstruction. JAK3 + EPHB2 is therefore an integration/downstream-priority decision, not a claim that the local chromatin model independently identified both genes.

---
## Biological question

This project asks whether transient LPS exposure leaves a persistent kinase-associated molecular state after washout, whether β-glucan can reverse or remodel that state, and whether persistent transcriptional changes are accompanied by persistent regulatory chromatin changes.

The analysis follows:

```text
acute response
    ↓
post-washout memory
    ↓
kinase prioritization
    ↓
β-glucan rescue / restimulation
    ↓
integrated RNA trajectory
    ↓
chromatin persistence
    ↓
final RNA + chromatin integration
```

---

## Dataset

Primary datasets:

- **GSE85246** — SuperSeries
- **GSE85245** — ChIP-seq subseries
- Human monocyte/macrophage LPS tolerance / innate-memory system
- RNA-seq
- H3K27ac ChIP-seq
- H3K4me1 ChIP-seq
- H3K4me3 ChIP-seq
- hg19 chromatin coordinates

The project is a downstream re-analysis of processed RNA estimates and GEO BigWig tracks.

A later audit identified three RPMI Day-6 replicate-2 tracks labelled `NotNormalized`:

```text
GSM2262963 — H3K27ac
GSM2263007 — H3K4me1
GSM2263015 — H3K4me3
```

Because RPMI Day 6 is the baseline for the Day-6 promoter comparisons, this issue was explicitly evaluated in the corrected chromatin audit.

---
## Reproducible workflows

### Historical RNA / exploratory chromatin workflow

```text
01 → 02 → 03 → 04 → 05 → 06 → 07 → 07b → 08
```

Run with:

```bash
Rscript --vanilla run_all.R
```

Steps 07/07b/08 are retained for analysis provenance but their original final chromatin ranking has been superseded.

### Corrected Day-6 promoter audit

```text
12a → 12c → 12d → 12e → 12j → 12m
```

Run with:

```bash
Rscript --vanilla run_corrected_chromatin_audit.R
```

The most recent complete corrected audit reported:

```text
Successful steps: 6 / 6
Downstream priority: JAK3, EPHB2
```

---
## Analysis stages

### 01 — Setup and data audit

Checks the RNA sample structure, matched pairs, processed MMSEQ files, KinHub/OpenKinome reference list, and required analysis inputs.

### 02 — Initial Day-1 LPS response

Characterizes the acute LPS response using matched Day-1 samples.

Because Day 1 contains only two matched pairs, this step is used as trajectory context rather than as a mandatory filter for later memory candidates.

### 03 — Persistent Day-6 kinase memory

Defines the primary persistent kinase set directly from the Day-6 LPS-history versus RPMI-history comparison.

Final formal Day-6 memory set:

```text
24 kinases
```

### 04 — β-Glucan rescue and restimulation

Evaluates:

- matched Day-6 restimulation responses;
- history-by-restimulation interaction;
- β-glucan movement toward RPMI-like baseline state;
- β-glucan movement toward the naive/RPMI-like restimulation response.

The interaction analysis did not produce genome-wide FDR-significant memory kinases, so rescue/restimulation evidence is integrated directionally rather than treated as an independent significance result.

### 05 — Integrated kinase trajectory

Combines Day-1 response, Day-6 persistent memory, β-glucan rescue, and restimulation evidence.

This produces:

```text
10 HIGH RNA-priority kinases
13 INTERMEDIATE
1 EXPLORATORY
```

Priority is an integrated evidence label, not a new significance test.

### 06 — Final RNA figures

Main output:

```text
figures/06_final_kinase/06D_FINAL_integrated_kinase_summary.png
```

### 07 / 07b / 08 — Historical exploratory chromatin analysis

Steps 07, 07b, and 08 evaluated directional H3K27ac/H3K4me1 patterns at promoter and non-promoter loci.

The original analysis prioritized:

```text
MAP3K8
BMPR1A
JAK3
```

This ranking is retained as **historical exploratory provenance** and no longer defines the final formal promoter-supported candidate set.

### 12 — Corrected formal Day-6 promoter audit

The corrected audit uses:

```text
Promoter window:       TSS ±2 kb
Signal:                mean BigWig promoter signal
Transformation:        log2(signal + 1) locally
Design:                unpaired
Statistical framework: limma + empirical Bayes
Testing scope:         KinHub kinase genes
```

The independent review confirmed log transformation, kinase-restricted testing, and an unpaired design.

The local audit does not reproduce every independent corrected FDR, particularly the EPHB2 result. These two evidence sources are therefore kept explicitly separate.

The final generated Step-12m table carries **JAK3 and EPHB2** forward for downstream validation.
## Current candidate hierarchy

### Tier 1 — downstream external validation

#### JAK3

- persistent RNA-supported candidate;
- strong positive local promoter signal;
- independent corrected H3K27ac FDR = 0.0065;
- independent corrected H3K4me3 FDR = 0.0141.

JAK3 is the strongest current cross-layer candidate.

#### EPHB2

- not part of the original 24 Day-6 kinase-FDR subset;
- treated as a broader washout-associated candidate;
- independent corrected H3K27ac FDR = 0.030;
- the significant chromatin result was not reproduced locally.

EPHB2 is therefore carried forward on the basis of the integrated RNA review plus the independent corrected promoter evidence.

### Tier 2 — RNA-supported, promoter chromatin not confirmed

#### MET

MET retains RNA-level interest. Its H3K27ac result changed from FDR 0.032 to **0.051** after correction of the RPMI Day-6 normalization issue.

#### MAP3K8

MAP3K8 retains RNA-level interest but has no significant corrected formal promoter support.

#### BMPR1A

BMPR1A retains RNA-level interest but has no significant corrected formal promoter support.

---
## Key figures

### RNA evidence

```text
figures/06_final_kinase/
├── 06A_integrated_kinase_evidence_matrix.png
├── 06B_Day6_memory_kinase_effects.png
├── 06C_Day1_vs_Day6_kinase_effects.png
└── 06D_FINAL_integrated_kinase_summary.png
```

### Current corrected integration

```text
figures/12_kinase_reevaluation/
└── 12m_FINAL_candidate_summary.png
```

This figure represents the current downstream integration decision.

### Historical exploratory chromatin figures

```text
figures/08_final_chromatin/
├── 08A_corrected_chromatin_support_all24.png
├── 08B_HIGH_RNA_Day6_effect_with_chromatin_class.png
├── 08C_TOP_candidates_chromatin_direction.png
└── 08D_FINAL_RNA_chromatin_integration.png
```

The Step-08 figures are retained for provenance and do not represent the current formal promoter-supported ranking.

---
## Key result tables

RNA results:

```text
results/03_memory_kinases/
results/05_integrated_trajectory/
```

Corrected promoter audit:

```text
results/12_kinase_reevaluation/
  12_run_status.csv
  12j_all_kinase_corrected_limma.csv
  12j_target_results.csv
  12m_local_formal_promoter_results.csv
  12m_independent_corrected_benchmark.csv
  12m_FINAL_candidate_status.csv
  12m_reproducibility_summary.csv
```

The main current candidate table is:

```text
results/12_kinase_reevaluation/12m_FINAL_candidate_status.csv
```

Historical exploratory chromatin outputs remain under:

```text
results/07b_unbiased_distal_chromatin/
results/08_final_chromatin/
```

---
## HTML report

The current rendered report is:

```text
docs/analysis_report.html
```

Source:

```text
analysis_report.Rmd
```

Render with:

```bash
Rscript --vanilla render_report.R
```

The report now separates the historical Step-08 chromatin analysis from the corrected Step-12 promoter audit.

---
## IGV validation

Previous MAP3K8, BMPR1A, and JAK3 IGV inspections are retained as qualitative historical checks.

IGV visualization is not treated as a statistical test and does not override the corrected formal promoter analysis.

Future locus-level visualization should prioritize:

```text
JAK3
EPHB2
```

with particular attention to Day-6 H3K27ac and H3K4me3 promoter signal in hg19.

---
## Repository structure

```text
GSE85246_innate_memory_analysis/
├── README.md
├── run_all.R
├── analysis_report.Rmd
├── render_report.R
├── scripts/
├── metadata/
├── figures/
├── results/
├── docs/
└── data/
```

Large raw and processed data files are intentionally excluded from version control.

---

## Large files not tracked by Git

The repository excludes:

```text
data/raw/
data/processed/
*.tar
*.gz
*.fastq
*.fastq.gz
*.bam
*.bai
*.bw
*.bigWig
```

This includes the large MMSEQ files and ChIP-seq BigWig tracks.

---

## Important interpretation notes

- Day-1 RNA analysis contains only two matched pairs.
- β-glucan rescue/restimulation is interpreted primarily as supporting trajectory evidence.
- Kinase-family FDR and genome-wide FDR are distinct quantities.
- Day-6 formal promoter analysis uses an unpaired design.
- Three RPMI Day-6 replicate-2 BigWigs are labelled `NotNormalized` in GEO.
- H3K27ac and H3K4me3 are the principal marks for current formal promoter interpretation.
- H3K4me1 remains useful primarily for enhancer/non-promoter interpretation.
- Directional chromatin support is not equivalent to formal statistical significance.
- Steps 07/07b/08 are retained as historical exploratory analyses.
- The local reconstruction does not reproduce all independent corrected FDR values.
- JAK3 + EPHB2 is a downstream integration decision, not a claim that the local limma model independently identified both genes.
- Concordant RNA and chromatin changes do not establish causal epigenetic inheritance.

---
## Reproducibility

The repository contains two related workflows.

Historical Steps 01–08:

```bash
Rscript --vanilla run_all.R
```

Corrected promoter audit:

```bash
Rscript --vanilla run_corrected_chromatin_audit.R
```

The latest corrected audit completed:

```text
Successful steps: 6 / 6
Downstream priority: JAK3, EPHB2
```

The exact independent corrected promoter FDRs are stored separately from the local limma results and are not claimed as locally reproduced statistics.

## R environment

Exact R and package versions used by this analysis are recorded in `renv.lock`.

After cloning the repository, restore the package environment from within R:

    install.packages("renv")
    renv::restore()

