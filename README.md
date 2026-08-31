# GSE85246 Innate Immune Memory Analysis

Integrated RNA-seq and ChIP-seq re-analysis of **GSE85246**, focused on persistent kinase-associated transcriptional states after LPS exposure, β-glucan rescue/restimulation, and accompanying H3K27ac/H3K4me1 chromatin patterns.

## Live Analysis Report

[Open the complete interactive HTML report](https://mehwish-nasir.github.io/GSE85246_innate_memory_analysis/)


## Main result

The analysis produced the following evidence hierarchy:

```text
24 Day-6 memory-associated kinases
        ↓
10 HIGH RNA-priority kinases
        ↓
8/10 with corrected chromatin support
        ↓
3/10 with corrected dual-mark support
        ↓
MAP3K8, BMPR1A, JAK3
```

The three strongest integrated candidates are:

- **MAP3K8**
- **BMPR1A**
- **JAK3**

These candidates combine persistent Day-6 RNA-state evidence, directional rescue/restimulation support, and corrected H3K27ac/H3K4me1 evidence.

> Chromatin results are interpreted as directional/descriptive same-locus support. They do not establish causal epigenetic inheritance.

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

Primary dataset:

- **GSE85246**
- Human monocyte/macrophage innate immune memory / LPS tolerance system
- RNA-seq
- H3K27ac ChIP-seq
- H3K4me1 ChIP-seq

The project is a **downstream re-analysis** of processed expression estimates and normalized chromatin tracks. It is not intended as an exact reproduction of the original publication's differential-expression pipeline.

---

## Canonical pipeline

The current reproducible workflow is:

```text
scripts/01_setup_and_data.R
        ↓
scripts/02_initial_LPS_response.R
        ↓
scripts/03_persistent_kinase_memory.R
        ↓
scripts/04_BG_rescue_and_restimulation.R
        ↓
scripts/05_integrated_kinase_trajectory.R
        ↓
scripts/06_make_final_kinase_figure.R
        ↓
scripts/07_epigenetic_persistence.R
        ↓
scripts/07b_unbiased_distal_chromatin.R
        ↓
scripts/08_make_final_chromatin_figures.R
```

Run the complete pipeline from the project root with:

```bash
Rscript --vanilla run_all.R
```

Each analysis step is launched in a separate clean R process. The run stops immediately if any step fails.

Run status is written to:

```text
results/run_all_status.csv
```

Logs are written under:

```text
logs/run_all/
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

### 07 — Initial chromatin persistence analysis

Quantifies H3K27ac and H3K4me1 signal across promoter and non-promoter windows.

The initial exploratory distal implementation scanned many windows per gene and could therefore inflate apparent support through post-hoc window selection.

### 07b — Corrected unbiased distal chromatin analysis

Corrects the distal-window selection problem.

For each gene, one non-promoter window is selected using **condition-blind overall chromatin abundance** before testing Day-1, Day-6, or β-glucan direction.

Corrected results:

```text
Promoter any-mark support:       15 / 24
Promoter dual-mark support:       6 / 24

Fixed non-promoter any support:  14 / 24
Fixed non-promoter dual support:  3 / 24

Any corrected chromatin support: 18 / 24
Any corrected dual-mark support:  7 / 24

HIGH RNA + any chromatin:         8 / 10
HIGH RNA + dual chromatin:        3 / 10
```

### 08 — Final RNA + chromatin integration

Final integrated candidates:

```text
MAP3K8
BMPR1A
JAK3
```

Main figure:

```text
figures/08_final_chromatin/08D_FINAL_RNA_chromatin_integration.png
```

---

## Final candidates

### MAP3K8

- HIGH RNA priority
- Day-6 persistent RNA up-regulation
- strong integrated directional RNA support
- dual-mark promoter chromatin support
- additional H3K27ac support at the fixed non-promoter locus

Its strongest chromatin evidence is promoter-centered.

### BMPR1A

- HIGH RNA priority
- Day-6 persistent RNA down-regulation
- strong integrated directional RNA support
- H3K27ac and H3K4me1 support at the promoter
- H3K27ac and H3K4me1 directional support at the fixed non-promoter locus

The H3K27ac rescue pattern is particularly strong. H3K4me1 rescue is supportive but shows overshoot and should be interpreted cautiously.

### JAK3

- HIGH RNA priority
- Day-6 persistent RNA up-regulation
- strong integrated directional RNA support
- dual-mark promoter support
- H3K27ac support at the fixed proximal/non-promoter locus
- H3K4me1 at the fixed proximal locus does not satisfy the final rescue criterion

The selected non-promoter JAK3 window lies close to the promoter boundary and is therefore best described as a **proximal non-promoter regulatory window**, not a validated distal enhancer.

---

## Key figures

### Integrated RNA evidence

```text
figures/06_final_kinase/
├── 06A_integrated_kinase_evidence_matrix.png
├── 06B_Day6_memory_kinase_effects.png
├── 06C_Day1_vs_Day6_kinase_effects.png
└── 06D_FINAL_integrated_kinase_summary.png
```

### Corrected chromatin integration

```text
figures/08_final_chromatin/
├── 08A_corrected_chromatin_support_all24.png
├── 08B_HIGH_RNA_Day6_effect_with_chromatin_class.png
├── 08C_TOP_candidates_chromatin_direction.png
└── 08D_FINAL_RNA_chromatin_integration.png
```

---

## Key result tables

```text
results/03_memory_kinases/
results/04_BG_rescue_restimulation/
results/05_integrated_trajectory/
results/07b_unbiased_distal_chromatin/
results/08_final_chromatin/
```

Especially useful files:

```text
results/05_integrated_trajectory/05_integrated_kinase_trajectory.csv

results/07b_unbiased_distal_chromatin/
  07b_corrected_gene_level_chromatin_summary.csv

results/08_final_chromatin/
  08_TOP_integrated_candidates.csv
  08_FINAL_professor_table.csv
  08_TOP_candidate_IGV_loci_hg19.csv
```

---

## HTML report

The complete rendered report is:

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

---

## IGV validation

Final candidate loci were visually inspected in IGV using **hg19** and normalized H3K27ac/H3K4me1 BigWig tracks.

Coordinates used for the final candidates are stored in:

```text
results/08_final_chromatin/08_TOP_candidate_IGV_loci_hg19.csv
```

IGV inspection is used as qualitative locus-level validation and is not treated as a statistical test.

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
- β-glucan rescue RNA data are limited and interpreted directionally.
- No memory kinase reached genome-wide FDR significance for the formal history-by-restimulation interaction.
- Day-1 and Day-6 ChIP samples are not longitudinally paired donors.
- Same-locus Day-1/Day-6 chromatin patterns do not prove causal inheritance of H3K27ac or H3K4me1.
- Non-promoter windows are operational computational regions, not experimentally validated enhancers.
- Step 07b supersedes the exploratory distal classification from Step 07 for final biological interpretation.

---

## Reproducibility

To rerun the complete canonical workflow:

```bash
Rscript --vanilla run_all.R
```

A successful complete run should report:

```text
Successful steps: 9 / 9
All key final outputs are present.

Final top integrated candidates: 3
MAP3K8, BMPR1A, JAK3
```



## R environment

Exact R and package versions used by this analysis are recorded in `renv.lock`.

After cloning the repository, restore the package environment from within R:

    install.packages("renv")
    renv::restore()

