# TCGA-BRCA — JAK3 and EPHB2

## Overview

This analysis evaluates **JAK3** and **EPHB2** expression and clinical associations in TCGA Breast Invasive Carcinoma (TCGA-BRCA).

The analysis includes:

- tumor versus normal expression
- matched tumor-normal comparison
- normal-derived aberrant-expression prevalence
- pathologic stage association
- tumor-grade availability
- overall survival
- proportional-hazards diagnostics

No requirement for statistical significance was imposed.

---

## Cohort

RNA-seq STAR-count files:

- Primary Tumor: **1,111**
- Solid Tissue Normal: **113**
- Total eligible files: **1,224**

After collapsing multiple tumor aliquots to the patient level:

- Tumor patients: **1,095**
- Normal patients: **113**

Matched tumor-normal patients:

- **113**

Genes:

- JAK3 — ENSG00000105639
- EPHB2 — ENSG00000133216

Expression is reported as TPM.

---

## Tumor versus normal

### JAK3

Tumor expression was only modestly higher than normal:

- tumor median TPM: **5.31**
- normal median TPM: **4.78**
- tumor-normal FDR: **0.0551**

The difference therefore did **not** satisfy the FDR < 0.05 criterion.

Matched analysis was also nonsignificant:

- matched patients: **113**
- tumor higher: **56**
- tumor lower: **57**
- paired FDR: **0.353**

Approximately **20.2%** of tumors exceeded the normal 95th-percentile threshold, but because the overall tumor-normal comparison was nonsignificant, this prevalence is treated as **descriptive only**.

### EPHB2

EPHB2 was clearly elevated in BRCA tumors:

- tumor median TPM: **3.18**
- normal median TPM: **1.70**
- tumor-normal comparison: highly significant

Matched analysis strongly supported the same direction:

- matched patients: **113**
- tumor higher: **95**
- tumor lower: **18**
- paired FDR: **9.80 × 10^-13**

Using the normal 95th percentile as the high-expression threshold:

- **441 / 1,095 tumors**
- **40.3%**

were classified as EPHB2-high.

---

## Pathologic stage

Formal stage testing used Stage I-IV after collapsing substages.

### JAK3

No stage association was detected:

- Kruskal-Wallis FDR: **0.936**
- Spearman stage-trend FDR: **0.934**

### EPHB2

EPHB2 showed increasing median expression and high-expression prevalence across later stages, but the formal tests did not remain significant after multiple-testing correction:

- Kruskal-Wallis FDR: **0.374**
- stage-trend FDR: **0.0695**

Therefore the stage pattern is considered **descriptive/suggestive rather than statistically significant**.

Stage 0 contained only four tumors and was retained for descriptive purposes only.

---

## Tumor grade

Tumor grade could not be evaluated reliably.

Among 1,098 TCGA-BRCA clinical cases:

- usable tumor-grade records: **0**

The GDC `tumor_grade` field was essentially unpopulated for this cohort.

Grade was therefore recorded as **not assessable**, rather than inferred from another clinical variable.

---

## Overall survival

Corrected GDC survival extraction produced:

- usable OS patients: **1,095**
- deaths: **151**
- censored: **944**

The primary survival predictor was continuous:

`log2(TPM + 1)`

### Model specification

An initial model used gene expression + age + early/advanced pathologic stage.

Proportional-hazards testing showed that stage violated the PH assumption:

- JAK3 model stage PH P = **0.0021**
- EPHB2 model stage PH P = **0.0015**

The final primary Cox model therefore used:

`Survival ~ gene expression + age + strata(Early/Advanced stage)`

This allows different baseline hazards for early and advanced stage.

Sex was not included in the primary adjusted model because there was only **one male death** in the complete-case cohort, making the sex effect poorly estimable.

Grade was unavailable.

### JAK3

Final stage-stratified continuous model:

- HR: **0.912**
- FDR: **0.261**
- gene PH P: **0.768**
- global PH P: **0.725**

No significant survival association was detected.

Normal-derived high-expression group sensitivity analysis:

- HR: **0.874**
- FDR: **0.546**

Also nonsignificant.

### EPHB2

Final stage-stratified continuous model:

- HR: **1.209**
- FDR: **0.154**
- gene PH P: **0.967**
- global PH P: **0.695**

The continuous primary survival association was therefore nonsignificant.

The normal-derived EPHB2-high group showed an adverse direction:

- adjusted stage-stratified HR: **1.405**
- FDR: **0.101**

but this also did not remain significant after FDR correction.

Descriptively, the earlier unadjusted groups had median OS of approximately:

- reference: **12.20 years**
- EPHB2-high: **9.48 years**

These descriptive values are not treated as evidence of an independently significant prognostic association.

---

## Overall interpretation

**EPHB2 is reproducibly overexpressed in TCGA-BRCA.** This finding is supported by both unmatched and matched tumor-normal analyses, and approximately 40% of tumors exceed the normal-derived high-expression threshold.

However, EPHB2 does not show a statistically significant association with pathologic stage or overall survival after appropriate adjustment, proportional-hazards correction, and multiple-testing correction.

**JAK3 shows no convincing BRCA-wide tumor-normal dysregulation, stage association, or survival association.**

No causal interpretation is implied.

---

## Figure

`figures/01_tumor_normal_JAK3_EPHB2.png`

shows JAK3 and EPHB2 tumor-normal expression distributions using TCGA-BRCA STAR-count RNA-seq TPM values.
