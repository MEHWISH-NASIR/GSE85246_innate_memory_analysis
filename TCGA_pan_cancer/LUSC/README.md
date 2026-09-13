# TCGA-LUSC — JAK3 and EPHB2

## Overview

This analysis evaluates **JAK3** and **EPHB2** expression and clinical associations in TCGA Lung Squamous Cell Carcinoma (TCGA-LUSC).

The workflow includes:

- tumor versus normal expression
- matched tumor-normal analysis
- normal-derived aberrant-expression prevalence
- pathologic stage association
- tumor-grade availability
- overall survival
- proportional-hazards diagnostics

Statistical significance was not required for retention or reporting.

---

## Cohort

RNA-seq STAR-count files:

- Primary Tumor: **511**
- Solid Tissue Normal: **51**
- Total eligible files: **562**

After collapsing multiple tumor aliquots to the patient level:

- Tumor patients: **501**
- Normal patients: **51**

Matched tumor-normal patients:

- **51**

Genes:

- JAK3 — ENSG00000105639
- EPHB2 — ENSG00000133216

Expression is reported as TPM.

---

## Tumor versus normal

### JAK3

JAK3 did not show significant tumor-normal dysregulation.

- tumor median TPM: **9.73**
- normal median TPM: **9.18**
- tumor-normal FDR: **0.135**

Matched analysis agreed:

- matched patients: **51**
- tumor higher: **28**
- tumor lower: **23**
- paired FDR: **0.497**

Using the normal 95th-percentile threshold:

- **119 / 501 tumors**
- **23.8%**

were above the threshold.

Because the overall tumor-normal comparison was nonsignificant, this prevalence is treated as **descriptive only**.

### EPHB2

EPHB2 was strongly elevated in LUSC tumors.

- tumor median TPM: **6.75**
- normal median TPM: **2.18**
- tumor-normal FDR: approximately **4.44 × 10^-16**

Matched analysis strongly supported the same direction:

- matched patients: **51**
- tumor higher: **44**
- tumor lower: **7**
- paired FDR: **2.48 × 10^-8**

Using the normal 95th-percentile threshold:

- **318 / 501 tumors**
- **63.5%**

were classified as EPHB2-high.

---

## Pathologic stage

Formal stage testing used Stage I-IV after collapsing substages.

### JAK3

No significant stage association was detected:

- Kruskal-Wallis FDR: **0.636**
- stage-trend FDR: **0.889**

### EPHB2

Despite strong tumor-associated upregulation, EPHB2 showed no significant stage association:

- Kruskal-Wallis FDR: **0.636**
- stage-trend FDR: **0.889**

EPHB2-high prevalence was:

- Stage I: **60.9%**
- Stage II: **63.8%**
- Stage III: **58.7%**
- Stage IV: **80.0%**

Stage IV contained only **5 patients**, so that percentage should be interpreted cautiously.

---

## Tumor grade

Tumor grade could not be assessed.

All **504** TCGA-LUSC clinical cases lacked usable GDC `tumor_grade` information.

Grade was therefore recorded as **not assessable**, rather than inferred from another clinical field.

---

## Overall survival

Corrected GDC survival extraction produced:

- usable OS patients: **498**
- deaths: **215**
- censored: **283**

The initial adjusted model included:

`gene expression + age + sex + early/advanced pathologic stage`

The complete adjusted cohort contained:

- **475 patients**
- **210 deaths**

Both sex categories and both early/advanced stage groups contained sufficient numbers of events.

---

## Proportional-hazards diagnostics

The ordinary adjusted Cox models satisfied the proportional-hazards assumption.

### JAK3 model

- gene PH P: **0.981**
- global PH P: **0.642**

### EPHB2 model

- gene PH P: **0.076**
- global PH P: **0.226**

No stage stratification or time-dependent correction was therefore required.

---

## JAK3 survival

Final adjusted continuous model:

- HR: **1.008**
- adjusted FDR: **0.911**

No significant survival association was detected.

Normal-derived JAK3-high group:

- HR: **1.192**
- Cox FDR: **0.281**
- log-rank FDR: **0.281**

This was also nonsignificant.

Descriptive median OS:

- reference: approximately **4.75 years**
- JAK3-high: approximately **3.26 years**

These descriptive differences were not statistically supported.

---

## EPHB2 survival

Final adjusted continuous model:

- HR: **1.068**
- adjusted FDR: **0.510**

No significant continuous survival association was detected.

Normal-derived EPHB2-high group:

- HR: **1.278**
- Cox FDR: **0.182**
- log-rank FDR: **0.180**

This grouped analysis was also nonsignificant.

Descriptive median OS:

- reference: approximately **5.29 years**
- EPHB2-high: approximately **3.66 years**

The shorter median survival in the high-expression group is therefore treated as descriptive rather than evidence of an independent prognostic association.

---

## Overall interpretation

**EPHB2 is strongly and reproducibly overexpressed in TCGA-LUSC.**

This is supported by both unmatched and matched tumor-normal analyses, and approximately **63.5%** of tumors exceed the normal-derived high-expression threshold.

However, EPHB2 does not show a statistically significant association with pathologic stage or overall survival.

**JAK3 shows no convincing LUSC-wide tumor-normal dysregulation, stage association, or survival association.**

No causal interpretation is implied.

---

## Figure

`figures/01_tumor_normal_JAK3_EPHB2.png`

shows TCGA-LUSC JAK3 and EPHB2 tumor-normal expression distributions.
