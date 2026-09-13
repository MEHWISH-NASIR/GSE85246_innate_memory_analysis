# TCGA-LUAD — JAK3 and EPHB2

## Overview

This analysis evaluates **JAK3** and **EPHB2** expression and clinical associations in TCGA Lung Adenocarcinoma (TCGA-LUAD).

The workflow includes tumor-normal expression, matched tumor-normal analysis, normal-derived aberrant-expression prevalence, pathologic stage, tumor-grade availability, and overall survival.

Statistical significance was not required for retention or reporting.

---

## Cohort

RNA-seq STAR-count files:

- Primary Tumor: **540**
- Solid Tissue Normal: **59**
- Total eligible files: **599**

After collapsing multiple tumor aliquots:

- Tumor patients: **517**
- Normal patients: **59**

Matched tumor-normal patients:

- **58**

Genes:

- JAK3 — ENSG00000105639
- EPHB2 — ENSG00000133216

Expression is reported as TPM.

---

## Tumor versus normal

### JAK3

JAK3 was significantly elevated in LUAD tumors.

- tumor median TPM: approximately **12.09**
- normal median TPM: **6.39**
- tumor-normal FDR: approximately **4.44 × 10^-16**

Matched analysis supported the same direction:

- matched patients: **58**
- tumor higher: **49**
- tumor lower: **9**
- paired FDR: **5.54 × 10^-9**

Using the normal 95th-percentile threshold:

- **173 / 517 tumors**
- **33.5%**

were classified as JAK3-high.

### EPHB2

EPHB2 was also strongly elevated in LUAD tumors.

- tumor median TPM: **6.82**
- normal median TPM: **2.01**

The tumor-normal P-value was below the numerical display precision used in the R output and should not be interpreted as literally zero.

Matched analysis strongly supported tumor-associated elevation:

- matched patients: **58**
- tumor higher: **53**
- tumor lower: **5**
- paired FDR: **3.73 × 10^-10**

Using the normal 95th-percentile threshold:

- **318 / 517 tumors**
- **61.5%**

were classified as EPHB2-high.

---

## Pathologic stage

Formal testing used Stage I-IV after collapsing substages.

### JAK3

No significant stage association was detected:

- Kruskal-Wallis FDR: **0.583**
- stage-trend FDR: **0.172**

### EPHB2

EPHB2 showed higher expression and high-expression prevalence in later-stage tumors, particularly Stage III-IV, but the formal associations did not remain significant after multiple-testing correction:

- Kruskal-Wallis FDR: **0.124**
- stage-trend FDR: **0.0897**

The pattern is therefore considered descriptive/suggestive rather than statistically significant.

EPHB2-high prevalence was:

- Stage I: **56.8%**
- Stage II: **58.2%**
- Stage III: **69.8%**
- Stage IV: **68.2%**

---

## Tumor grade

Tumor grade was not assessable.

All **585** TCGA-LUAD clinical records lacked usable GDC `tumor_grade` information.

Grade was therefore recorded as unavailable rather than inferred from another clinical variable.

---

## Overall survival

Corrected GDC survival extraction produced:

- usable OS patients: **513**
- deaths: **184**
- censored: **329**

The primary survival predictor was continuous:

`log2(TPM + 1)`

### Model development

The initial adjusted model included:

`gene expression + age + sex + early/advanced stage`

The EPHB2 model had a marginal global proportional-hazards result:

- global PH P: approximately **0.049**

No individual term showed a formal PH violation, although stage was the closest clinical term and was a strong prognostic variable.

A conservative stage-stratified sensitivity model was therefore used:

`Survival ~ gene expression + age + sex + strata(Early/Advanced stage)`

The corrected models satisfied proportional-hazards diagnostics.

---

## JAK3 survival

Final stage-stratified continuous model:

- HR: **0.918**
- 95% CI: **0.768–1.098**
- FDR: **0.349**
- gene PH P: **0.973**
- global PH P: **0.162**

No significant JAK3 survival association was detected.

Normal-derived JAK3-high grouped analysis:

- HR: **0.820**
- FDR: **0.247**

also showed no significant association.

---

## EPHB2 survival

Final continuous stage-stratified model:

- HR: **1.152**
- 95% CI: **1.012–1.312**
- raw P: **0.0328**
- FDR: **0.0656**
- gene PH P: **0.756**
- global PH P: **0.126**

The continuous primary association therefore did **not** remain significant after multiple-testing correction.

The predefined normal-derived EPHB2-high group showed significant adverse survival in the age- and sex-adjusted, stage-stratified sensitivity analysis:

- HR: **1.520**
- 95% CI: **1.092–2.114**
- FDR: **0.0259**
- group PH P: **0.878**
- global PH P: **0.192**

In the descriptive unadjusted Kaplan-Meier groups:

- reference median OS: approximately **6.35 years**
- EPHB2-high median OS: approximately **3.44 years**

The grouped result is treated as **secondary evidence**, because continuous expression was the prespecified primary survival analysis.

---

## Overall interpretation

Both **JAK3 and EPHB2 are strongly and reproducibly upregulated in TCGA-LUAD**, supported by unmatched and matched tumor-normal analyses.

Approximately **33.5%** of tumors exceed the normal-derived JAK3 high-expression threshold, while approximately **61.5%** exceed the corresponding EPHB2 threshold.

Neither gene shows a statistically significant association with pathologic stage after multiple-testing correction.

JAK3 shows no evidence of an overall-survival association.

EPHB2 shows an adverse survival direction in continuous analyses, but the primary continuous result does not remain significant after FDR correction. The predefined EPHB2-high group does remain significantly associated with poorer survival after age and sex adjustment and stage stratification, providing secondary prognostic evidence.

No causal interpretation is implied.

---

## Figure

`figures/01_tumor_normal_JAK3_EPHB2.png`

shows the TCGA-LUAD tumor-normal expression distributions for JAK3 and EPHB2.
