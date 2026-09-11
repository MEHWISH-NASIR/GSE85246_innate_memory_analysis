# TCGA-BLCA analysis of JAK3 and EPHB2

## Objective

This analysis evaluates the two priority kinase candidates, **JAK3** and **EPHB2**, in TCGA Bladder Urothelial Carcinoma (BLCA).

The analysis includes:

- tumor versus normal expression
- matched tumor-normal comparison
- prevalence of aberrant expression
- pathological stage association
- tumor grade association
- overall survival analysis

## Dataset

RNA-seq STAR-count data and clinical information were obtained from the NCI Genomic Data Commons.

Expression extraction included:

- 412 primary-tumor samples
- 19 solid-tissue normal samples

After collapsing duplicate tumor aliquots to the patient level:

- 406 primary-tumor patients
- 19 normal samples
- 19 matched tumor-normal patients

Expression was evaluated using TPM and log2(TPM + 1) where appropriate.

---

# JAK3

## Tumor versus normal

JAK3 did not show significant overall tumor-associated dysregulation in BLCA.

- Tumor median TPM: 6.31
- Normal median TPM: 5.80
- tumor-normal FDR: 0.392

Matched tumor-normal analysis was also not significant:

- 10/19 matched tumors had higher expression
- 9/19 had lower expression
- paired FDR: 0.533

Using the normal 95th percentile as a descriptive high-expression threshold:

- 21/406 tumors exceeded the threshold
- prevalence: 5.2%

Because the overall tumor-normal comparison was not significant, this prevalence estimate is treated as descriptive rather than evidence of BLCA-wide JAK3 overexpression.

## Pathological stage

JAK3 showed no significant association with pathological stage.

- Kruskal-Wallis FDR: 0.695
- stage-trend FDR: 0.645

Therefore, the data do not support a progressive stage-related change in JAK3 expression in BLCA.

## Tumor grade

JAK3 expression was significantly higher in High Grade than Low Grade tumors.

- Low Grade median TPM: 3.30
- High Grade median TPM: 6.18
- FDR: 3.25e-04

### Interpretation

JAK3 is not generally overexpressed in BLCA relative to normal bladder tissue, but within tumors its expression is significantly higher in High Grade disease.

## Overall survival

Continuous JAK3 expression was not associated with overall survival.

Age-, sex-, and stage-adjusted Cox regression:

- HR: 0.91
- 95% CI: 0.76-1.10
- FDR: 0.354

The proportional-hazards assumption was satisfied.

The normal-derived high-JAK3 group was also not significantly associated with survival:

- HR: 0.78
- FDR: 0.483
- log-rank FDR: 0.482

### JAK3 conclusion

JAK3 shows no convincing BLCA-wide tumor-associated upregulation or survival association, although higher expression is associated with High Grade disease.

---

# EPHB2

## Tumor versus normal

EPHB2 was significantly elevated in BLCA tumors.

- Tumor median TPM: 7.00
- Normal median TPM: 2.71
- tumor-normal FDR: 1.51e-04

Matched tumor-normal analysis strongly supported the same direction:

- 17/19 matched tumors had higher EPHB2 expression
- paired FDR: 0.0028

Using the normal 95th percentile as the high-expression threshold:

- 198/406 tumors exceeded the threshold
- prevalence: 48.8%

Because BLCA contains only 19 normal samples, the percentile-based prevalence estimate should be interpreted with appropriate caution.

## Pathological stage

EPHB2 showed a significant positive monotonic stage trend.

Median TPM:

- Stage I: 3.96
- Stage II: 5.75
- Stage III: 6.80
- Stage IV: 8.61

Stage-trend statistics:

- Spearman rho: 0.167
- FDR: 0.0126

The overall Kruskal-Wallis comparison across stages was not significant after FDR correction:

- FDR: 0.116

Therefore, the evidence supports a positive monotonic stage trend rather than a general difference among all stage groups.

## Tumor grade

EPHB2 was significantly higher in High Grade tumors.

- Low Grade median TPM: 2.26
- High Grade median TPM: 6.82
- FDR: 3.25e-04

The prevalence of high EPHB2 was:

- Low Grade: 15.8%
- High Grade: 46.6%

### Interpretation

EPHB2 shows strong BLCA-associated upregulation, supported by both unpaired and matched tumor-normal analyses. Expression also increases with pathological stage and is substantially higher in High Grade tumors.

## Overall survival

Continuous EPHB2 expression was not significant after multiple-testing correction in the adjusted Cox model.

Age-, sex-, and stage-adjusted Cox regression:

- HR: 1.19
- 95% CI: 1.00-1.42
- raw P: 0.0495
- FDR: 0.099

The proportional-hazards assumption was satisfied.

However, the predefined normal-derived high-EPHB2 group showed significantly poorer overall survival:

- HR: 1.44
- FDR: 0.0345
- log-rank FDR: 0.0333
- median OS, reference group: 4.57 years
- median OS, high-EPHB2 group: 2.18 years

This grouped survival result is treated as secondary evidence and is distinguished from the nonsignificant continuous adjusted Cox result.

---

# Survival-model limitation

Tumor grade was not included in the final BLCA Cox adjustment because the complete-case Low Grade subgroup had zero deaths:

- High Grade: 125 alive, 92 deaths
- Low Grade: 18 alive, 0 deaths

This produced complete separation and an unstable grade coefficient in the standard Cox model.

The final stable multivariable survival model therefore adjusted for:

- age
- sex
- pathological stage

Grade remained fully evaluated in the separate grade-expression analysis.

---

# Overall BLCA conclusion

The BLCA analysis distinguishes the two kinase candidates.

**JAK3**

- no significant tumor-normal difference
- no significant matched tumor-normal difference
- no significant stage trend
- significantly higher expression in High Grade tumors
- no significant overall-survival association

**EPHB2**

- significant tumor-associated upregulation
- 17/19 matched tumors showed higher expression
- approximately 49% of tumors exceeded the normal 95th-percentile threshold
- significant positive stage trend
- significantly higher expression in High Grade disease
- continuous adjusted survival association did not remain significant after FDR correction
- the predefined high-expression group showed significantly poorer survival

These findings describe cancer-context associations and do not establish a causal role for either kinase in BLCA.
