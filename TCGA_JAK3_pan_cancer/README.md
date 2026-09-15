# TCGA JAK3 Pan-Cancer Analysis

## Objective

Systematically evaluate whether the innate-memory kinase candidate JAK3 shows robust dysregulation across TCGA cancer types.

## Analyses

For each TCGA cancer type:

- Tumor versus normal expression
- Matched tumor-normal analysis
- Early-stage tumor versus normal
- Stage-wise expression trends
- Grade-wise expression trends where available
- Aberrant/high expression prevalence
- Survival analysis
- Immune infiltration controls

## Immune correction

Because JAK3 is strongly associated with immune biology, bulk RNA-seq signals are evaluated using:

- PTPRC/CD45
- CD3D/CD3E
- MS4A1/CD79A
- NKG7
- LST1/TYROBP
- FCER1G
- CD68

A composite immune score is generated from immune marker expression.

## Current validation cancer

KIRC

Initial findings:

- JAK3 is strongly elevated in tumor compared with normal kidney.
- Stage I tumors show strong elevation.
- JAK3 correlates with immune infiltration.
- Tumor and Stage I associations remain significant after immune adjustment.

