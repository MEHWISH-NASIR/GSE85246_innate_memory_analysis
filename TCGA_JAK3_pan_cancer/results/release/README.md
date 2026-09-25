# JAK3 TCGA Pan-Cancer Analysis — Final Results

This directory contains the canonical result tables used in the
rendered TCGA pan-cancer report.

## Scope

JAK3 was evaluated across all 33 TCGA cancer projects, with analyses
performed where the required data were available.

The analysis includes:

- tumor versus normal expression
- matched tumor-normal expression
- Stage-I tumor versus normal expression
- stage-wise expression patterns
- grade-wise expression patterns
- prevalence of high/aberrant expression
- survival associations as secondary evidence
- leukocyte-fraction associations
- CIBERSORT immune-cell analyses
- PTPRC/CD45 adjustment
- 10-marker immune-expression adjustment
- sample-size and data-availability limitations

## Main integrated files

- `JAK3_master_pan_cancer_summary.csv`
  Complete integrated 33-cancer table.

- `JAK3_master_key_flags_corrected.csv`
  Descriptive statistical flags for filtering and report generation.

- `JAK3_complete_33cancer_analysis_availability.csv`
  Data availability and sample-size audit across all TCGA projects.

## Interpretation note

TCGA RNA-seq is bulk-tissue expression. Associations that remain after
immune-content adjustment are described as residual or adjusted
associations and should not be interpreted as proof of malignant-cell-
intrinsic expression.
