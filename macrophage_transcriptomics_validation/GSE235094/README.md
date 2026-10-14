# GSE235094

JAK3 Day6 persistence analysis after inactivated SARS-CoV-2 training in primary human monocytes.

## Design

- Homo sapiens
- Primary human monocytes
- Initial exposure to inactivated SARS-CoV-2 for 24 h
- Washout followed by 5-day rest
- Day6 transcriptome before secondary LPS stimulation

## Primary comparison

Trained, medium-only Day6 samples:

- Sample3
- Sample7
- Sample10
- Sample14
- Sample16
- Sample26

versus untrained, medium-only Day6 samples:

- Sample6
- Sample13
- Sample20
- Sample21
- Sample22

Sample5 was excluded because no public GEO sample mapping was available.

Sample8 was excluded because its public metadata were internally inconsistent regarding LPS restimulation.

## Analysis

Fractional gene-level estimated counts were analyzed using:

- TMM library normalization
- limma-voom
- genome-wide Benjamini-Hochberg FDR correction

## JAK3 result

- logFC = -1.709
- 95% CI = -3.249 to -0.169
- P = 0.0324
- FDR = 0.1436

Mean JAK3 voom log2 CPM:

- Untrained = 4.119
- Trained = 2.691

## Interpretation

JAK3 expression was lower in trained Day6 monocytes.

The nominal P value was significant and the 95% confidence interval excluded zero, but the result did not remain significant after genome-wide multiple-testing correction.

## Classification

Nominal persistent JAK3 downregulation, but not FDR-supported.

This dataset is retained as supportive monocyte evidence rather than strict monocyte-derived macrophage Tier-1 evidence.
