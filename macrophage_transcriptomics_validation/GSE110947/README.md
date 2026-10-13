# GSE110947

JAK3 Day6 post-rest validation in primary human monocyte-derived macrophages.

## Design

- Homo sapiens
- Five healthy donors: HD34, HD37, HD48, HD49, HD51
- Initial exposure for 24 h:
  - RPMI control
  - Heme
  - Beta-glucan
- Five-day differentiation/rest period
- Day6 resting macrophage RNA analyzed before secondary LPS challenge
- Deposited MMSEQ gene-level expression files

## JAK3 results

### Heme vs RPMI
- logFC = -0.505
- 95% CI = -1.132 to 0.121
- P = 0.101
- FDR = 0.359
- 5/5 donors showed lower JAK3

Interpretation: consistent negative tendency, but not statistically supported after genome-wide multiple-testing correction.

### Beta-glucan vs RPMI
- logFC = 0.045
- 95% CI = -0.582 to 0.671
- P = 0.875
- FDR = 0.967
- 2/5 donors positive and 3/5 negative

Interpretation: no evidence of persistent JAK3 deregulation.

## Overall classification

No statistically supported persistent JAK3 deregulation was detected in GSE110947.

Heme produced a highly consistent donor-wise negative tendency that requires independent replication.
