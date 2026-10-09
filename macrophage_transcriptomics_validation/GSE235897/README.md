# GSE235897

JAK3 post-rest validation in human monocyte-derived macrophages after MMA training.

## Design

- 3 biological replicates per condition
- Day6 naive macrophages
- Day6 MMA-trained macrophages
- 24 h training followed by 5-day rest

## Formal analysis

Deposited MMSEQ gene-level expression estimates were jointly normalized using `readmmseq()`.

Because GEO does not explicitly identify replicate 1/2/3 as matched donors, formal inference used an unpaired limma model.

## JAK3 result

Day6 MMA vs Day6 naive:

- effect: 0.038
- 95% CI: -1.744 to 1.820
- P = 0.959
- FDR = 0.9999

Replicate-index directions:

- 1 positive
- 2 negative

## Classification

**No evidence of persistent JAK3 deregulation after MMA training.**
