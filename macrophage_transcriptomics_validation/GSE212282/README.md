# GSE212282

JAK3 Day6 post-rest validation across multiple trained-immunity conditions in primary human monocyte-derived macrophages.

## Design

- Homo sapiens
- Three donors: A, B, C
- Initial exposure conditions:
  - RPMI
  - beta-glucan
  - DMI
  - itaconate
  - beta-glucan + DMI
- Day6 macrophage transcriptome after five days of differentiation/rest
- No secondary LPS restimulation included in the primary analysis

## JAK3 results

### Beta-glucan vs RPMI
- logFC = 0.464
- 95% CI = -0.245 to 1.172
- FDR = 0.527
- 2/3 donors positive

### DMI vs RPMI
- logFC = 0.036
- 95% CI = -0.673 to 0.745
- FDR = 0.997
- 2/3 donors positive

### Itaconate vs RPMI
- logFC = 0.550
- 95% CI = -0.158 to 1.259
- FDR = 0.805
- 3/3 donors positive

### Beta-glucan + DMI vs RPMI
- logFC = 0.113
- 95% CI = -0.596 to 0.822
- FDR = 0.915
- 2/3 donors positive

## Classification

No statistically supported persistent JAK3 deregulation was detected in any Day6 comparison.

Itaconate showed a consistent positive JAK3 tendency across all three donors, but the effect was not statistically supported after genome-wide multiple-testing correction.
