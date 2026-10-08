# GSE262098

JAK3 post-washout transcriptomic validation following BCG training in primary human monocytes/macrophages.

## Design

- 3 independent biological units: exp1, exp2, exp3
- RPMI control
- BCG treatment
- 24 h treatment state
- Day 6 post-washout state

The primary analysis compared:

- Day 1 BCG vs Day 1 RPMI
- Day 6 BCG vs Day 6 RPMI
- direct Day 6 vs Day 1 BCG-effect interaction

## Analysis

Deposited unnormalized HTSeq gene counts were analyzed using DESeq2.

Model:

`~ donor + time + condition + time:condition`

JAK3 Ensembl ID:

`ENSG00000105639`

## Main results

### Day 1

- log2FC: 1.826
- FDR: 0.00513
- donor consistency: 3/3 positive

### Day 6

- log2FC: 2.365
- FDR: 0.00191
- donor consistency: 3/3 positive

### Direct trajectory test

Day 6 effect minus Day 1 effect:

- log2FC difference: 0.539
- FDR: 0.818

## Interpretation

GSE262098 provides strong independent evidence that JAK3 remains positively deregulated at Day 6 after BCG exposure and washout.

The Day 6 effect is numerically larger than the Day 1 effect, but the direct interaction test does not support a statistically significant increase in effect magnitude over time.

## Classification

**Strong positive persistent JAK3 validation after BCG washout.**

## Contents

- `scripts/` — reproducible analysis
- `results/` — DESeq2 and donor-level outputs
- `figures/` — paired and effect-size plots
- `report.qmd` — report source
- `report.html` — rendered report
