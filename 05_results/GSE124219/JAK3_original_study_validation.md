# GSE124219: JAK3 ATAC-seq validation

## Study design
- Human CD14+ monocytes
- 13 paired donors
- Before BCG vaccination (D0) versus Month 3 (D90)
- Genome build: hg19
- Source: GEO GSE124219

## JAK3-associated peaks

### Peak 1
- Coordinates: chr19:17958539-17958917
- Mean normalized signal change: -19.29%
- Donors decreased: 10/13
- Wilcoxon p: 0.039795
- Two-peak BH q: 0.079590

### Peak 2
- Coordinates: chr19:17961261-17961674
- Mean normalized signal change: -1.46%
- Donors decreased: 7/13
- Wilcoxon p: 0.892578
- Two-peak BH q: 0.892578

## Original publication validation
Source: Cirovic et al., Cell (2020)
Supplementary Table S10: mmc5.xlsx

- Core differentially accessible regions: 47
- Core regions on chr19: 4
- Peak 1 overlap: No
- Peak 2 overlap: No

## Interpretation
Peak 1 displays an exploratory decrease in
normalized accessibility after BCG vaccination.

Neither JAK3-associated peak overlaps the
original study's reported core DA regions.

No genome-wide significant JAK3 accessibility
change has been established.

## Limitations
- Normalized counts used for exploratory testing.
- BH adjustment covers only two targeted peaks.
- Original full differential-accessibility statistics
  for these peaks have not been recovered.
- Precise transcript-specific promoter annotation
  requires final verification.
- Independent BCG model; not directly matched to
  the heme RNA-seq dataset.

## Analysis outputs
- JAK3_peak_summary.csv
- JAK3_donor_paired_values.csv
- Figures in 06_figures/GSE124219/
- Script: 03_scripts/gse124219_jak3_paired.py
