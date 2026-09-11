# TCGA Pan-Cancer Analysis of JAK3 and EPHB2

This section evaluates the two priority kinase candidates, **JAK3** and **EPHB2**, across TCGA cancer types.

The purpose is not to require either gene to be significant in every cancer. Each cancer is evaluated independently using the same general framework, and significant, nonsignificant, unavailable, and technically limited results are all retained.

## General analysis framework

Where appropriate data are available, each TCGA cancer is evaluated for:

1. tumor versus normal expression
2. matched tumor-normal expression
3. prevalence of aberrant expression
4. pathological stage association
5. tumor grade association
6. overall survival association

Cancer-specific clinical variables and limitations are handled according to the available TCGA data rather than forcing identical models where they are not appropriate.

## Cancer-level analyses

| TCGA project | Cancer | Status | Main pattern |
|---|---|---|---|
| TCGA-KIRC | Kidney renal clear cell carcinoma | Completed | Strong JAK3 upregulation; EPHB2 downregulated relative to normal |
| TCGA-BLCA | Bladder urothelial carcinoma | Completed | JAK3 largely NS; EPHB2 strongly upregulated |
| Other TCGA projects | Pending | To be analysed | — |

Each cancer has its own folder containing reproducible scripts, compact result tables, interpretation, and figures where available.
