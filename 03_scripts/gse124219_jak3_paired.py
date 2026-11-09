from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

base = Path(".")
source = base / "data/atac/GSE124219/GSE124219_ATAC_DESeq-norm.csv.gz"
out = base / "05_results/GSE124219"
figdir = base / "06_figures/GSE124219"

targets = [
    "chr19_17958539_17958917",
    "chr19_17961261_17961674"
]

df = pd.read_csv(source, sep=r"\s+", index_col=0)

assert df.shape == (47575, 26), f"Unexpected shape: {df.shape}"
assert df.index.is_unique, "Duplicate peak IDs found"

metadata = []
for col in df.columns:
    m = re.fullmatch(r"S(\d+)T(d0|m3)", col)
    assert m, f"Unexpected sample name: {col}"
    metadata.append((col, int(m.group(1)), m.group(2)))

donors = sorted({x[1] for x in metadata})
assert len(donors) == 13

for d in donors:
    assert f"S{d}Td0" in df.columns
    assert f"S{d}Tm3" in df.columns

assert all(p in df.index for p in targets)
assert df.loc[targets].notna().all().all()

summary = []
paired_rows = []

for number, peak in enumerate(targets, start=1):
    before = np.array([
        df.loc[peak, f"S{d}Td0"] for d in donors
    ], dtype=float)

    after = np.array([
        df.loc[peak, f"S{d}Tm3"] for d in donors
    ], dtype=float)

    delta = after - before
    n = len(delta)

    wil = stats.wilcoxon(
        delta, alternative="two-sided",
        zero_method="wilcox", method="auto"
    )
    tt = stats.ttest_rel(after, before)

    sd = np.std(delta, ddof=1)
    se = sd / np.sqrt(n)
    ci = stats.t.interval(
        0.95, df=n-1,
        loc=np.mean(delta), scale=se
    )

    for donor, a, b, change in zip(
        donors, before, after, delta
    ):
        paired_rows.append({
            "peak": peak,
            "donor": donor,
            "day0": a,
            "month3": b,
            "paired_difference": change
        })

    summary.append({
        "peak": peak,
        "n_pairs": n,
        "day0_mean": np.mean(before),
        "month3_mean": np.mean(after),
        "percent_change_of_means":
            100 * (np.mean(after) / np.mean(before) - 1),
        "mean_paired_difference": np.mean(delta),
        "median_paired_difference": np.median(delta),
        "ci95_lower": ci[0],
        "ci95_upper": ci[1],
        "cohens_dz": np.mean(delta) / sd if sd else np.nan,
        "donors_decreased": int(np.sum(delta < 0)),
        "donors_increased": int(np.sum(delta > 0)),
        "wilcoxon_p": wil.pvalue,
        "paired_t_p": tt.pvalue
    })

    fig, ax = plt.subplots(figsize=(5.5, 5))

    for a, b in zip(before, after):
        ax.plot([0, 1], [a, b],
                marker="o", alpha=0.65,
                linewidth=1.2)

    ax.plot([0, 1],
            [np.mean(before), np.mean(after)],
            color="black", marker="D",
            linewidth=2.5, label="Mean")

    ax.set_xticks([0, 1], ["Day 0", "Month 3"])
    ax.set_ylabel("Normalized ATAC peak counts")
    ax.set_title(f"GSE124219 | JAK3-associated peak {number}")
    ax.legend()
    fig.tight_layout()

    fig.savefig(figdir / f"JAK3_peak{number}_paired.png",
                dpi=300)
    fig.savefig(figdir / f"JAK3_peak{number}_paired.pdf")
    plt.close(fig)

result = pd.DataFrame(summary)

# BH correction across ONLY the two targeted peaks.
p = result["wilcoxon_p"].to_numpy()
order = np.argsort(p)
q = np.empty(len(p))
q[order] = np.minimum.accumulate(
    (p[order] * len(p) / np.arange(1, len(p)+1))[::-1]
)[::-1].clip(max=1)

result["BH_q_two_peaks_only"] = q
result["analysis_type"] = "Exploratory; normalized counts"

result.to_csv(out / "JAK3_peak_summary.csv", index=False)
pd.DataFrame(paired_rows).to_csv(
    out / "JAK3_donor_paired_values.csv", index=False
)

print("\nGSE124219 ANALYSIS COMPLETE")
print("Matrix shape:", df.shape)
print("Donor pairs:", len(donors))
print("\n", result.to_string(index=False))
print("\nResults:", out)
print("Figures:", figdir)
print("Genome-wide FDR NOT calculated.")
