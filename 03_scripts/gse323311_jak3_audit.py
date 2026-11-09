from pathlib import Path
from zipfile import ZipFile
import pandas as pd
import re

root = Path(".")
archive = root / "data/atac/GSE323311/Data.zip"
out = root / "05_results/GSE323311"
out.mkdir(parents=True, exist_ok=True)

prefix = "Data/DESeq/D6vD11/"

# GRCh38 JAK3 region with 5-kb flanks
chrom = "chr19"
start = 17819782
end = 17852982

records = []
retained_records = []

with ZipFile(archive) as z:
    for condition in ["bg", "mdp", "ifn"]:
        retained_name = prefix + f"{condition}_retained.bed"

        with z.open(retained_name) as f:
            bed = pd.read_csv(
                f, sep=r"\s+", header=None,
                comment="#"
            )

        bed = bed.iloc[:, :3].copy()
        bed.columns = ["chrom", "start", "end"]
        bed["start"] = pd.to_numeric(bed["start"])
        bed["end"] = pd.to_numeric(bed["end"])

        retained = bed[
            (bed["chrom"] == chrom) &
            (bed["start"] < end) &
            (bed["end"] > start)
        ]

        print(
            f"{condition.upper()} retained JAK3 overlaps:",
            len(retained)
        )

        for _, r in retained.iterrows():
            retained_records.append({
                "condition": condition.upper(),
                "chrom": r["chrom"],
                "start": int(r["start"]),
                "end": int(r["end"])
            })

        for day in [6, 11]:
            name = prefix + f"{condition}_res_d{day}.csv"

            with z.open(name) as f:
                df = pd.read_csv(f)

            peak_col = df.columns[0]
            coords = df[peak_col].astype(str).str.extract(
                r"^(chr[^\s]+)\s+(\d+)-(\d+)$"
            )
            coords.columns = ["chrom", "start", "end"]

            if coords.isna().any().any():
                raise ValueError(
                    f"Unparsed genomic coordinates in {name}"
                )

            coords["start"] = coords["start"].astype(int)
            coords["end"] = coords["end"].astype(int)

            mask = (
                (coords["chrom"] == chrom) &
                (coords["start"] < end) &
                (coords["end"] > start)
            )

            hits = df.loc[mask].copy()
            hit_coords = coords.loc[mask]

            print(
                f"{condition.upper()} Day {day}:",
                len(hits), "candidate peaks"
            )

            for idx, row in hits.iterrows():
                c = hit_coords.loc[idx]

                overlap_retained = (
                    (retained["chrom"] == c["chrom"]) &
                    (retained["start"] < c["end"]) &
                    (retained["end"] > c["start"])
                ).any()

                padj = row["padj"]

                if pd.isna(padj):
                    status = "padj_missing"
                elif padj < 0.05:
                    status = "FDR_significant"
                else:
                    status = "not_FDR_significant"

                records.append({
                    "dataset": "GSE323311",
                    "donor": "ND410b",
                    "genome": "hg38",
                    "condition": condition.upper(),
                    "day": day,
                    "peak": row[peak_col],
                    "chrom": c["chrom"],
                    "start": int(c["start"]),
                    "end": int(c["end"]),
                    "baseMean": row["baseMean"],
                    "log2FoldChange": row["log2FoldChange"],
                    "lfcSE": row["lfcSE"],
                    "pvalue": row["pvalue"],
                    "padj": padj,
                    "significance_status": status,
                    "overlaps_retained_bed": bool(
                        overlap_retained
                    )
                })

result = pd.DataFrame(records)
if result.empty:
    raise RuntimeError("No JAK3 candidate peaks found")

result = result.sort_values(
    ["condition", "peak", "day"]
)

result.to_csv(
    out / "JAK3_D6_D11_DESeq_results.csv",
    index=False
)

pd.DataFrame(
    retained_records,
    columns=["condition", "chrom", "start", "end"]
).to_csv(
    out / "JAK3_retained_overlaps.csv",
    index=False
)

summary = result.groupby(
    ["condition", "day"]
).agg(
    candidate_peaks=("peak", "nunique"),
    FDR_significant=(
        "significance_status",
        lambda x: (x == "FDR_significant").sum()
    ),
    missing_padj=("padj", lambda x: x.isna().sum()),
    retained_overlaps=(
        "overlaps_retained_bed", "sum"
    )
).reset_index()

summary.to_csv(
    out / "JAK3_condition_summary.csv",
    index=False
)

print("\nGSE323311 JAK3 EXTRACTION COMPLETE")
print(summary.to_string(index=False))
print("\nFiles saved to:", out)
print("Original ZIP unchanged.")
