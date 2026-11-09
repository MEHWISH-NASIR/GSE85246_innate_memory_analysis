from pathlib import Path
from zipfile import ZipFile
import pandas as pd

archive = Path("data/atac/GSE323311/Data.zip")
out = Path("05_results/GSE323311")
prefix = "Data/DESeq/ReadCutoff20/"

chrom, left, right = "chr19", 17819782, 17852982
donors = {"A": "ND410b", "B": "ND650", "C": "ND500b"}
records = []
bed_summary = []

with ZipFile(archive) as z:
    assert z.testzip() is None, "ZIP integrity failed"

    for condition in ("bg", "mdp", "ifn"):
        for group, donor in donors.items():
            filename = f"{condition}_res_{group}.csv"

            with z.open(prefix + filename) as f:
                df = pd.read_csv(f)

            ids = df.iloc[:, 0].astype(str)
            coords = ids.str.extract(
                r"^(chr[^\s]+)\s+(\d+)-(\d+)$"
            )
            assert coords.notna().all().all(), filename

            starts = coords[1].astype(int)
            ends = coords[2].astype(int)

            hits = df[
                (coords[0] == chrom) &
                (starts < right) &
                (ends > left)
            ]

            for _, row in hits.iterrows():
                records.append({
                    "dataset": "GSE323311",
                    "genome": "hg38",
                    "day": 6,
                    "condition": condition.upper(),
                    "donor": donor,
                    "peak": row.iloc[0],
                    "baseMean": row["baseMean"],
                    "log2FC": row["log2FoldChange"],
                    "pvalue": row["pvalue"],
                    "padj": row["padj"],
                    "FDR_status": (
                        "unavailable" if pd.isna(row["padj"])
                        else "significant" if row["padj"] < 0.05
                        else "not_significant"
                    )
                })

        for suffix in (
            "sig_two.bed",
            "sig_two_UP.bed",
            "sig_two_DOWN.bed"
        ):
            filename = f"{condition}_{suffix}"

            with z.open(prefix + filename) as f:
                bed = pd.read_csv(
                    f, sep=r"\s+", header=None, comment="#"
                )

            assert bed.shape[1] >= 3, filename
            starts = pd.to_numeric(bed.iloc[:, 1])
            ends = pd.to_numeric(bed.iloc[:, 2])

            mask = (
                (bed.iloc[:, 0] == chrom) &
                (starts < right) &
                (ends > left)
            )

            bed_summary.append({
                "condition": condition.upper(),
                "file": filename,
                "total_regions": len(bed),
                "JAK3_window_overlaps": int(mask.sum())
            })

results = pd.DataFrame(records)
summary = pd.DataFrame(bed_summary)

assert len(results) == 18, (
    f"Expected 18 JAK3 result rows; found {len(results)}"
)
assert len(summary) == 9

results.to_csv(
    out / "JAK3_3donor_DESeq_results.csv",
    index=False
)
summary.to_csv(
    out / "JAK3_shared_region_audit.csv",
    index=False
)

print("\nGSE323311 CROSS-DONOR AUDIT COMPLETE")
print("Donor-specific JAK3 rows:", len(results))
print("FDR-significant:", (results.FDR_status == "significant").sum())
print("Missing FDR:", (results.FDR_status == "unavailable").sum())
print("\nShared BED results:")
print(summary.to_string(index=False))
print("\nResults saved:", out)
print("Original data unchanged.")
