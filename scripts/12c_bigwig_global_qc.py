from pathlib import Path
import csv
import sys

try:
    import pyBigWig
except ImportError:
    raise SystemExit(
        "pyBigWig is not installed.\n"
        "Run: pip install pyBigWig"
    )

root = Path("data/chip/day6_formal")
out = Path(
    "results/12_kinase_reevaluation/"
    "12c_bigwig_global_qc.csv"
)

files = sorted(root.rglob("*.bw"))

print("\n========================================")
print("12c — BIGWIG GLOBAL QC")
print("========================================\n")

print("BigWigs found:", len(files))

if len(files) != 18:
    raise SystemExit(
        f"Expected 18 BigWigs but found {len(files)}"
    )

rows = []

for i, f in enumerate(files, 1):

    print(f"[{i:02d}/18] {f}")

    try:
        bw = pyBigWig.open(str(f))

        if bw is None:
            raise RuntimeError("pyBigWig could not open file")

        h = bw.header()

        n_bases = h["nBasesCovered"]
        sum_data = h["sumData"]

        mean_covered = (
            sum_data / n_bases
            if n_bases > 0
            else float("nan")
        )

        chroms = bw.chroms()
        genome_size = sum(chroms.values())

        genome_mean_including_zero = (
            sum_data / genome_size
            if genome_size > 0
            else float("nan")
        )

        rows.append({
            "file": str(f),
            "mark": f.parent.name,
            "sample": f.name,
            "normalization_label":
                "NotNormalized"
                if "NotNormalized" in f.name
                else "Normalized",
            "file_size_MB":
                round(f.stat().st_size / 1024**2, 3),
            "nBasesCovered": n_bases,
            "sumData": sum_data,
            "mean_over_covered_bases":
                mean_covered,
            "genome_size": genome_size,
            "genome_mean_including_zero":
                genome_mean_including_zero,
            "minVal": h["minVal"],
            "maxVal": h["maxVal"]
        })

        bw.close()

    except Exception as e:
        print("ERROR:", e)
        sys.exit(1)

out.parent.mkdir(
    parents=True,
    exist_ok=True
)

with out.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=rows[0].keys()
    )
    writer.writeheader()
    writer.writerows(rows)

print("\nAll 18 BigWigs opened successfully.")

print("\n========================================")
print("GLOBAL SIGNAL SUMMARY")
print("========================================\n")

for r in rows:
    print(
        f"{r['mark']:8s} "
        f"{r['sample'][:34]:34s} "
        f"{r['normalization_label']:13s} "
        f"sum={r['sumData']:.3f} "
        f"globalMean={r['genome_mean_including_zero']:.6f}"
    )

print("\nSaved:")
print(out)

print("\n========================================")
print("12c COMPLETE")
print("========================================")
