import csv
import subprocess
from pathlib import Path
from urllib.parse import urlparse

manifest = Path(
    "results/12_kinase_reevaluation/"
    "12a_day6_chromatin_file_manifest.csv"
)

root = Path("data/chip/day6_formal")

if not manifest.exists():
    raise SystemExit(f"Missing manifest: {manifest}")

with manifest.open(newline="") as fh:
    rows = list(csv.DictReader(fh))

print(f"Files in manifest: {len(rows)}")

for i, row in enumerate(rows, start=1):

    mark = row["Mark"]
    condition = row["Condition"]
    url = row["URL"]

    source_name = Path(
        urlparse(url).path
    ).name

    suffix = ".NotNormalized.bw" if "NotNormalized" in source_name else ".Normalized.bw"

    outfile = (
        root /
        mark /
        f"{condition}_{mark}{suffix}"
    )

    outfile.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    print(
        f"\n[{i:02d}/{len(rows):02d}] "
        f"{mark} {condition}"
    )

    print(f"URL: {url}")
    print(f"OUT: {outfile}")

    if outfile.exists() and outfile.stat().st_size > 0:
        print("Already present — skipping.")
        continue

    cmd = [
        "curl",
        "-L",
        "--fail",
        "--retry", "3",
        "--retry-delay", "5",
        "-o", str(outfile),
        url
    ]

    subprocess.run(
        cmd,
        check=True
    )

print("\nAll downloads completed.")
