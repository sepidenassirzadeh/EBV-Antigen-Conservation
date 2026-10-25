from pathlib import Path
from collections import Counter
import csv

# Define project paths
BASE = Path(__file__).resolve().parent
ALIGNMENTS = BASE / "results" / "alignments"
MAPPING = BASE / "results" / "tables" / "antigen_region_mapping.csv"

PROTEINS = ["EBNA1", "EBNA3A", "BZLF1"]

# Summarize protein identity screening results
for protein in PROTEINS:

    audit_file = ALIGNMENTS / f"{protein}_identity_audit.csv"

    with open(audit_file, newline="", encoding="utf-8") as file:
        rows = list(csv.DictReader(file))

    counts = Counter(row["status"] for row in rows)

    print(f"\n--- {protein}: Identity QC ---")
    print("Total records:", len(rows))

    for status, count in counts.most_common():
        print(f"{status}: {count}")

# Load antigen-region mapping results
with open(MAPPING, newline="", encoding="utf-8") as file:
    mapping = list(csv.DictReader(file))

# Verify that the expected status column exists
if not mapping or "region_status" not in mapping[0]:
    raise ValueError("Missing region_status column in mapping table.")

# Display region-mapping statuses for each protein
print("\n--- Region Mapping ---")

for protein in PROTEINS:

    protein_rows = [
        row for row in mapping
        if row["protein"] == protein
    ]

    counts = Counter(
        row["region_status"] for row in protein_rows
    )

    print(f"\n{protein}:")
    print("Total mapped records:", len(protein_rows))

    for status, count in counts.most_common():
        print(f"{status}: {count}")
