
from pathlib import Path
from collections import Counter
import csv

# Define project paths
BASE = Path(__file__).resolve().parent
ALIGNMENTS = BASE / "results" / "alignments"
MAPPING = BASE / "results" / "tables" / "antigen_region_mapping.csv"

PROTEINS = ["EBNA1", "EBNA3A", "BZLF1"]

# Summarize initial protein identity screening
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

# Identify the mapping-status column from actual values
known_statuses = {
    "EXACT_MATCH",
    "VARIANT",
    "INCOMPLETE_REGION"
}

status_columns = [
    column
    for column in mapping[0]
    if any(
        row[column] in known_statuses
        for row in mapping
    )
]

print("\n--- Region Mapping ---")
print("CSV columns:", list(mapping[0]))

# Display mapping status frequencies by protein
if len(status_columns) == 1:
    status_col = status_columns[0]

    for protein in PROTEINS:
        rows = [
            row for row in mapping
            if row.get("protein") == protein
        ]

        print(f"\n{protein}:")
        print("Total mapped records:", len(rows))

        counts = Counter(row[status_col] for row in rows)

        for status, count in counts.most_common():
            print(f"{status}: {count}")

else:
    print("Mapping status column needs manual review.")
    print("Candidate columns:", status_columns)
