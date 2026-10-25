
from pathlib import Path
from collections import Counter
import csv
import matplotlib.pyplot as plt
import numpy as np

# Define project paths
BASE = Path(__file__).resolve().parent
ALIGNMENTS = BASE / "results" / "alignments"
MAPPING = BASE / "results" / "tables" / "antigen_region_mapping.csv"
FIGURES = BASE / "results" / "figures"

FIGURES.mkdir(parents=True, exist_ok=True)

proteins = ["EBNA1", "EBNA3A", "BZLF1"]

# Load the antigen-region mapping table
with open(MAPPING, newline="", encoding="utf-8") as file:
    mapping = list(csv.DictReader(file))

summary = []

# Count records at each analysis stage
for protein in proteins:

    audit_file = ALIGNMENTS / f"{protein}_identity_audit.csv"

    with open(audit_file, newline="", encoding="utf-8") as file:
        audit = list(csv.DictReader(file))

    region_rows = [
        row for row in mapping
        if row["protein"] == protein
    ]

    region_counts = Counter(
        row["region_status"] for row in region_rows
    )

    initial = len(audit)

    screened = sum(
        row["status"] == "INCLUDED"
        for row in audit
    )

    comparable = (
        region_counts["EXACT_MATCH"]
        + region_counts["VARIANT"]
    )

    exact = region_counts["EXACT_MATCH"]

    # Verify that all records are represented
    assert len(region_rows) == initial, protein
    assert comparable <= screened <= initial, protein

    summary.append([initial, screened, comparable, exact])

# Plot record counts across the analysis stages
stages = [
    "Retrieved",
    "Identity/QC passed",
    "Region comparable",
    "Exact reference match"
]

x = np.arange(len(stages))
width = 0.24

fig, ax = plt.subplots(figsize=(11, 6))

for i, protein in enumerate(proteins):
    counts = summary[i]
    positions = x + (i - 1) * width

    bars = ax.bar(
        positions,
        counts,
        width,
        label=protein
    )

    ax.bar_label(bars, padding=3, fontsize=9)

# Format and save the QC flow figure
ax.set_xticks(x)
ax.set_xticklabels(stages)
ax.set_ylabel("Number of accession records")
ax.set_title("EBV Antigen Analysis:Record Screening and Reference Matching ")
ax.legend()
ax.set_ylim(0, max(row[0] for row in summary) * 1.15)

fig.tight_layout()

output = FIGURES / "figure3_qc_flow.png"
fig.savefig(output, dpi=300)
plt.close(fig)

print("QC flow figure saved:", output)

# Display summary counts
for protein, counts in zip(proteins, summary):
    print(protein, dict(zip(stages, counts)))
