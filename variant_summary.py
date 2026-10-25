
from pathlib import Path
from collections import Counter
import csv

# Define project paths
BASE = Path(__file__).resolve().parent
input_file = BASE / "results/tables/antigen_region_mapping.csv"
output_file = BASE / "results/tables/variant_summary.csv"

# Load mapped antigen regions
with open(input_file, newline="", encoding="utf-8") as file:
    rows = list(csv.DictReader(file))

proteins = ["EBNA1", "EBNA3A", "BZLF1"]
summary = []

# Summarize peptide variants for each protein
for protein in proteins:

    # Select fully comparable antigen regions
    eligible = [
        row for row in rows
        if row["protein"] == protein
        and row["region_status"] in ("EXACT_MATCH", "VARIANT")
    ]

    # Ensure comparable records are available
    if not eligible:
        raise ValueError(
            f"No comparable antigen regions found for {protein}"
        )

    # Identify the reference antigen peptide
    reference = eligible[0]["reference_region"]

    # Verify peptide lengths and reference consistency
    if any(
        len(row["observed_region"]) != len(reference)
        or row["reference_region"] != reference
        for row in eligible
    ):
        raise ValueError(
            f"Inconsistent antigen regions for {protein}"
        )

    # Count the frequency of each observed peptide
    peptide_counts = Counter(
        row["observed_region"] for row in eligible
    )

    print(f"\n--- {protein} ---")
    print("Currently mapped records:", len(eligible))
    print("Distinct antigen peptides:", len(peptide_counts))

    # Identify amino acid substitutions in each peptide
    for peptide, count in peptide_counts.most_common():

        differences = [
            f"{reference[i]}{i + 1}{peptide[i]}"
            for i in range(len(reference))
            if reference[i] != peptide[i]
        ]

        change = ", ".join(differences) if differences else "Reference"

        print(f"{peptide}: {count} records | {change}")

        summary.append({
            "protein": protein,
            "peptide": peptide,
            "record_count": count,
            "difference_from_reference": change
        })

# Save the peptide frequency summary
with open(output_file, "w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(
        file,
        fieldnames=[
            "protein",
            "peptide",
            "record_count",
            "difference_from_reference"
        ]
    )

    writer.writeheader()
    writer.writerows(summary)

print("\nVariant summary saved successfully!")
