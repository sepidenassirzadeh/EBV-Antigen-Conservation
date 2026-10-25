
from pathlib import Path
from collections import Counter
import csv

# Locate the antigen-region mapping table
BASE = Path(__file__).resolve().parent
input_file = BASE / "results/tables/antigen_region_mapping.csv"
output_file = BASE / "results/tables/variant_summary.csv"

# Load the mapped antigen regions
with open(input_file, newline="", encoding="utf-8") as file:
    rows = list(csv.DictReader(file))

proteins = ["EBNA1", "EBNA3A", "BZLF1"]
summary = []

# Summarize exact matches and variants for each protein
for protein in proteins:

    # Select records with fully mapped, unambiguous regions
    eligible = [
        row for row in rows
        if row["protein"] == protein
        and row["region_status"] in ("EXACT_MATCH", "VARIANT")
    ]

    # Count the frequency of each observed peptide
    peptide_counts = Counter(
        row["observed_region"] for row in eligible
    )

    print(f"\n--- {protein} ---")
    print("Currently mapped records:", len(eligible))
    print("Distinct antigen peptides:", len(peptide_counts))

    for peptide, count in peptide_counts.most_common():

        # Identify substitutions relative to the reference
        reference = eligible[0]["reference_region"]

        differences = [
            f"{reference[i]}{i+1}{peptide[i]}"
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

# Save peptide frequency summary
with open(output_file, "w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(
        file,
        fieldnames=[
            "protein", "peptide",
            "record_count", "difference_from_reference"
        ]
    )
    writer.writeheader()
    writer.writerows(summary)

print("\nVariant summary saved successfully!")
