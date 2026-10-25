from pathlib import Path
from collections import Counter
from math import log2
import csv

# Define input and output paths
BASE = Path(__file__).resolve().parent
INPUT = BASE / "results/tables/antigen_region_mapping.csv"
OUT = BASE / "results/tables"
OUT.mkdir(parents=True, exist_ok=True)

# Load previously mapped antigen regions
with open(INPUT, newline="", encoding="utf-8") as file:
    rows = list(csv.DictReader(file))

proteins = ["EBNA1", "EBNA3A", "BZLF1"]

summary = []
position_results = []

# Analyze conservation for each protein
for protein in proteins:

    protein_rows = [
        row for row in rows
        if row["protein"] == protein
    ]

    # Select fully mapped, unambiguous antigen regions
    eligible = [
        row for row in protein_rows
        if row["region_status"] in ("EXACT_MATCH", "VARIANT")
    ]

    # Stop if no comparable antigen regions are available
    if not eligible:
     raise ValueError(
        f"No eligible antigen regions found for {protein}"
    )

    reference = eligible[0]["reference_region"]

    # Verify reference consistency and peptide lengths
    if any(
        row["reference_region"] != reference
        or len(row["observed_region"]) != len(reference)
        for row in eligible
    ):
        raise ValueError(
            f"Inconsistent peptide regions for {protein}"
        )

    # Count all eligible NCBI records
    total_records = len(eligible)

    exact_records = sum(
        row["observed_region"] == reference
        for row in eligible
    )

    record_match_percent = 100 * exact_records / total_records

    # Deduplicate identical full protein sequences
    unique_rows = {
        row["alignment_id"]: row
        for row in eligible
    }

    unique_sequences = list(unique_rows.values())

    exact_unique = sum(
        row["observed_region"] == reference
        for row in unique_sequences
    )

    unique_match_percent = (
        100 * exact_unique / len(unique_sequences)
    )

    # Save exact-match conservation statistics
    summary.append({
        "protein": protein,
        "reference_region": reference,
        "total_records": len(protein_rows),
        "eligible_records": total_records,
        "exact_records": exact_records,
        "record_match_percent": round(record_match_percent, 2),
        "unique_full_sequences": len(unique_sequences),
        "exact_unique_sequences": exact_unique,
        "unique_match_percent": round(unique_match_percent, 2)
    })

    # Calculate amino-acid frequencies at each position
    for i, reference_aa in enumerate(reference):

        record_counts = Counter(
            row["observed_region"][i]
            for row in eligible
        )

        unique_counts = Counter(
            row["observed_region"][i]
            for row in unique_sequences
        )

        # Calculate Shannon entropy using unique protein sequences
        n_unique = len(unique_sequences)

        entropy = -sum(
            (count / n_unique) * log2(count / n_unique)
            for count in unique_counts.values()
        )

        # Store every observed amino acid at this position
        for aa, count in sorted(unique_counts.items()):

            position_results.append({
                "protein": protein,
                "position_in_region": i + 1,
                "reference_amino_acid": reference_aa,
                "observed_amino_acid": aa,
                "record_count": record_counts.get(aa, 0),
                "unique_sequence_count": count,
                "unique_frequency_percent": round(
                    100 * count / n_unique, 2
                ),
                "shannon_entropy_bits": round(entropy, 4)
            })

    # Display conservation summary
    print(f"\n--- {protein} ---")
    print("Eligible records:", total_records)
    print("Exact matches:", exact_records)
    print("Record match (%):", round(record_match_percent, 2))
    print("Unique protein sequences:", len(unique_sequences))
    print("Unique sequence match (%):",
          round(unique_match_percent, 2))

# Save the conservation summary table
with open(OUT / "conservation_summary.csv",
          "w", newline="", encoding="utf-8") as file:

    writer = csv.DictWriter(
        file, fieldnames=list(summary[0].keys())
    )
    writer.writeheader()
    writer.writerows(summary)

# Save position-level amino-acid frequencies
with open(OUT / "position_conservation.csv",
          "w", newline="", encoding="utf-8") as file:

    writer = csv.DictWriter(
        file, fieldnames=list(position_results[0].keys())
    )
    writer.writeheader()
    writer.writerows(position_results)

print("\nConservation analysis completed!")
