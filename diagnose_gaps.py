
from Bio import SeqIO
from pathlib import Path
from collections import Counter
import csv

# Define paths and reference positions
BASE = Path(__file__).resolve().parent

TARGETS = {
    "EBNA1": ("YP_401677.1", 386, 405),
    "EBNA3A": ("YP_401669.1", 325, 333),
    "BZLF1": ("YP_401673.1", 54, 64),
}

# Load the existing region mapping table
mapping_file = BASE / "results/tables/antigen_region_mapping.csv"

with open(mapping_file, newline="", encoding="utf-8") as file:
    mapping = list(csv.DictReader(file))

# Examine alignment gaps for each protein
for protein, (ref_id, start, end) in TARGETS.items():

    alignment_file = (
        BASE / "results/alignments"
        / f"{protein}_screened_aligned.fasta"
    )

    alignment = SeqIO.to_dict(
        SeqIO.parse(alignment_file, "fasta")
    )

    reference = str(alignment[ref_id].seq)

    # Find alignment columns corresponding to reference positions
    columns = []
    reference_position = 0

    for col, amino_acid in enumerate(reference):
        if amino_acid != "-":
            reference_position += 1

            if start <= reference_position <= end:
                columns.append(col)

    assert len(columns) == end - start + 1

    counts = Counter()
    examples = []

    # Examine records previously classified as incomplete
    incomplete = [
        row for row in mapping
        if row["protein"] == protein
        and row["region_status"] == "INCOMPLETE_REGION"
    ]

    for row in incomplete:
        alignment_id = row["alignment_id"]
        seq = str(alignment[alignment_id].seq)

        # Identify the first and last aligned residues
        occupied = [
            i for i, amino_acid in enumerate(seq)
            if amino_acid != "-"
        ]

        if not occupied:
            category = "NO_ALIGNED_RESIDUES"

        else:
            first = occupied[0]
            last = occupied[-1]

            # Separate terminal gaps from internal gaps
            terminal_gaps = [
                col for col in columns
                if seq[col] == "-"
                and (col < first or col > last)
            ]

            internal_gaps = [
                col for col in columns
                if seq[col] == "-"
                and first <= col <= last
            ]

            if terminal_gaps and internal_gaps:
                category = "MIXED_GAPS"
            elif terminal_gaps:
                category = "TERMINAL_ONLY"
            elif internal_gaps:
                category = "INTERNAL_ONLY"
            else:
                category = "NO_TARGET_GAPS"

        counts[category] += 1

        # Keep a few examples for inspection
        if len(examples) < 5:
            examples.append((
                row["accession"],
                category,
                "".join(seq[col] for col in columns)
            ))

    # Display the diagnostic summary
    print(f"\n--- {protein} ---")
    print("Incomplete-region records:", len(incomplete))

    for category, count in counts.most_common():
        print(f"{category}: {count}")

    print("Examples:")
    for accession, category, peptide in examples:
        print(accession, category, peptide)
