
from Bio import SeqIO
from pathlib import Path

# Load the selected EBNA3A reference protein
BASE = Path(__file__).resolve().parent
record = SeqIO.read(BASE / "EBNA-3A.fasta", "fasta")

sequence = str(record.seq).upper()

# Verify the reference accession and length
print("Reference:", record.id)
print("Protein length:", len(sequence))

# Extract the two reported coordinate ranges
print("\nPositions 193-201:", sequence[192:201])
print("Positions 325-333:", sequence[324:333])

# Locate the published and reference peptide variants
for peptide in ["FLRGRAYGL", "FLRGRAYGI"]:
    position = sequence.find(peptide)

    if position == -1:
        print(f"{peptide}: Not found")
    else:
        print(
            f"{peptide}: positions "
            f"{position + 1}-{position + len(peptide)}"
        )

# Check the amino acid substitution
reference_peptide = sequence[324:333]
published_peptide = "FLRGRAYGL"

differences = [
    (i + 1, a, b)
    for i, (a, b) in enumerate(
        zip(reference_peptide, published_peptide)
    )
    if a != b
]

print("\nDifferences:", differences)
