
from pathlib import Path
from Bio import SeqIO
import csv

# Locate the existing BZLF1 files
BASE_DIR = Path(__file__).resolve().parent

raw_file = BASE_DIR / "data/raw/BZLF1.fasta"
msa_file = BASE_DIR / "results/alignments/BZLF1_aligned.fasta"
metadata_file = BASE_DIR / "metadata/BZLF1_inventory.csv"

# Load original sequences and source descriptions
records = list(SeqIO.parse(raw_file, "fasta"))

with open(metadata_file, encoding="utf-8", newline="") as file:
    metadata = {
        row["accession"]: row["description"]
        for row in csv.DictReader(file)
    }

# Check the distribution of original protein lengths
lengths = [len(record.seq) for record in records]

print("Total records:", len(records))
print("Shortest sequence:", min(lengths))
print("Longest sequence:", max(lengths))
print("Sequences longer than 300 aa:", sum(n > 300 for n in lengths))
print("Sequences longer than 500 aa:", sum(n > 500 for n in lengths))

# Display the longest records and their annotations
print("\nLongest records:")

for record in sorted(records, key=lambda r: len(r.seq), reverse=True)[:5]:
    print(record.id, len(record.seq))
    print(metadata.get(record.id, "Description unavailable")[:120])

# Inspect the aligned reference and gap count
aligned = list(SeqIO.parse(msa_file, "fasta"))

ref = next(
    record for record in aligned
    if record.id == "YP_401673.1"
)

print("\nMSA length:", len(ref.seq))
print("Reference gaps:", str(ref.seq).count("-"))
print("Reference residues:", len(str(ref.seq).replace("-", "")))
