from pathlib import Path
import csv

# Load the mapped antigen regions
BASE = Path(__file__).resolve().parent
file_path = BASE / "results/tables/antigen_region_mapping.csv"

with open(file_path, newline="", encoding="utf-8") as file:
    rows = list(csv.DictReader(file))

# Find EBNA3A records with the unusual peptide
for row in rows:
    if (row["protein"] == "EBNA3A"
            and row["observed_region"] == "LLRGRAYGQ"):

        print("Accession:", row["accession"])
        print("Alignment ID:", row["alignment_id"])
        print("Observed peptide:", row["observed_region"])
        print("--------------------")

from Bio import SeqIO
from pathlib import Path

# Load the original EBNA3A protein sequences
base = Path(__file__).resolve().parent
fasta_file = base / "data/raw/EBNA3A.fasta"

records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# Inspect the two unusual EBNA3A records
target_ids = ["ATO60671.1", "ATO60382.1"]

for accession in target_ids:
    record = records[accession]

    print("\nAccession:", record.id)
    print("Description:", record.description)
    print("Length:", len(record.seq))
    print("Peptide position:", str(record.seq).find("LLRGRAYGQ"))

# Check whether both records contain identical full sequences
same_sequence = (
    records[target_ids[0]].seq == records[target_ids[1]].seq
)

print("\nIdentical protein sequences:", same_sequence)