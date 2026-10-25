
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from pathlib import Path
from collections import Counter
import csv

# Define project files
BASE = Path(__file__).resolve().parent
RAW = BASE / "data/raw/BZLF1.fasta"
AUDIT = BASE / "results/alignments/BZLF1_identity_audit.csv"
OUTPUT = BASE / "results/tables/BZLF1_identity_review.csv"

# Load the validated BZLF1 reference
reference = SeqIO.read(BASE / "BZLF1.fasta", "fasta")
ref_seq = str(reference.seq).upper()

# Load downloaded sequences and previous QC results
records = SeqIO.to_dict(SeqIO.parse(RAW, "fasta"))

with open(AUDIT, newline="", encoding="utf-8") as file:
    flagged = [
        row for row in csv.DictReader(file)
        if row["status"] == "REVIEW_IDENTITY"
    ]

# Configure local protein alignment
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5

results = []

# Compare each flagged protein with the BZLF1 reference
for row in flagged:
    accession = row["accession"]
    sequence = str(records[accession].seq).upper()

    alignment = aligner.align(ref_seq, sequence)[0]

    # Get aligned amino-acid pairs
    pairs = [
        (int(i), int(j))
        for i, j in zip(*alignment.indices)
        if i >= 0 and j >= 0
    ]

    # Calculate reference coverage and amino-acid identity
    coverage = len({i for i, j in pairs}) / len(ref_seq)

    identity = (
        sum(ref_seq[i] == sequence[j] for i, j in pairs)
        / len(pairs)
        if pairs else 0
    )

    # Flag strong matches for further identity verification
    status = (
        "LIKELY_BZLF1"
        if coverage >= 0.80 and identity >= 0.80
        else "REVIEW_IDENTITY"
    )

    results.append({
        "accession": accession,
        "description": row["description"],
        "length": len(sequence),
        "reference_coverage_percent": round(coverage * 100, 2),
        "aligned_identity_percent": round(identity * 100, 2),
        "screening_status": status
    })

# Save screening results without deleting any records
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with open(OUTPUT, "w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(file, fieldnames=list(results[0].keys()))
    writer.writeheader()
    writer.writerows(results)

# Display the screening summary
print("Total reviewed:", len(results))
print("Screening:", dict(Counter(
    row["screening_status"] for row in results
)))

# Show examples of potentially misclassified BZLF1 records
print("\nPotential BZLF1 records:")
for row in results:
    if row["screening_status"] == "LIKELY_BZLF1":
        print(
            row["accession"],
            row["description"][:65],
            row["reference_coverage_percent"],
            row["aligned_identity_percent"]
        )
        if sum(r["screening_status"] == "LIKELY_BZLF1"
               for r in results[:results.index(row)+1]) >= 5:
            break

print("\nIdentity review saved!")
