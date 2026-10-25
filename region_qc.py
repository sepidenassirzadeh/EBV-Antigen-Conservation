from Bio import SeqIO
from Bio.Align import PairwiseAligner
from pathlib import Path
import csv

# Define project folders
BASE_DIR = Path(__file__).resolve().parent
RAW_DIR = BASE_DIR / "data" / "raw"
META_DIR = BASE_DIR / "metadata"
META_DIR.mkdir(exist_ok=True)

# Define reference sequences and target regions
targets = {
    "EBNA1": {
        "file": "EBNA1.fasta",
        "id": "YP_401677.1",
        "start": 386,
        "end": 405
    },
    "EBNA3A": {
        "file": "EBNA-3A.fasta",
        "id": "YP_401669.1",
        "start": 325,
        "end": 333
    },
    "BZLF1": {
        "file": "BZLF1.fasta",
        "id": "YP_401673.1",
        "start": 54,
        "end": 64
    }
}

# Configure local pairwise alignment for the QC screening
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5

qc_results = []

# Process each target protein separately
for protein, info in targets.items():

    # Load and verify the reference sequence
    reference = SeqIO.read(
        BASE_DIR / info["file"], "fasta"
    )

    if reference.id != info["id"]:
        raise ValueError(f"Incorrect reference for {protein}")

    start = info["start"] - 1
    end = info["end"]
    expected_length = end - start

    reference_region = str(reference.seq[start:end])

    if len(reference_region) != expected_length:
        raise ValueError(f"Invalid reference region: {protein}")

    # Load downloaded candidate sequences
    records = list(SeqIO.parse(
        RAW_DIR / f"{protein}.fasta", "fasta"
    ))

    # Reuse results for identical full protein sequences
    cache = {}

    for record in records:
        sequence = str(record.seq).upper()

        if sequence not in cache:

            # Align each distinct sequence to its reference
            alignment = aligner.align(
                str(reference.seq), sequence
            )[0]

            # Map reference indices to candidate indices
            mapping = {
                int(ref_pos): int(query_pos)
                for ref_pos, query_pos in zip(
                    *alignment.indices
                )
                if ref_pos >= 0 and query_pos >= 0
            }

            # Extract aligned residues at reference positions
            positions = [
                mapping.get(pos)
                for pos in range(start, end)
            ]

            covered = sum(pos is not None for pos in positions)

            region = "".join(
                sequence[pos] if pos is not None else "-"
                for pos in positions
            )

            # Screen region coverage and unresolved residues
            if covered < expected_length:
                status = "INCOMPLETE_REGION"
            elif any(aa not in "ACDEFGHIKLMNPQRSTVWY"
                     for aa in region):
                status = "AMBIGUOUS_REGION"
            else:
                status = "REVIEW_COVERED"

            cache[sequence] = (
                covered, region, status, alignment.score
            )

        covered, region, status, score = cache[sequence]

        # Store the QC screening result for each accession
        qc_results.append({
            "protein": protein,
            "accession": record.id,
            "reference_id": reference.id,
            "reference_region": reference_region,
            "mapped_region": region,
            "covered_positions": covered,
            "total_positions": expected_length,
            "coverage_percent": round(
                100 * covered / expected_length, 2
            ),
            "alignment_score": round(score, 2),
            "qc_status": status
        })

    print(f"{protein}: checked {len(records)} records")

# Export the QC screening table
output_path = META_DIR / "region_qc.csv"

with open(
    output_path, "w", newline="", encoding="utf-8"
) as file:

    writer = csv.DictWriter(
        file, fieldnames=list(qc_results[0].keys())
    )
    writer.writeheader()
    writer.writerows(qc_results)

# Summarize the screening results
print("\nQC summary:")

for protein in targets:
    rows = [
        row for row in qc_results
        if row["protein"] == protein
    ]

    print(f"\n{protein}:")
    for status in [
        "REVIEW_COVERED",
        "INCOMPLETE_REGION",
        "AMBIGUOUS_REGION"
    ]:
        count = sum(
            row["qc_status"] == status for row in rows
        )
        print(f"  {status}: {count}")

print("\nQC table saved:", output_path)
