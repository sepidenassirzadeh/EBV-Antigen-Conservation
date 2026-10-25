
from Bio import SeqIO
from pathlib import Path
from collections import Counter
import csv

# Define project folders
BASE = Path(__file__).resolve().parent
ALIGN_DIR = BASE / "results" / "alignments"
OUT_DIR = BASE / "results" / "tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Define reference coordinates and expected peptides
targets = {
    "EBNA1": ("YP_401677.1", 386, 405, "SQSSSSGSPPRRPPPGRRPF"),
    "EBNA3A": ("YP_401669.1", 325, 333, "FLRGRAYGI"),
    "BZLF1": ("YP_401673.1", 54, 64, "EPLPQGQLTAY")
}

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
all_results = []

# Process each protein alignment
for protein, (ref_id, start, end, expected) in targets.items():

    # Load the screened multiple sequence alignment
    msa_path = ALIGN_DIR / f"{protein}_screened_aligned.fasta"
    aligned = {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(msa_path, "fasta")
    }

    # Check that all aligned sequences have equal length
    if len(set(map(len, aligned.values()))) != 1:
        raise ValueError(f"Invalid alignment length: {protein}")

    # Identify reference positions in alignment coordinates
    reference = aligned[ref_id]
    columns = []
    position = 0

    for column, aa in enumerate(reference):
        if aa != "-":
            position += 1

            if start <= position <= end:
                columns.append(column)

    # Verify the mapped reference antigen region
    ref_region = "".join(reference[i] for i in columns)

    if ref_region != expected:
        raise ValueError(f"Reference region mismatch: {protein}")

    # Load the previous protein-identity screening results
    audit_path = ALIGN_DIR / f"{protein}_identity_audit.csv"

    with open(audit_path, newline="", encoding="utf-8") as file:
        audit = list(csv.DictReader(file))

    # Extract the target region for every included accession
    for row in audit:

        alignment_id = row["alignment_id"]
        original_status = row["status"]
        region = ""

        if original_status != "INCLUDED":
            region_status = "NOT_ASSESSED_" + original_status

        else:
            if alignment_id not in aligned:
                raise ValueError(f"Missing alignment ID: {alignment_id}")

            sequence = aligned[alignment_id]
            region = "".join(sequence[i] for i in columns)

            # Classify gaps, ambiguity, and sequence variation
            if "-" in region:
                region_status = "INCOMPLETE_REGION"
            elif any(aa not in STANDARD_AA for aa in region):
                region_status = "AMBIGUOUS_REGION"
            elif region == expected:
                region_status = "EXACT_MATCH"
            else:
                region_status = "VARIANT"

        # Save one result per original NCBI accession
        all_results.append({
            "protein": protein,
            "accession": row["accession"],
            "alignment_id": alignment_id,
            "reference_region": expected,
            "observed_region": region,
            "identity_screen": original_status,
            "region_status": region_status
        })

    # Display the screening summary
    protein_rows = [
        row for row in all_results
        if row["protein"] == protein
    ]

    print(f"\n{protein}:")
    print(dict(Counter(
        row["region_status"] for row in protein_rows
    )))

# Save the antigen-region mapping table
output = OUT_DIR / "antigen_region_mapping.csv"

with open(output, "w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(
        file,
        fieldnames=list(all_results[0].keys())
    )
    writer.writeheader()
    writer.writerows(all_results)

print("\nMapping table saved:", output)

