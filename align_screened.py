
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfamsa import Aligner, Sequence
from pathlib import Path
from collections import Counter
import csv
import re

# Define project directories
BASE = Path(__file__).resolve().parent
RAW = BASE / "data" / "raw"
META = BASE / "metadata"
OUT = BASE / "results" / "alignments"

OUT.mkdir(parents=True, exist_ok=True)

# Define reference proteins and identity patterns
targets = {
    "EBNA1": {
        "reference": "EBNA1.fasta",
        "id": "YP_401677.1",
        "pattern": r"\bebna[- ]?1\b|\bnuclear antigen[- ]?1\b",
        "region_length": 20
    },
    "EBNA3A": {
        "reference": "EBNA-3A.fasta",
        "id": "YP_401669.1",
        "pattern": r"\bebna[- ]?3a\b|\bebna[- ]?3\b|\bnuclear antigen[- ]?3a\b",
        "region_length": 9
    },
    "BZLF1": {
        "reference": "BZLF1.fasta",
        "id": "YP_401673.1",
        "pattern": r"\bbzlf[- ]?1\b|\bzebra\b|\bzta\b",
        "region_length": 11
    }
}

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Initialize the multiple sequence aligner
aligner = Aligner(threads=2)

# Load provisionally identified BZLF1 accessions
# Load BZLF1 accessions confirmed by GenBank gene annotation
validation_file = (
    BASE / "results/tables/BZLF1_gene_validation.csv"
)

with open(validation_file, newline="", encoding="utf-8") as file:
    rescued_ids = {
        row["accession"]
        for row in csv.DictReader(file)
        if row["validation_status"] == "GENE_CONFIRMED"
    }

print("Gene-confirmed BZLF1 records:", len(rescued_ids))
# Process each target protein
for protein, info in targets.items():

    print(f"\nProcessing {protein}...")

    # Load and verify the reference protein
    reference = SeqIO.read(BASE / info["reference"], "fasta")

    if reference.id != info["id"]:
        raise ValueError(f"Incorrect reference: {protein}")

    ref_seq = str(reference.seq).upper()
    unique = {ref_seq: reference.id}

    # Load NCBI protein descriptions
    with open(
        META / f"{protein}_inventory.csv",
        newline="",
        encoding="utf-8"
    ) as file:
        descriptions = {
            row["accession"]: row["description"]
            for row in csv.DictReader(file)
        }

    # Load downloaded protein sequences
    records = list(
        SeqIO.parse(RAW / f"{protein}.fasta", "fasta")
    )

    audit = []

    # Screen each individual record
    for record in records:

        seq = str(record.seq).upper()
        description = descriptions.get(record.id, "")

        # Check the annotated protein identity
        name_match = bool(
            re.search(info["pattern"], description, re.I)
        )

        # Check whether the sequence can enter the MSA
        if len(seq) < info["region_length"]:
            status = "TOO_SHORT"
            evidence = "LENGTH_SCREEN"

        elif set(seq) - STANDARD_AA:
            status = "REVIEW_NONSTANDARD"
            evidence = "NONSTANDARD_RESIDUES"

        elif name_match:
            status = "INCLUDED"
            evidence = "ANNOTATION_MATCH"

        elif protein == "BZLF1" and record.id in rescued_ids:
            status = "INCLUDED"
            evidence = "GENBANK_GENE_CONFIRMED"

        else:
            status = "REVIEW_IDENTITY"
            evidence = "IDENTITY_UNRESOLVED"

        alignment_id = ""

        # Add eligible sequences without duplicating identical proteins
        if status == "INCLUDED":

            if seq not in unique:
                unique[seq] = record.id

            alignment_id = unique[seq]

        # Record the decision for every accession
        audit.append({
            "protein": protein,
            "accession": record.id,
            "description": description,
            "length": len(seq),
            "status": status,
            "identity_evidence": evidence,
            "alignment_id": alignment_id
        })

    # Save the complete identity audit
    audit_file = OUT / f"{protein}_identity_audit.csv"

    with open(audit_file, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=list(audit[0].keys())
        )
        writer.writeheader()
        writer.writerows(audit)

    # Display screening results
    print("Total records:", len(records))
    print("Screening:", dict(Counter(
        row["status"] for row in audit
    )))

    provisional = sum(
        row["identity_evidence"] == "SIMILARITY_PROVISIONAL"
        for row in audit
    )

    print("Provisionally included:", provisional)
    print("Unique sequences for MSA:", len(unique))

    # Convert unique sequences to PyFAMSA format
    sequences = [
        Sequence(seq_id.encode(), seq.encode())
        for seq, seq_id in unique.items()
    ]

    # Perform multiple sequence alignment
    msa = aligner.align(sequences)

    aligned_records = [
        SeqRecord(
            Seq(row.sequence.decode()),
            id=row.id.decode(),
            description=""
        )
        for row in msa
    ]

    # Verify that the reference sequence was preserved
    aligned_reference = next(
        row for row in aligned_records
        if row.id == reference.id
    )

    if str(aligned_reference.seq).replace("-", "") != ref_seq:
        raise ValueError(f"Reference mismatch: {protein}")

    # Save the updated alignment
    output = OUT / f"{protein}_screened_aligned.fasta"

    SeqIO.write(aligned_records, output, "fasta")

    print("New alignment length:", len(aligned_records[0].seq))

print("\nIdentity screening and alignment completed!")
