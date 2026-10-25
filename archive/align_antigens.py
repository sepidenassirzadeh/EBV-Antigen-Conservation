
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfamsa import Aligner, Sequence
from pathlib import Path
import csv

# Set up input and output folders
BASE_DIR = Path(__file__).resolve().parent
RAW_DIR = BASE_DIR / "data" / "raw"
OUT_DIR = BASE_DIR / "results" / "alignments"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Define validated reference sequences
targets = {
    "EBNA1": ("EBNA1.fasta", "YP_401677.1"),
    "EBNA3A": ("EBNA-3A.fasta", "YP_401669.1"),
    "BZLF1": ("BZLF1.fasta", "YP_401673.1")
}

# Use standard amino acids for this alignment run
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Initialize the multiple sequence aligner
aligner = Aligner(threads=2)

for protein, (ref_file, ref_id) in targets.items():

    print(f"\nProcessing {protein}...")

    # Load and validate the reference sequence
    reference = SeqIO.read(BASE_DIR / ref_file, "fasta")
    ref_seq = str(reference.seq).upper()

    if reference.id != ref_id:
        raise ValueError(f"Incorrect reference for {protein}")

    if set(ref_seq) - STANDARD_AA:
        raise ValueError(f"Non-standard reference: {protein}")

    # Load downloaded sequences
    records = list(
        SeqIO.parse(RAW_DIR / f"{protein}.fasta", "fasta")
    )

    # Keep one alignment copy of each distinct sequence
    unique_sequences = {ref_seq: ref_id}
    accession_map = []
    excluded_count = 0

    for record in records:
        seq = str(record.seq).upper()

        # Identify residues unsupported in this alignment run
        nonstandard = sorted(set(seq) - STANDARD_AA)

        if nonstandard:
            excluded_count += 1

            accession_map.append({
                "accession": record.id,
                "alignment_id": "",
                "status": "EXCLUDED_FROM_MSA",
                "reason": "Non-standard residues: " + ",".join(nonstandard)
            })
            continue

        # Assign a representative ID to identical sequences
        if seq not in unique_sequences:
            unique_sequences[seq] = record.id

        accession_map.append({
            "accession": record.id,
            "alignment_id": unique_sequences[seq],
            "status": "INCLUDED_IN_MSA",
            "reason": ""
        })

    print("Total records:", len(records))
    print("Unique alignment sequences:", len(unique_sequences))
    print("Excluded from this MSA:", excluded_count)

    # Convert unique sequences to PyFAMSA format
    sequences = [
        Sequence(seq_id.encode(), seq.encode())
        for seq, seq_id in unique_sequences.items()
    ]

    # Perform multiple sequence alignment
    msa = aligner.align(sequences)

    # Convert alignment results to FASTA records
    aligned_records = [
        SeqRecord(
            Seq(row.sequence.decode()),
            id=row.id.decode(),
            description=""
        )
        for row in msa
    ]

    # Verify the reference sequence after alignment
    aligned_reference = next(
        r for r in aligned_records if r.id == ref_id
    )

    if str(aligned_reference.seq).replace("-", "") != ref_seq:
        raise ValueError(f"Reference mismatch: {protein}")

    # Save the aligned protein sequences
    output_fasta = OUT_DIR / f"{protein}_aligned.fasta"
    SeqIO.write(aligned_records, output_fasta, "fasta")

    # Save inclusion, exclusion, and accession mapping
    output_csv = OUT_DIR / f"{protein}_accession_map.csv"

    with open(output_csv, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=[
                "accession", "alignment_id", "status", "reason"
            ]
        )
        writer.writeheader()
        writer.writerows(accession_map)

    print("Alignment saved:", output_fasta.name)
    print("Alignment length:", len(aligned_records[0].seq))

print("\nAll three alignments completed!")
