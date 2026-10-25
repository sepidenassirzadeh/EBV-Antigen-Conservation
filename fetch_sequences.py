from Bio import Entrez, SeqIO
from pathlib import Path
from datetime import date
from http.client import RemoteDisconnected
from urllib.error import URLError, HTTPError
import csv
import os
import time

# Configure NCBI access
email = os.environ.get("NCBI_EMAIL")

if not email:
    raise ValueError("Please set the NCBI_EMAIL environment variable.")

Entrez.email = email

# Define project directories
BASE = Path(__file__).resolve().parent
RAW_DIR = BASE / "data" / "raw"
META_DIR = BASE / "metadata"

RAW_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)

# Define target proteins
PROTEINS = ["EBNA1", "EBNA3A", "BZLF1"]

# Retrieve NCBI records with retries
def fetch_batch(ids):
    for attempt in range(5):
        try:
            with Entrez.efetch(
                db="protein",
                id=",".join(ids),
                rettype="gb",
                retmode="text"
            ) as handle:
                records = list(SeqIO.parse(handle, "genbank"))

            if len(records) != len(ids):
                raise ValueError("Incomplete batch response")

            return records

        except (
            RemoteDisconnected,
            URLError,
            HTTPError,
            ConnectionResetError,
            TimeoutError,
            ValueError
        ) as error:
            print(f"Retry {attempt + 1}/5: {error}")

            if attempt == 4:
                raise

            time.sleep(3 * (attempt + 1))

# Process each target protein
for protein in PROTEINS:

    fasta_path = RAW_DIR / f"{protein}.fasta"
    csv_path = META_DIR / f"{protein}_inventory.csv"

    # Protect existing input datasets
    if fasta_path.exists() or csv_path.exists():
        print(f"{protein}: existing dataset found; skipping download.")
        continue

    # Search for candidate protein records
    query = (
        f'"Human gammaherpesvirus 4"[Organism] '
        f'AND {protein}[All Fields]'
    )

    with Entrez.esearch(
        db="protein",
        term=query,
        retmax=10000
    ) as handle:
        results = Entrez.read(handle)

    ids = list(results["IdList"])
    total = int(results["Count"])

    if len(ids) != total:
        raise ValueError(f"Incomplete ID retrieval: {protein}")

    print(f"\n{protein}: {total} candidate records")

    # Download all records before saving output files
    records = []

    for start in range(0, len(ids), 10):
        batch = ids[start:start + 10]
        records.extend(fetch_batch(batch))

    # Verify downloaded accession uniqueness
    accessions = [record.id for record in records]

    if len(accessions) != len(set(accessions)):
        raise ValueError(f"Duplicate accessions: {protein}")

    if len(records) != total:
        raise ValueError(f"Record count mismatch: {protein}")

    # Prepare metadata for every downloaded record
    metadata = []

    for record in records:
        metadata.append({
            "target_protein": protein,
            "accession": record.id,
            "description": record.description,
            "length": len(record.seq),
            "partial_label": (
                "partial" in record.description.lower()
            ),
            "retrieval_date": date.today().isoformat()
        })

    # Write downloaded sequences
    SeqIO.write(records, fasta_path, "fasta")

    # Write accession metadata
    with open(
        csv_path,
        "w",
        newline="",
        encoding="utf-8"
    ) as file:
        writer = csv.DictWriter(
            file,
            fieldnames=list(metadata[0].keys())
        )
        writer.writeheader()
        writer.writerows(metadata)

    print(f"{protein}: saved {len(records)} records")

print("\nSequence retrieval completed!")
