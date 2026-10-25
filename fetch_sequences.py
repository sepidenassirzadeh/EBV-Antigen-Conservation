
from Bio import Entrez, SeqIO
from pathlib import Path
from datetime import date
from http.client import RemoteDisconnected
from urllib.error import URLError, HTTPError
import csv
import time
import os


# Set up NCBI access and project directories
Entrez.email = os.environ.get("NCBI_EMAIL")

if not Entrez.email:
    raise ValueError("Please set the NCBI_EMAIL environment variable.")

# Define the target proteins
proteins = ["EBNA1", "EBNA3A", "BZLF1"]

# Search NCBI and retrieve all matching protein IDs
for protein in proteins:

    query = f'"Human gammaherpesvirus 4"[Organism] AND {protein}[All Fields]'

    with Entrez.esearch(
        db="protein",
        term=query,
        retmax=10000
    ) as handle:
        results = Entrez.read(handle)

    ids = results["IdList"]
    total = int(results["Count"])

    print(f"\n{protein}: {total} search results")

    if len(ids) != total:
        raise ValueError(f"Incomplete ID retrieval for {protein}")

    fasta_path = RAW_DIR / f"{protein}.fasta"
    csv_path = META_DIR / f"{protein}_inventory.csv"

    # Load previously saved accessions to resume interrupted downloads
    saved_ids = set()

    if csv_path.exists():
        with open(csv_path, newline="", encoding="utf-8") as file:
            saved_ids = {
                row["accession"]
                for row in csv.DictReader(file)
            }

    # Fetch records in batches with retries
    for start in range(0, len(ids), 10):

        batch = ids[start:start + 10]

        for attempt in range(5):
            try:
                with Entrez.efetch(
                    db="protein",
                    id=",".join(batch),
                    rettype="gb",
                    retmode="text"
                ) as handle:
                    batch_records = list(
                        SeqIO.parse(handle, "genbank")
                    )

                if len(batch_records) != len(batch):
                    raise ValueError("Incomplete batch response")

                break

            except (RemoteDisconnected, URLError, HTTPError,
                    ConnectionResetError, TimeoutError, ValueError) as error:

                print(f"Retry {attempt + 1}/5: {error}")

                if attempt == 4:
                    raise

                time.sleep(3 * (attempt + 1))

        # Keep only records that have not already been saved
        new_records = [
            record for record in batch_records
            if record.id not in saved_ids
        ]

        if new_records:

            # Append validated records to FASTA
            with open(fasta_path, "a") as fasta_file:
                SeqIO.write(new_records, fasta_file, "fasta")

            # Append metadata to the inventory CSV
            file_exists = csv_path.exists() and csv_path.stat().st_size > 0

            with open(csv_path, "a", newline="", encoding="utf-8") as file:
                writer = csv.DictWriter(
                    file,
                    fieldnames=[
                        "target_protein",
                        "accession",
                        "description",
                        "length",
                        "partial_label",
                        "retrieval_date"
                    ]
                )

                if not file_exists:
                    writer.writeheader()

                for record in new_records:
                    writer.writerow({
                        "target_protein": protein,
                        "accession": record.id,
                        "description": record.description,
                        "length": len(record.seq),
                        "partial_label": "partial" in record.description.lower(),
                        "retrieval_date": date.today().isoformat()
                    })

                    saved_ids.add(record.id)

        print(f"{protein}: saved {len(saved_ids)} records")

print("\nSequence retrieval completed!")

