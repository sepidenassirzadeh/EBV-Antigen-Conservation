from Bio import Entrez, SeqIO
from pathlib import Path
import csv
import os

# Configure NCBI access
email = os.environ.get("NCBI_EMAIL")

if not email:
    raise ValueError("Please set the NCBI_EMAIL environment variable.")

Entrez.email = email

# Define project paths
BASE = Path(__file__).resolve().parent
OUTPUT_DIR = BASE / "metadata"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Save the exploratory inventory separately
OUTPUT = OUTPUT_DIR / "EBNA3A_inventory_review.csv"

# Prevent accidental overwriting
if OUTPUT.exists():
    raise FileExistsError(
        f"Review inventory already exists: {OUTPUT}"
    )

# Search NCBI for EBNA3A candidate records
query = (
    '"Human gammaherpesvirus 4"[Organism] '
    'AND EBNA3A[All Fields]'
)

with Entrez.esearch(
    db="protein",
    term=query,
    retmax=10000
) as handle:
    results = Entrez.read(handle)

ids = list(results["IdList"])
total = int(results["Count"])

# Check that all search IDs were retrieved
if len(ids) != total:
    raise ValueError(
        f"Incomplete search results: {len(ids)}/{total}"
    )

print("Total results:", total)

# Collect protein metadata
records_info = []

for start in range(0, len(ids), 25):
    batch = ids[start:start + 25]

    with Entrez.efetch(
        db="protein",
        id=",".join(batch),
        rettype="gb",
        retmode="text"
    ) as handle:

        for record in SeqIO.parse(handle, "genbank"):
            records_info.append({
                "accession": record.id,
                "description": record.description,
                "length": len(record.seq),
                "partial_label": (
                    "partial" in record.description.lower()
                )
            })

# Verify the number of retrieved records
if len(records_info) != total:
    raise ValueError(
        f"Retrieved {len(records_info)} of {total} records"
    )

# Save the review inventory
with open(
    OUTPUT,
    "w",
    newline="",
    encoding="utf-8"
) as file:
    writer = csv.DictWriter(
        file,
        fieldnames=[
            "accession",
            "description",
            "length",
            "partial_label"
        ]
    )

    writer.writeheader()
    writer.writerows(records_info)

print("Records saved:", len(records_info))
print("Review inventory saved:", OUTPUT)
