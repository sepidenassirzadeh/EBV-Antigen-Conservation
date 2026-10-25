from Bio import Entrez, SeqIO
from pathlib import Path
import csv

# Set the email address required for NCBI access
Entrez.email = os.environ.get("NCBI_EMAIL")

if not Entrez.email:
    raise ValueError("Please set the NCBI_EMAIL environment variable.")

# Search NCBI for EBNA3A protein records
query = '"Human gammaherpesvirus 4"[Organism] AND EBNA3A[All Fields]'

with Entrez.esearch(db="protein", term=query, retmax=200) as handle:
    results = Entrez.read(handle)

ids = results["IdList"]

print("Total results:", results["Count"])
print("IDs retrieved:", len(ids))

# Create a folder for metadata files
Path("metadata").mkdir(exist_ok=True)

# Create an empty list to store protein information
records_info = []

# Fetch protein records in batches of 25
for start in range(0, len(ids), 25):
    batch = ids[start:start + 25]

    with Entrez.efetch(
        db="protein",
        id=",".join(batch),
        rettype="gb",
        retmode="text"
    ) as handle:

        # Extract metadata from each protein record
        for record in SeqIO.parse(handle, "genbank"):
            records_info.append({
                "accession": record.id,
                "description": record.description,
                "length": len(record.seq),
                "partial_label": "partial" in record.description.lower()
            })

# Save all retrieved records in a CSV inventory
with open(
    "metadata/ebna3a_inventory.csv",
    "w",
    newline="",
    encoding="utf-8"
) as file:

    writer = csv.DictWriter(
        file,
        fieldnames=["accession", "description", "length", "partial_label"]
    )

    writer.writeheader()
    writer.writerows(records_info)

print("Records saved:", len(records_info))
print("Inventory saved successfully!")
