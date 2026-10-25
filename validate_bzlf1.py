
from Bio import Entrez, SeqIO
from pathlib import Path
from collections import Counter
import csv
import os
import time

# Define project paths
BASE = Path(__file__).resolve().parent
INPUT = BASE / "results/tables/BZLF1_identity_review.csv"
OUTPUT = BASE / "results/tables/BZLF1_gene_validation.csv"

Entrez.email = os.environ.get("NCBI_EMAIL")

if not Entrez.email:
    raise ValueError("Please set the NCBI_EMAIL environment variable.")

# Load provisionally identified proteins
with open(INPUT, newline="", encoding="utf-8") as file:
    accessions = [
        row["accession"]
        for row in csv.DictReader(file)
        if row["screening_status"] == "LIKELY_BZLF1"
    ]

results = []

# Retrieve GenBank annotations in batches
for start in range(0, len(accessions), 20):
    batch = accessions[start:start + 20]

    with Entrez.efetch(
        db="protein",
        id=",".join(batch),
        rettype="gp",
        retmode="text"
    ) as handle:
        records = list(SeqIO.parse(handle, "genbank"))

    for record in records:

        # Extract annotated gene names
        genes = set()

        for feature in record.features:
            genes.update(feature.qualifiers.get("gene", []))

        genes = {gene.strip().upper() for gene in genes}

        # Classify gene-annotation evidence
        if "BZLF1" in genes:
            status = "GENE_CONFIRMED"
        elif genes:
            status = "OTHER_GENE_ANNOTATED"
        else:
            status = "GENE_UNRESOLVED"

        results.append({
            "accession": record.id,
            "description": record.description,
            "annotated_genes": "; ".join(sorted(genes)),
            "validation_status": status
        })

    print("Retrieved:", min(start + 20, len(accessions)))
    time.sleep(0.4)

# Check that all requested records were returned
returned = {row["accession"] for row in results}
missing = sorted(set(accessions) - returned)

if missing:
    print("Missing accessions:", missing)

# Save gene-validation evidence
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with open(OUTPUT, "w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(file, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

# Display validation summary
print("\nGene validation results:")
print(dict(Counter(
    row["validation_status"] for row in results
)))

print("\nSaved:", OUTPUT)
