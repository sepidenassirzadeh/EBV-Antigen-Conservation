from Bio import SeqIO
from Bio.Align import PairwiseAligner
from pathlib import Path

# Locate downloaded protein sequences
BASE_DIR = Path(__file__).resolve().parent
fasta_path = BASE_DIR / "data" / "raw" / "EBNA3A.fasta"

records = list(SeqIO.parse(fasta_path, "fasta"))


# Load the previously validated EBNA3A reference
reference_path = BASE_DIR / "EBNA-3A.fasta"
reference = SeqIO.read(reference_path, "fasta")
# Select a candidate with a different amino acid sequence
candidate = next(
    (
        record for record in records
        if str(record.seq) != str(reference.seq)
    ),
    None
)

# Ensure a different candidate was found
if candidate is None:
    raise ValueError("No different EBNA3A sequence found")

# Display the selected candidate
print("Candidate ID:", candidate.id)
print("Candidate length:", len(candidate.seq))

# Verify the reference accession before alignment
if reference.id != "YP_401669.1":
    raise ValueError("Incorrect EBNA3A reference")

print("Reference loaded:", reference.id)
print("Reference length:", len(reference.seq))
# Configure local protein sequence alignment
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5

# Align the candidate sequence to the reference
alignment = aligner.align(
    reference.seq,
    candidate.seq
)[0]

print("Reference:", reference.id)
print("Candidate:", candidate.id)
print("Reference length:", len(reference.seq))
print("Candidate length:", len(candidate.seq))
print("Alignment score:", alignment.score)

# Map reference positions 325-333 to the candidate
region_positions = range(324, 333)

position_map = {
    int(ref_pos): int(query_pos)
    for ref_pos, query_pos in zip(*alignment.indices)
    if ref_pos >= 0
}

mapped_positions = [
    position_map.get(pos, -1)
    for pos in region_positions
]

# Report whether the entire antigen region is covered
if -1 in mapped_positions:
    print("Region status: INCOMPLETE")
else:
    region = "".join(
        str(candidate.seq[pos])
        for pos in mapped_positions
    )
    print("Reference region: FLRGRAYGI")
    print("Candidate region:", region)
    print("Region length:", len(region))
    print("Region status: MAPPED")

# Check how many reference-region positions are aligned
covered_positions = sum(pos != -1 for pos in mapped_positions)

# Report the coverage of the target antigen region
print("Covered positions:", covered_positions, "/ 9")
print("Region coverage:", round(covered_positions / 9 * 100, 2), "%")

# Identify unmapped positions in the reference antigen region
for ref_pos, candidate_pos in zip(region_positions, mapped_positions):
    if candidate_pos == -1:
        print("Unmapped reference position:", ref_pos + 1)

# Display the aligned sequences around the target region
print("\nAlignment:")
print(alignment)

# Show the reference and candidate alignment boundaries
print("Aligned blocks:", alignment.aligned)

# Show the reference position 325 and its corresponding mapping
print("Reference amino acid 325:", reference.seq[324])
print("Mapped candidate position:", position_map.get(324))