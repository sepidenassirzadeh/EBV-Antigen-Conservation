
from pathlib import Path
import csv
import matplotlib.pyplot as plt

# Define input and output directories
BASE = Path(__file__).resolve().parent
TABLES = BASE / "results" / "tables"
FIGURES = BASE / "results" / "figures"

FIGURES.mkdir(parents=True, exist_ok=True)

# Load conservation summary
with open(
    TABLES / "conservation_summary.csv",
    newline="",
    encoding="utf-8"
) as file:
    summary = list(csv.DictReader(file))

# Verify that all three target proteins are present
expected_proteins = {"EBNA1", "EBNA3A", "BZLF1"}

if (
    len(summary) != 3
    or {row["protein"] for row in summary} != expected_proteins
):
    raise ValueError(
        "Conservation summary must contain exactly three target proteins."
    )

# Extract exact-match percentages
proteins = [row["protein"] for row in summary]

record_rates = [
    float(row["record_match_percent"])
    for row in summary
]

unique_rates = [
    float(row["unique_match_percent"])
    for row in summary
]

# Verify that percentage values are valid
if any(
    not 0 <= rate <= 100
    for rate in record_rates + unique_rates
):
    raise ValueError("Invalid exact-match percentage detected.")

# Figure 1: Compare record-level and unique-sequence match rates
x = list(range(len(proteins)))
width = 0.35

fig, ax = plt.subplots(figsize=(9, 5))

# Plot exact-match percentages for eligible NCBI records
ax.bar(
    [i - width / 2 for i in x],
    record_rates,
    width,
    label="All eligible records"
)

# Plot exact-match percentages for unique protein sequences
ax.bar(
    [i + width / 2 for i in x],
    unique_rates,
    width,
    label="Unique protein sequences"
)

# Format the first figure
ax.set_xticks(x, proteins)
ax.set_ylim(0, 110)
ax.set_ylabel("Exact reference-peptide match (%)")
ax.set_title("EBV antigen-region exact-match comparison")
ax.legend()

fig.tight_layout()

# Save the first figure
fig.savefig(
    FIGURES / "figure1_exact_match.png",
    dpi=300
)

plt.close(fig)

# Load position-level conservation results
with open(
    TABLES / "position_conservation.csv",
    newline="",
    encoding="utf-8"
) as file:
    positions = list(csv.DictReader(file))

# Extract Shannon entropy values for EBNA3A
entropy_by_position = {}

for row in positions:

    if row["protein"] == "EBNA3A":
        position = int(row["position_in_region"])
        entropy = float(row["shannon_entropy_bits"])

        # Verify consistent entropy values within each position
        if (
            position in entropy_by_position
            and entropy_by_position[position] != entropy
        ):
            raise ValueError(
                f"Inconsistent entropy at EBNA3A position {position}"
            )

        entropy_by_position[position] = entropy

# Verify all nine EBNA3A peptide positions are present
if set(entropy_by_position) != set(range(1, 10)):
    raise ValueError(
        "Missing or unexpected EBNA3A antigen positions."
    )

# Verify entropy values are nonnegative
if any(value < 0 for value in entropy_by_position.values()):
    raise ValueError("Negative Shannon entropy detected.")

# Prepare positions and entropy values for plotting
ordered_positions = sorted(entropy_by_position)

entropy_values = [
    entropy_by_position[position]
    for position in ordered_positions
]

# Figure 2: Plot amino-acid variability across EBNA3A
fig, ax = plt.subplots(figsize=(9, 5))

ax.bar(ordered_positions, entropy_values)

# Format the second figure
ax.set_xticks(ordered_positions)
ax.set_xlabel("Position within EBNA3A antigen region")
ax.set_ylabel("Shannon entropy (bits)")
ax.set_title("EBNA3A amino-acid diversity by position")

fig.tight_layout()

# Save the second figure
fig.savefig(
    FIGURES / "figure2_ebna3a_entropy.png",
    dpi=300
)

plt.close(fig)

print("Both figures saved successfully!")
print("Output folder:", FIGURES)
