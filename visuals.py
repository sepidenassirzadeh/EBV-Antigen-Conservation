
from pathlib import Path
import csv
import matplotlib.pyplot as plt

# Define input and output folders
BASE = Path(__file__).resolve().parent
TABLES = BASE / "results" / "tables"
FIGURES = BASE / "results" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

# Load conservation summary
with open(TABLES / "conservation_summary.csv",
          newline="", encoding="utf-8") as file:
    summary = list(csv.DictReader(file))

# Figure 1: Compare exact reference-peptide match percentages
proteins = [row["protein"] for row in summary]
record_rates = [
    float(row["record_match_percent"]) for row in summary
]
unique_rates = [
    float(row["unique_match_percent"]) for row in summary
]

x = list(range(len(proteins)))
width = 0.35

fig, ax = plt.subplots(figsize=(9, 5))

ax.bar([i - width / 2 for i in x],
       record_rates, width, label="All eligible records")

ax.bar([i + width / 2 for i in x],
       unique_rates, width, # Label deduplicated protein records accurately
       label="Unique protein sequences")

ax.set_xticks(x, proteins)
ax.set_ylim(0, 110)
ax.set_ylabel("Exact reference-peptide match (%)")
ax.set_title("EBV antigen-region exact-match comparison")
ax.legend()
fig.tight_layout()

fig.savefig(FIGURES / "figure1_exact_match.png", dpi=300)
plt.close(fig)

# Load position-level conservation results
with open(TABLES / "position_conservation.csv",
          newline="", encoding="utf-8") as file:
    positions = list(csv.DictReader(file))

# Figure 2: Extract one entropy value per EBNA3A position
entropy_by_position = {}

for row in positions:
    if row["protein"] == "EBNA3A":
        position = int(row["position_in_region"])
        entropy_by_position[position] = float(
            row["shannon_entropy_bits"]
        )

ordered_positions = sorted(entropy_by_position)
entropy_values = [
    entropy_by_position[p] for p in ordered_positions
]

fig, ax = plt.subplots(figsize=(9, 5))

ax.bar(ordered_positions, entropy_values)
ax.set_xticks(ordered_positions)
ax.set_xlabel("Position within EBNA3A antigen region")
ax.set_ylabel("Shannon entropy (bits)")
ax.set_title("EBNA3A amino-acid diversity by position")

fig.tight_layout()
fig.savefig(FIGURES / "figure2_ebna3a_entropy.png", dpi=300)
plt.close(fig)

print("Both figures saved successfully!")
print("Output folder:", FIGURES)
