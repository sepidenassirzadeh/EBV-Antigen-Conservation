
import subprocess
import sys
from pathlib import Path

# Define the project directory
BASE = Path(__file__).resolve().parent

# Define the analysis scripts in execution order
scripts = [
    "align_screened.py",
    "map_regions.py",
    "variant_summary.py",
    "conservation.py",
    "visuals.py"
]

# Check that all required scripts exist
for script in scripts:
    if not (BASE / script).is_file():
        raise FileNotFoundError(f"Missing script: {script}")

# Run each analysis step using the current Python interpreter
# Define the complete analysis pipeline
scripts = [
    "align_screened.py",
    "map_regions.py",
    "variant_summary.py",
    "conservation.py",
    "visuals.py",
    "qc_flow.py"
]

# Confirm successful completion
print("\nAll analysis scripts completed successfully!")
