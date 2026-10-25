import subprocess
import sys
from pathlib import Path

# Define the project directory
BASE = Path(__file__).resolve().parent

# Define scripts in the correct execution order
SCRIPTS = [
    "align_screened.py",
    "map_regions.py",
    "variant_summary.py",
    "conservation.py",
    "visuals.py",
    "qc_flow.py"
]

# Check that all analysis scripts exist
for script in SCRIPTS:
    script_path = BASE / script

    if not script_path.is_file():
        raise FileNotFoundError(
            f"Required script not found: {script}"
        )

# Execute every analysis step
for script in SCRIPTS:
    print(f"\n{'=' * 45}", flush=True)
    print(f"Running: {script}", flush=True)
    print(f"{'=' * 45}", flush=True)

    subprocess.run(
        [sys.executable, str(BASE / script)],
        cwd=BASE,
        check=True
    )

# Report success only after all scripts finish
print("\nAll analysis scripts completed successfully!")
