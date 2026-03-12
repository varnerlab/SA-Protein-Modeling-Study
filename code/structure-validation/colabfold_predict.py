#!/usr/bin/env python3
"""
ColabFold AlphaFold2 structure prediction for SA protein study.

Run this script on Google Colab (with GPU runtime) to predict structures
for SA-generated and stored protein sequences using AlphaFold2.

Usage on Colab:
  1. Upload the colabfold_input/ directory to Colab
  2. Run this script (or paste into cells)
  3. Download the colabfold_output/ directory when done

The script processes one family at a time with checkpointing,
so it can be restarted across multiple Colab sessions.
"""

import os
import subprocess
import sys
import json
import glob

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

INPUT_DIR = "colabfold_input"
OUTPUT_DIR = "colabfold_output"
NUM_MODELS = 1          # only need 1 model for pLDDT/TM comparison
NUM_RECYCLES = 3        # default recycling
USE_AMBER = False       # skip relaxation to save time

FAMILIES = [
    "PF00397",   # WW (shortest, ~23 residues — fastest)
    "PF00096",   # zf-C2H2 (~23 residues)
    "PF00018",   # SH3 (~60 residues)
    "PF00014",   # Kunitz (~58 residues)
    "PF00076",   # RRM (~70 residues)
    "PF00595",   # PDZ (~90 residues)
    "PF00069",   # Pkinase (~262 residues — longest, do last)
]

CATEGORIES = ["sa_generation", "stored"]

# ─────────────────────────────────────────────────────────────────────────────
# Install ColabFold (run once per Colab session)
# ─────────────────────────────────────────────────────────────────────────────

def install_colabfold():
    """Install ColabFold if not already present."""
    try:
        import colabfold
        print("ColabFold already installed.")
        return
    except ImportError:
        pass

    print("Installing ColabFold...")
    subprocess.check_call([
        sys.executable, "-m", "pip", "install",
        "colabfold[alphafold] @ git+https://github.com/sokrypton/ColabFold.git",
        "--quiet"
    ])
    print("ColabFold installed successfully.")

# ─────────────────────────────────────────────────────────────────────────────
# Run predictions
# ─────────────────────────────────────────────────────────────────────────────

def count_completed(output_subdir):
    """Count how many sequences already have predictions."""
    if not os.path.isdir(output_subdir):
        return 0
    # ColabFold outputs *_unrelaxed_rank_001_*.pdb for each sequence
    pdbs = glob.glob(os.path.join(output_subdir, "*_unrelaxed_rank_001_*.pdb"))
    return len(pdbs)

def count_input(fasta_path):
    """Count sequences in a FASTA file."""
    if not os.path.isfile(fasta_path):
        return 0
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count

def run_colabfold_batch(fasta_path, output_subdir):
    """Run colabfold_batch on a FASTA file."""
    os.makedirs(output_subdir, exist_ok=True)

    n_input = count_input(fasta_path)
    n_done = count_completed(output_subdir)

    if n_done >= n_input and n_input > 0:
        print(f"  Already complete: {n_done}/{n_input} predictions")
        return

    print(f"  Running ColabFold: {n_input} sequences ({n_done} already done)")

    cmd = [
        "colabfold_batch",
        fasta_path,
        output_subdir,
        "--num-models", str(NUM_MODELS),
        "--num-recycle", str(NUM_RECYCLES),
        "--model-type", "alphafold2_ptm",
    ]

    if not USE_AMBER:
        cmd.append("--amber")  # this flag is for enabling amber; omitting enables unrelaxed only
        # Actually, ColabFold runs unrelaxed by default; --amber adds relaxation
        # So we do NOT add --amber to skip relaxation
        cmd = [c for c in cmd if c != "--amber"]

    print(f"  Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def main():
    # Install ColabFold
    install_colabfold()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    total_done = 0
    total_seqs = 0

    for fam_id in FAMILIES:
        print(f"\n{'='*60}")
        print(f"  {fam_id}")
        print(f"{'='*60}")

        for category in CATEGORIES:
            fasta_path = os.path.join(INPUT_DIR, fam_id, f"{category}.fasta")
            output_subdir = os.path.join(OUTPUT_DIR, fam_id, category)

            if not os.path.isfile(fasta_path):
                print(f"  {category}: FASTA not found, skipping")
                continue

            n_input = count_input(fasta_path)
            total_seqs += n_input

            print(f"\n  {category} ({n_input} sequences):")
            run_colabfold_batch(fasta_path, output_subdir)

            n_done = count_completed(output_subdir)
            total_done += n_done
            print(f"  {category}: {n_done}/{n_input} complete")

    print(f"\n{'='*60}")
    print(f"  TOTAL: {total_done}/{total_seqs} predictions complete")
    print(f"  Output: {OUTPUT_DIR}/")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
