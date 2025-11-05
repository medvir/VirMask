#!/usr/bin/env python3
"""
Contaminant detection and viral database masking.

Steps:
1. Run Extract_overlaps.py to find sequences to mask.
2. Run Mask_contaminants.py to mask DB based on 
recurrent contaminants and homologous sequences.
"""

# Import standard libraries
import argparse
import os
import subprocess
import sys

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Detection and masking of recurrent contaminants within the viral database."
)
parser.add_argument(
    "--mNGS_runs",
    type=str, 
    default="./mNGS_to_analyse.txt",
    help="Path to a .txt file listing all mNGS runs to consider (one run per row)")
parser.add_argument(
    "--virmet_db",
    type=str, 
    default="/data/virmet_databases_update/",
    help="Path to the VirMet database (default: /data/virmet_databases_update/)")
parser.add_argument(
    "--across_sample_freq",
    type=float,
    default=0.2,
    help="Across-sample frequency threshold (default: 0.2 = 20%%)")
parser.add_argument(
    "--within_acc_dominance",
    type=float,
    default=0.2,
    help="Within-accession dominance threshold (default: 0.2 = 20%%)")
parser.add_argument(
    "--min_fragment_size",
    type=int,
    default=21,
    help="Minimum size of the fragment to be masked (default: 21 bp)")
parser.add_argument(
    "--max_fragment_size",
    type=int,
    default=499,
    help="Maximum size of the fragment to be masked (default: 499 bp)")
parser.add_argument(
    "--max_genome_perc",
    type=float,
    default=0.05,
    help="Maximum percentage of genome to mask (default: 0.05 = 5%%)")
parser.add_argument(
    "--controls",
    type=str,
    default="Tunavirus T1, phage MS2",
    help=(
        "Pattern to identify viruses that should be prevented from masking. "
        "If more than one is provided, add them within quotes and separated by commas. "
        "Example: 'Tunavirus T1, phage MS2'. "
        "This avoids masking your internal controls, which appear in all samples."
        )
)
parser.add_argument(
    "--threads",
    type=int,
    default=24,
    help="Number of threads to use (default: 24)"
)
args = parser.parse_args()

# Define script directory
script_dir = os.path.dirname(os.path.abspath(__file__))

# Ensure that stdout and stderr are visible while the code is running
def run_command_live(cmd, log_path):
    """Run a command and show output to both console and log file."""
    with open(log_path, "w") as log:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        # Stream output line by line
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
            log.flush()
        process.wait()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, cmd)

# Run a script
def run_script(script_name, extra_args=None):
    script_path = os.path.join(script_dir, script_name)
    if not os.path.isfile(script_path):
        print(f"Error: Script {script_path} not found.")
        sys.exit(1)

    cmd = ["python3", script_path]
    if extra_args:
        cmd.extend(extra_args)

    print(f"\nRunning: {' '.join(cmd)}")
    log_path = f"{os.path.splitext(script_name)[0]}.log"
    try:
        run_command_live(cmd, log_path)
    except subprocess.CalledProcessError as e:
        print(f"Error: {script_name} failed (exit code {e.returncode}). See log: {log_path}")
        sys.exit(e.returncode)

# Extract contaminants
extract_args = [
    "--mNGS_runs", str(args.mNGS_runs),
    "--virmet_db", str(args.virmet_db),
    "--across_sample_freq", str(args.across_sample_freq),
    "--within_acc_dominance", str(args.within_acc_dominance),
    "--controls", str(args.controls),
]
run_script("Extract_contaminants.py", extract_args)

# Mask contaminants
mask_args = [
    "--virmet_db", str(args.virmet_db),
    "--min_fragment_size", str(args.min_fragment_size),
    "--max_fragment_size", str(args.max_fragment_size),
    "--max_genome_perc", str(args.max_genome_perc),
    "--controls", str(args.controls),
    "--threads", str(args.threads),
]
run_script("Mask_contaminants.py", mask_args)

print("\nMasking of recurrent contaminants completed successfully!")