#!/usr/bin/env python3
"""
Contaminant detection and viral database masking.

Steps:
1. Run Mask_human_regions.py to detect and mask human regions.
2. Run Summarise_human_regions.py to provide a summary of
the human regions that have been masked.
"""

# Import standard libraries
import argparse
import os
import subprocess
import sys

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Detection and masking of human sequences within the viral database."
)
parser.add_argument(
    "--virmet_db",
    type=str, 
    default="/data/virmet_databases_update/",
    help="Path to the VirMet database (default: /data/virmet_databases_update/)"
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

# Run a regular Python script
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

# Run via Snakemake
def run_snakemake(script_name, threads, config):
    script_path = os.path.join(script_dir, script_name)
    if not os.path.isfile(script_path):
        print(f"Error: Script {script_path} not found.")
        sys.exit(1)

    cmd = [
        "snakemake",
        "-s", script_path,
        "--cores", str(threads),
        "--latency-wait", str(300),
        "--config", f"virmet_db={config}"
    ]

    print(f"\nRunning Snakemake: {' '.join(cmd)}")
    log_path = f"{os.path.splitext(script_name)[0]}.log"
    try:
        run_command_live(cmd, log_path)
    except subprocess.CalledProcessError as e:
        print(f"Error: Snakemake failed (exit code {e.returncode}). See log: {log_path}")
        sys.exit(e.returncode)

# Mask human contaminants
run_snakemake(
    "Mask_human_regions.py",
    threads=args.threads,
    config=args.virmet_db
)

# Summarise human contaminants
summary_args = ["--virmet_db", str(args.virmet_db)]
run_script("Summarise_human_regions.py", summary_args)

print("\nMasking of human sequences completed successfully!")