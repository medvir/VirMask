#!/usr/bin/env python3
"""
This script identifies and masks contaminant-prone regions in a viral reference database.
1. Reads candidate contaminant intervals from "Recurrent_contaminant_regions.tsv".
2. Extracts the corresponding subsequences from the reference viral FASTA.
3. Filters out very short (e.g., <21 bp) or long (e.g., >499 bp) sequences.
4. Performs a BLASTN search of these subsequences against the full viral database.
5. Excludes subsequences that align to internal control viruses (e.g., MS2, T1).
6. Merges overlapping BLAST hits for each viral genome into consolidated intervals.
7. Masks the reference sequences at these intervals by replacing nucleotides with 'N',
   limiting masking to ≤x% (e.g., 5%) of the genome length (except for genomes <1 kb).
8. Saves:
   - The updated masking table ("Recurrent_contaminant_regions.tsv").
   - The homologous regions ("Sequence_similarity_search.tsv").
   - The masked viral reference FASTA ("masked_viral_database.fasta").
   - The extracted masked fragments ("Homologous_regions_masked.fasta") for independent validation.
"""

# Import standard libraries
import argparse
import copy
import os
import re
import subprocess

# Import third-party libraries
import pandas as pd

# Biopython imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Mask viral database based on recurrent contaminant regions."
)
parser.add_argument(
    "--virmet_db",
    type=str, 
    default="/data/virmet_databases_update/",
    help="Path to the VirMet database (default: /data/virmet_databases_update/)"
)
parser.add_argument(
    "--min_fragment_size",
    type=int,
    default=21,
    help="Minimum size of the fragment to be masked (default: 21 bp)"
)
parser.add_argument(
    "--max_fragment_size",
    type=int,
    default=499,
    help="Maximum size of the fragment to be masked (default: 499 bp)"
)
parser.add_argument(
    "--max_genome_perc",
    type=float,
    default=0.05,
    help="Maximum percentage of genome to mask (default: 0.05 = 5%%)"
)
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

# Define function to get clean ACC
def clean_accession(acc):
    match = re.search(r'\|?([A-Z0-9_]+\.\d+)\|?$', acc)
    return match.group(1) if match else acc

# Define function to merge overlapping intervals for a single group
def merge_intervals(group):
    # Normalize intervals: always (min, max)
    intervals = [
        (min(row.sstart, row.send), max(row.sstart, row.send), 
        row.qseqid) for idx, row in group.iterrows()]
    
    # Sort intervals by start
    intervals.sort(key=lambda x: x[0])
    
    # Define first interval
    merged = []
    current_start = intervals[0][0]
    current_end = intervals[0][1]
    current_qseqids =  [intervals[0][2]]
    
    # Add all regions contained within an interval
    for start, end, qseqid in intervals[1:]:
        if start <= current_end:
            if end > current_end:
                current_end = end
            current_qseqids.append(qseqid)
        else:
            merged.append((current_start, current_end, ",".join(current_qseqids)))
            current_start, current_end, current_qseqids = start, end, [qseqid]
    
    # Add last interval
    merged.append((current_start, current_end, ",".join(current_qseqids)))
    
    merged_df = pd.DataFrame(merged, columns=['sstart', 'send', 'qseqid'])
    
    # Add the other columns (except pident, qcovs)
    merged_df['sseqid'] = group.iloc[0].sseqid
    merged_df['stitle'] = group.iloc[0].stitle
    merged_df['sscinames'] = group.iloc[0].sscinames
    merged_df['slen'] = group.iloc[0].slen
    
    # Reorder columns
    merged_df = merged_df[[
        'sseqid', 'sstart', 'send', 'qseqid', 'stitle', 'sscinames', 'slen'
        ]]
    
    return merged_df


# Load masking table
mask_df = pd.read_csv("Recurrent_contaminant_regions.tsv", sep='\t')
mask_df["ACC"] = mask_df["ACC"].apply(clean_accession)

# Load reference FASTA
db_folder = "%s/viral_nuccore/viral_database.fasta" % args.virmet_db
ref_seqs = SeqIO.to_dict(SeqIO.parse(
    db_folder, 
    "fasta"
    ))

# Load target db
target_db = "%s/viral_nuccore/viral_db" % args.virmet_db

# Create a temporary FASTA file for all query subsequences
query_records = []
valid_rows = []
for _, row in mask_df.iterrows():
    acc = row["ACC"]
    start = row["Start"] - 1
    end = row["End"]

    if acc not in ref_seqs:
        print(f"Accession {acc} not found in reference FASTA")
        continue

    sub_seq = ref_seqs[acc].seq[start:end]

    # Remove sequences too short or too long to prevent masking undesired regions
    if args.min_fragment_size <= len(sub_seq) <= args.max_fragment_size:
        query_records.append(SeqRecord(sub_seq, id=f"{acc}_{start+1}_{end}", description=""))
        valid_rows.append(row)
    else:
        print(f"Skipping {acc} subsequence ({start+1}-{end}): length {len(sub_seq)} is too short or too long")

SeqIO.write(query_records, "query_seqs.fasta", "fasta")

# Replace mask_df with only valid rows
mask_df = pd.DataFrame(valid_rows).reset_index(drop=True)

# Run BLASTN
blast_output = "Sequence_similarity_search.tsv"

# Copy current environment and add BLASTDB to find the taxid names
env = os.environ.copy()
env["BLASTDB"] = args.virmet_db

blast_cmd = [
    "blastn", "-task", "blastn", "-query", "query_seqs.fasta", "-db", target_db,
    "-outfmt", "6 sseqid sstart send qseqid stitle sscinames pident qcovs slen", 
    "-max_target_seqs", "1000000",
    "-num_threads", str(args.threads),
    "-out", blast_output
]
result= subprocess.run(blast_cmd, env=env, capture_output=True, text=True)
print("STDOUT:", result.stdout)
print("STDERR:", result.stderr)

# Parse BLAST results
cols = ["sseqid", "sstart", "send", "qseqid", "stitle", "sscinames", "pident", "qcovs", "slen"]
blast_out = pd.read_csv(blast_output, sep='\t', header=None, names=cols)

blast_df_raw = blast_out[(blast_out["pident"] > 80) & (blast_out["qcovs"] > 80)].copy()

# Define patterns to exclude from masking (controls):
if args.controls == "Tunavirus T1, phage MS2":
    pattern_ic_ms2 = ("Emesvirus zinderi", "phage MS2", "virus MS2")
    pattern_ic_t1 = ("Tunavirus", "phage T1", "virus T1")
    all_patterns = pattern_ic_ms2 + pattern_ic_t1
else:
    all_patterns = tuple(item.strip() for item in args.controls.split(","))

# compile regex pattern (escaped and joined with |)
regex_pattern = "|".join(map(re.escape, all_patterns))

# Identify qseqids that match excluded genomes
excluded_qseqids = blast_df_raw.loc[
    blast_df_raw['sscinames'].str.contains(regex_pattern, case=False, na=False) |
    blast_df_raw['stitle'].str.contains(regex_pattern, case=False, na=False),
    'qseqid'
].unique()

# Filter raw file
blast_df = blast_df_raw[~blast_df_raw['qseqid'].isin(excluded_qseqids)].copy()

# Apply clean_accession again (just in case)
blast_df["sseqid"] = blast_df["sseqid"].apply(clean_accession)

# Map qseqid to mask_df rows
mask_df["qseqid"] = mask_df.apply(lambda r: f"{r['ACC']}_{r['Start']}_{r['End']}", axis=1)
# Drop rows with excluded qseqid
mask_df = mask_df[~mask_df["qseqid"].isin(excluded_qseqids)].copy()

# Save
mask_df = mask_df.drop(columns=["qseqid"])
mask_df.to_csv("Recurrent_contaminant_regions.tsv", sep='\t', index=False)

# For each sseqid in the outputs, merge all intervals
merged_df_list = []
for sseqid, group in blast_df.groupby('sseqid'):
    merged_intervals = merge_intervals(group)
    
    # Calculate total length to mask covered by merged intervals
    total_length = (merged_intervals['send'] - merged_intervals['sstart'] + 1).sum()
    slen = merged_intervals.iloc[0]['slen']
    
    # Check if total length < x% of slen to avoid masking the whole genome
    threshold = args.max_genome_perc
    if total_length <= threshold * slen or slen <= 1000:
        merged_df_list.append(merged_intervals)
        if slen <= 1000 and total_length > threshold * slen:
            print(f"Note! Masking {sseqid} is > {int(threshold*100)}% of genome ({total_length} bp)")
    else:
        print(f"Skipping {sseqid} because total length {total_length} > {int(threshold*100)}%")


# Combine all filtered sseqid results
final_df = pd.concat(merged_df_list, ignore_index=True)

final_df.to_csv(blast_output, sep='\t', index=False)

# Keep a copy of the original reference sequences
original_ref_seqs = copy.deepcopy(ref_seqs)

# Mask the database based on the final_df table
for target_id, group in final_df.groupby("sseqid"):
    if target_id not in ref_seqs:
        print(f"Target ID {target_id} not found in reference DB")
        continue

    record = ref_seqs[target_id]
    seq = str(record.seq)

    masked_seq = list(seq)

    for _, row in group.iterrows():
        s = min(row["sstart"], row["send"]) - 1  # BLAST is 1-based
        e = max(row["sstart"], row["send"])

        masked_seq[s:e] = ['N'] * (e - s)

    # Overwrite sequence
    new_seq = "".join(masked_seq)
    if new_seq != seq:
        record.seq = Seq(new_seq)

# Save masked version of FASTA2
SeqIO.write(ref_seqs.values(), "masked_viral_database.fasta", "fasta")

# Extract only masked regions
masked_records = []

for target_id, group in final_df.groupby("sseqid"):
    if target_id not in original_ref_seqs:
        continue

    record = original_ref_seqs[target_id]
    seq = str(record.seq)

    for _, row in group.iterrows():
        s = min(row["sstart"], row["send"]) - 1
        e = max(row["sstart"], row["send"])
        masked_fragment = seq[s:e]

        masked_records.append(
            SeqRecord(
                Seq(masked_fragment),
                id=f"{target_id}:{s+1}-{e}",
                description=f"masked region from {target_id} ({s+1}-{e})"
            )
        )

SeqIO.write(masked_records, "Homologous_regions_masked.fasta", "fasta")

# Clean
os.remove("query_seqs.fasta")