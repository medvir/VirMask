#!/usr/bin/env python3
"""
This script provides a summary of the masked human-associated viral sequences.
1. Collects viral read counts mapped to each human chromosome.
2. From individual files, it merges them into a single table.
3. Calculates the total reads per virus, retrieves virus titles from the viral
database, and outputs a combined CSV file with read counts and virus information.
"""

# Import standard libraries
import argparse
import os
import subprocess

# Import third-party libraries
import pandas as pd

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Summarise the human sequences present in the viral database."
    )
parser.add_argument(
    "--virmet_db",
    type=str, 
    default="/data/virmet_databases_update/",
    help="Path to the VirMet database (default: /data/virmet_databases_update/)"
)
args = parser.parse_args()

VIRAL_DB = os.path.join(args.virmet_db, "viral_nuccore", "viral_database.fasta")

# Do the summary
for i in range(1,26):
    if i == 23:
        i = 'X'
    elif i == 24:
        i = 'Y'
    elif i == 25:
        i = 'M'
    file_newDB = "chr%s/chr%s_ref_count_mapped_newDB.txt"%(i,i)
    col_read_name = 'reads_nb_chr%s'%i
    try:
        ref_newDB = pd.read_csv(file_newDB, sep=" ", names=[col_read_name,'virus'], skipinitialspace=True)
    except:
        print( "check if file %s exists. " %file_newDB)
    
    col_read_name = 'reads_nb_chr%s_y'%i
    if i == 1:
        chrs_merge_df = ref_newDB
    else:
        chrs_merge_df = chrs_merge_df.merge(ref_newDB, how='outer', on='virus')

chrs_merge_df['sum'] = chrs_merge_df.sum(axis=1, numeric_only=True)

title_list = []
for index, row in chrs_merge_df.iterrows():
    title = subprocess.check_output(["grep", row['virus'], VIRAL_DB])
    title_list.append(title)
    
string_list=[x.decode('utf-8') for x in title_list]
str_list =[x.strip().strip('\n') for x in string_list]
chrs_merge_df['title'] = str_list

first_column = chrs_merge_df.pop('virus')
chrs_merge_df.insert(0, 'virus', first_column)
chrs_merge_df.sort_values(by=['sum'], inplace=True)
output_file = 'chrs_mapped_ref_newDB.csv'
chrs_merge_df.to_csv(output_file, index=False)
