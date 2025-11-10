#!/usr/bin/env bash

set -e -o pipefail

###################
## Assign inputs ##
###################

# Define usage of the script
function print_usage {
  printf """
Usage: 
VirMask.sh -r <mNGS_runs> [options]
Options:
	[-h or --help]
	[-d or --virmet_db]
	[-s or --across_sample_freq]
	[-a or --within_acc_dominance]
	[-l or --min_fragment_size]
	[-u or --max_fragment_size]
	[-g or --max_genome_perc]
	[-c or --controls]
	[-o or --output]
	[-t or --threads]
"""
}
# Describe usage of the tool and provide help
function print_help {
  print_usage
  printf """
Options:
	-h, --help:
		Show this help message and exit.
	-d, --virmet_db:
		Path to the VirMet database.
		Default: /data/virmet_databases_update/
	-s, --across_sample_freq:
		Across-sample frequency threshold.
		Default: 0.2 (20%%)
	-a, --within_acc_dominance
		With-accesion dominance threshold.
		Default: 0.2 (20%%)
	-l, --min_fragment_size:
		Minimum size (bp) of the fragment to be masked.
		Default: 21
	-u, --max_fragment_size:
		Maximum size (bp) of the fragment to be masked.
		Default: 499
	-g, --max_genome_perc:
		Maximum percentage of genome to mask.
		Default: 0.05 (5%%)
	-c, --controls:
		Pattern to identify viruses that should be prevented from masking.
		If more than one is provided, add them within quotes and separated by commas.
		This avoids masking your internal controls, which appear in all samples.
		Default: 'Tunavirus T1, phage MS2'
	-o, --output:
		Path to the location where the outputs should be stored.
		Default: current folder, which is ./
	-t, --threads:
		Number of threads to use.
		Default: 24
Required arguments:
	-r, --mNGS_runs:
		Path to a text file (.txt) listing all mNGS runs to consider (one run per row).
"""
}
# Define defaults
virmet_db="/data/virmet_databases_update/"
across_sample_freq=0.2
within_acc_dominance=0.2
min_fragment_size=21
max_fragment_size=499
max_genome_perc=0.05
controls="Tunavirus T1, phage MS2"
output="./"
threads=24

# Define inputs
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--mNGS_runs) mNGS_runs="$2"; shift 2 ;;
        -d|--virmet_db) virmet_db="$2"; shift 2 ;;
        -s|--across_sample_freq) across_sample_freq="$2"; shift 2 ;;
        -a|--within_acc_dominance) within_acc_dominance="$2"; shift 2 ;;
        -l|--min_fragment_size) min_fragment_size="$2"; shift 2 ;;
        -u|--max_fragment_size) max_fragment_size="$2"; shift 2 ;;
        -g|--max_genome_perc) max_genome_perc="$2"; shift 2 ;;
        -c|--controls) controls="$2"; shift 2 ;;
        -o|--output) output="$2"; shift 2 ;;
        -t|--threads) threads="$2"; shift 2 ;;
        -h|--help) print_help; exit 0 ;;
        *) echo "Unknown option $1"; print_usage; exit 1 ;;
    esac
done

# Validate inputs
if [ -z "$mNGS_runs" ]; then
    echo "Error: -r or --mNGS_runs is required."
    print_usage
    exit 1
fi
if [ ! -f "$mNGS_runs" ]; then
    echo "Error: mNGS_runs file '$mNGS_runs' does not exist."
    exit 1
fi
# Validate threads, min_fragment_size and max_fragment_size
if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: --threads must be an integer"
    exit 1
fi
if ! [[ "$min_fragment_size" =~ ^[0-9]+$ ]]; then
    echo "Error: --min_fragment_size must be an integer"
    exit 1
fi
if ! [[ "$max_fragment_size" =~ ^[0-9]+$ ]]; then
    echo "Error: --max_fragment_size must be an integer"
    exit 1
fi

### Define paths and convert relative to absolute paths
script_dir=$( dirname "$(readlink -f "$0")" )
mNGS_runs=$( readlink -f "$mNGS_runs" )
virmet_db=$( readlink -f "$virmet_db" )
output=$( readlink -f "$output" )

### Print arguments
echo -e 'mNGS_runs: ' "$mNGS_runs"
echo -e 'virmet_db ' "$virmet_db"
echo -e 'across_sample_freq: ' "$across_sample_freq"
echo -e 'within_acc_dominance: ' "$within_acc_dominance"
echo -e 'min_fragment_size: ' "$min_fragment_size"
echo -e 'max_fragment_size: ' "$max_fragment_size"
echo -e 'max_genome_perc: ' "$max_genome_perc"
echo -e 'controls: ' "$controls"
echo -e 'output: ' "$output"
echo -e 'threads: ' "$threads"

# Set output directory
cur_dir_name="NC_db_new_$(date +'%Y_%m_%d')"
working_dir=${output}/${cur_dir_name}

if [ -d $working_dir ] 
then
    echo "Directory $working_dir exists. Delete or rename this directory."
    exit 1
fi

mkdir $working_dir
echo "Directory $working_dir created."

#Change to output directory for simplicity
initial_dir=$(pwd)
cd $working_dir

# Identify human contamination
${script_dir}/Human_contaminants/Human_masking.py --virmet_db "${virmet_db}" --threads ${threads}

# Move the database without human contamination (partly masked DB) to the DB folder
mkdir ${virmet_db}/viral_nuccore/initial_unmasked_db
mv ${virmet_db}/viral_nuccore/viral_database.fasta ${virmet_db}/viral_nuccore/initial_unmasked_db
mv ${virmet_db}/viral_nuccore/viral_db* ${virmet_db}/viral_nuccore/initial_unmasked_db
cp "viral_db_masked.fasta" ${virmet_db}/viral_nuccore/viral_database.fasta

# Index the partly masked DB
target_dir=${virmet_db}/viral_nuccore
makeblastdb -in ${target_dir}/viral_database.fasta \
    -dbtype nucl -hash_index \
    -title "Masked viral db indexed" \
    -out ${target_dir}/viral_db \
    -logfile ${target_dir}/blast.log \
    -parse_seqids -taxid_map ${target_dir}/viral_accn_taxid.dmp

# Find and mask recurrent contaminants
${script_dir}/Recurrent_contaminants/Recurrent_masking.py \
    --mNGS_runs "${mNGS_runs}" \
    --virmet_db "${virmet_db}" \
    --across_sample_freq "${across_sample_freq}" \
    --within_acc_dominance "${within_acc_dominance}" \
    --min_fragment_size "${min_fragment_size}" \
    --max_fragment_size "${max_fragment_size}" \
    --max_genome_perc "${max_genome_perc}" \
    --controls "${controls}" \
    --threads "${threads}"

# Remove intermediate files
rm -r virmet_output*
rm -r masked_to_chr*
rm -r chr*/

# Copy masked database to the database folder
cp masked_viral_database.fasta ${virmet_db}/viral_nuccore/

# Move back to the initial directory
cd $initial_dir

# Print final message
echo "VirMask finished successfully! Good luck using your masked database!"