# Remove Recurrent Contaminants

In viral mNGS diagnosis, a common challenge is the frequent classification of sequences
as viral that appear across unrelated samples. These spurious hits often arise
from non-viral elements such as rRNA, vector contamination and sequencing adapters
within the viral database, along with viral-derived sequences likely introduced
into the mNGS samples during sample preparation (e.g., CMV enhancers).

To improve diagnostic specificity and reduce false-positive identifications, we
recommend to systematically mask these recurrent contaminants.
This is the main goal of the Contaminants Masking step implemented into [`VirMask`](./Virmask.md).

In order to use the Contaminants Masking option alone (instead of 
[`VirMask`](./Virmask.md) as a whole) users need to run:

`Recurrent_masking.py --mNGS_runs <mNGS_to_analyse.txt> [options]`

Information on the input parameters can be obtained with the `-h` flag:

```
usage: Recurrent_masking.py [-h] --mNGS_runs MNGS_RUNS [options]

Detection and masking of recurrent contaminants within the viral database.

Options:
  -h, --help             show this help message and exit
  --mNGS_runs MNGS_RUNS
                         Path to a .txt file listing all mNGS runs to consider (one run per row)
  --virmet_db VIRMET_DB
                         Path to the VirMet database (default: /data/virmet_databases_update/)
  --across_sample_freq   ACROSS_SAMPLE_FREQ
                         Across-sample frequency threshold (default: 0.2 = 20%)
  --within_acc_dominance WITHIN_ACC_DOMINANCE
                         Within-accession dominance threshold (default: 0.2 = 20%)
  --min_fragment_size    MIN_FRAGMENT_SIZE
                         Minimum size of the fragment to be masked (default: 21 bp)
  --max_fragment_size    MAX_FRAGMENT_SIZE
                         Maximum size of the fragment to be masked (default: 499 bp)
  --max_genome_perc      MAX_GENOME_PERC
                         Maximum percentage of genome to mask (default: 0.05 = 5%)
  --controls CONTROLS    Pattern to identify viruses that should be prevented from masking. 
                         If more than one is provided, add them within quotes and separated by commas.
                         Example: 'Tunavirus T1, phage MS2'.
                         This avoids masking your internal controls, which appear in all samples.
  --threads THREADS      Number of threads to use (default: 24)
```

Please, note that the outputs will be stored at the working directory. Therefore,
it is highly recommended to create an output folder first and set it as the
current working directory:

```
mkdir contaminants_masking_outputs
cd contaminants_masking_outputs
Recurrent_masking.py --mNGS_runs mNGS_to_analyse.txt
```

Users who use [`VirMask`](./Virmask.md) pipeline don't need to worry about that. 
VirMask already offers the possibility to specify an output folder and creates a subfolder
there to store the outputs.

## Outputs

Running the previous command should provide the following output files:

* **Recurrent_contaminant_regions.tsv**: recurrent contaminants identified in the samples.
* **Sequence_similarity_search.tsv**: homologous regions that were ultimately masked.
* **Masked_viral_database.fasta**: the masked viral reference database (FASTA file).
* **Homologous_regions_masked.fasta**: the extracted masked fragments used for validation.

**Enjoy removing Recurrent Contaminants!**