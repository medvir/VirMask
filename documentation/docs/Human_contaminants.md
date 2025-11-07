# Remove Human Contaminants

Very often, viral reference databases (e.g., viral genomes downloaded from the NCBI RefSeq 
database) may contain non-viral contaminants, in particular host human DNA.
As a first step, we recommend masking all these human contaminants to dramatically reduce
the number of spurious alignments and false positives.

Even though there are multiple software available for that, we also implemented this option
into VirMask, aiming to provide users with an end-to-end tool for the masking process.

In order to use the Human Masking option alone (instead of 
[`VirMask`](./Virmask.md) as a whole) users need to run:

`Human_masking.py [options]`

Information on the input parameters can be obtained with the `-h` flag:

```
usage: Human_masking.py [-h] [--virmet_db VIRMET_DB] [--threads THREADS]

Detection and masking of human sequences within the viral database.

Options:
  -h, --help            show this help message and exit
  --virmet_db VIRMET_DB
                        Path to the VirMet database (default: /data/virmet_databases_update/)
  --threads THREADS     Number of threads to use (default: 24)
```

Please, note that the outputs will be stored at the working directory and will contain
multiple intermediate files. Therefore, it is highly recommended to create an output folder
first and set it as the current working directory:

```
mkdir human_masking_outputs
cd human_masking_outputs
Human_masking.py
```

Users who use [`VirMask`](./Virmask.md) pipeline don't need to worry about that. 
VirMask already offers the possibility to specify an output folder and creates a subfolder
there to store the outputs.

## Outputs

Running the previous command should provide several outputs, some of which are intermediate
files. The most important results are:

* **Viral_db_masked.fasta**: viral database with masked human sequences.
* **Chrs_mapped_ref_newDB.csv**: summary table of the masked human-associated viral sequences.

The other intermediate files are there for inspection only and can be removed by the users
at any time. They show the viral sequences that matched each of the human
chromosomes (including chromosomes 1 to 22, X, Y, M).

**Enjoy removing Human Sequences!**