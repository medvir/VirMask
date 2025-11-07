# Installing VirMask

VirMask code can be downloaded from our [GitHub](https://github.com/medvir/VirMask/) as follows:

`git clone https://github.com/medvir/VirMask.git $HOME/VirMask`

After downloading the code, we highly advise to add VirMask into your `$PATH`
in order to be able to call it easily from any folder:

```
# Add VirMask to your PATH
echo "export PATH=$HOME/VirMask:$PATH" >> $HOME/.bashrc
echo "export PATH=$HOME/VirMask/Human_contaminants:$PATH" >> $HOME/.bashrc
echo "export PATH=$HOME/VirMask/Recurrent_contaminants:$PATH" >> $HOME/.bashrc
source $HOME/.bashrc
```

The tool is designed for a Linux environment with several software 
and dependencies installed and accessible in your `$PATH`. See information below.

## Human Masking

To run the [`Human masking`](./Human_contaminants.md) function, users will need:

* **Python 3.8 or newer** — fully tested with Python 3.12.

* **Snakemake (≥7.0)** — fully tested with Snakemake 7.32.4.

* **samtools (≥1.10)** — fully tested with samtools 1.21.

* **bedtools (≥2.30)** — fully tested with bedtools v2.31.1.

* **bwa (≥0.7.17)** — fully tested with bwa 0.7.19-r1273

* **art_illumina** — fully tested with version 2.5.8.

* **awk, grep, sort, uniq, gzip** (standard UNIX tools, usually preinstalled).

Make sure all these tools can be called from the command line before running the scripts.

In addition, install the required Python packages using `pip`:

`pip install pandas`

With that, you are ready to start the human masking!

## Contaminants Masking

To run the [`Contaminants masking`](./Recurrent_contaminants.md) function, users will need:

* **Python 3.8 or newer** — fully tested with Python 3.12.

* **blastn (BLAST+ ≥2.12.0)** — fully tested with version 2.16.0+.

* **VirMet 2.0** — fully tested with VirMet 2.0.0.

* **awk, grep, sort, uniq, gzip** (standard UNIX tools, usually preinstalled).

Make sure all these tools can be called from the command line before running the scripts.

In addition, install the required Python packages using `pip`:

```
pip install pandas # if not done before
pip install biopython
```

With that, you are ready to start the masking of recurrent contaminants!

## Masking All Together: VirMask

To run [`VirMask`](./Virmask.md), including both human masking and recurrent contaminants masking,
users will need to have all software and dependencies listed above for each of
the two steps.

## Database

To run [`VirMask`](./Virmask.md), the tool expects users to have, at least,
a human database and a viral database (only for the
 [`Human masking`](./Human_contaminants.md)).

If they are also planning to use the [`Contaminants masking`](./Recurrent_contaminants.md),
they will need the whole [`VirMet`](https://github.com/medvir/VirMet) database,
which can be downloaded as follows:

```
virmet fetch --viral n
virmet fetch --human
virmet fetch --bact_fungal
virmet fetch --bovine
```
To make the VirMet database fully functional, they will also need to index it.

`virmet index --viral n --human --bact_fungal --bovine`

In either case, the database directory should have the following structure:

```
# Database needed for the Human Masking only:
virmet_databases_update/
├── viral_nuccore/
│   └── viral_database.fasta
└── human/
    └── fasta/
        └── GRCh38.fasta.gz

# Database needed for the Recurrent Contaminats Masking:
virmet_databases_update/
├── viral_nuccore/
│   ├── viral_database.fasta
│   ├── viral_accn_taxid.dmp
│   └── viral_seqs_info.tsv
├── human/
│   ├── fasta/
│   │   └── GRCh38.fasta.gz
│   └── bwa/
│       └── bwa_files
├── bovine/
│   ├── fasta/
│   │   └── ref_Bos_taurus.fasta.gz
│   └── bwa/
│       └── bwa_files
├── bact_fungi/
│   ├── library/
│   │   ├── bacteria
│   │   │   └── library.fna
│   │   └── fungi
│   │       └── library.fna
│   └── taxonomy/
│       └── taxdump.tar.gz
├── names.dmp.gz
└── nodes.dmp.gz
```
Please, note that the name of the database doesn't need to be necessarily
`virmet_databases_update`. However, this is the default name taken by
[`VirMask`](./Virmask.md) if not specified otherwise.

**Enjoy using VirMask!**