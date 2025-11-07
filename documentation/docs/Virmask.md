# Remove Human Sequences and Recurrent Contaminants

VirMask has been developed as a novel computational strategy to improve the specificity 
of viral mNGS diagnosis by masking problematic regions within viral reference databases.

Upon masking human sequences (likely attributable to the viral hosts), our approach mines 
prior mNGS runs to identify sequences that appear frequently across a high proportion 
of samples or that dominate the read counts for particular strains, suggesting they may
arise from technical artifacts or non-viral sources rather than genuine infections.
These sequences are then subjected to alignment-based similarity searches to identify
homologous regions in the reference database, which are subsequently masked as well.

VirMask's targeted masking reduces the likelihood that reads derived from non-specific
or contamination-prone regions will spuriously align to viral genomes.

In order to use VirMask, users need to run:

`VirMask.sh -r <mNGS_runs> [options]`

Information on the input parameters can be obtained with the `-h` flag:

```
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

Options:
        -h, --help:
                Show this help message and exit.
        -d, --virmet_db:
                Path to the VirMet database.
                Default: /data/virmet_databases_update/
        -s, --across_sample_freq:
                Across-sample frequency threshold.
                Default: 0.2 (20%)
        -a, --within_acc_dominance
                With-accesion dominance threshold.
                Default: 0.2 (20%)
        -l, --min_fragment_size:
                Minimum size (bp) of the fragment to be masked.
                Default: 21
        -u, --max_fragment_size:
                Maximum size (bp) of the fragment to be masked.
                Default: 499
        -g, --max_genome_perc:
                Maximum percentage of genome to mask.
                Default: 0.05 (5%)
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
```

Please, note that a subfolder will be created in the specified output directory, and all the
outputs will be stored there. If no output directory is provided, the current working directory
will be considered.

## Outputs

Running the previous command should provide the following output files:

* **Viral_db_masked.fasta**: viral database with masked human sequences only.
* **Chrs_mapped_ref_newDB.csv**: summary table of the masked human-associated viral sequences.
* **Recurrent_contaminant_regions.tsv**: recurrent contaminants identified in the samples.
* **Sequence_similarity_search.tsv**: homologous regions that were ultimately masked.
* **Masked_viral_database.fasta**: the fully masked viral reference database (FASTA file).
* **Homologous_regions_masked.fasta**: the extracted masked fragments used for validation.

Besides these important files, VirMask will also create several other intermediate
files are there for inspection only and can be removed by the users
at any time.

**Enjoy removing contaminants with VirMask!**

![](assets/logo.svg){: style="display:block; margin-left:auto; margin-right:auto; width:50%;" }