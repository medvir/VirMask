# Welcome to VirMask

VirMask is a novel computational tool created to systematically mask non-viral 
elements (e.g., human, rRNA sequences and other artifacts) within viral reference 
databases, as well as virus-derived sequences that are recurrently detected in 
mNGS samples, likely introduced during sample preparation or sequencing workflows.

The ultimate goal of the tool is to improve reliability of viral mNGS diagnosis.

To download the code, you can clone this repository:

`git clone https://github.com/medvir/VirMask.git $HOME/VirMask`

Information on how to install the tool can be found in the 
[`Installation`](https://medvir.github.io/VirMask/Installation/) section of the 
[`VirMask website`](https://medvir.github.io/VirMask/).

VirMask consists of two main parts:

* [`Human masking`](https://medvir.github.io/VirMask/Human_contaminants/): detects 
and masks viral refrence genomes that contain human sequences.

* [`Contaminants masking`](https://medvir.github.io/VirMask/Recurrent_contaminants/): 
identifies and masks recurrent contaminants within the viral database, such as rRNA 
seqences, vector contamination and sequencing adapters present in the viral database.

Users can run each of these two functions individually if only human sequences or 
recurrent contaminants need to be removed, or alternatively (and highly recommended), 
they can run VirMask as a whole:

* [`VirMask`](https://medvir.github.io/VirMask/Virmask/): spots and masks
human sequences and recurrent contaminants from the viral database.

Some help can be obtained with `VirMask.sh -h`.

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

**Enjoy using VirMask!**

<p align="center">
  <img src="documentation/docs/assets/logo.svg" width="50%">
</p>
