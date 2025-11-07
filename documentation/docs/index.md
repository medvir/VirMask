# Welcome to VirMask

VirMask is a novel computational tool created to systematically mask non-viral 
elements (e.g., human, rRNA sequences and other artifacts) within viral reference 
databases, as well as virus-derived sequences that are recurrently detected in 
mNGS samples, likely introduced during sample preparation or sequencing workflows.

The ultimate goal of the tool is to improve reliability of viral mNGS diagnosis.

To see and download the code, visit our [GitHub](https://github.com/medvir/VirMask/).
Information on how to install the tool can be found in the [`Installation`](./Installation.md) section.

VirMask consists of two main parts:

* [`Human masking`](./Human_contaminants.md): detects and masks viral refrence 
genomes that contain human sequences.

* [`Contaminants masking`](./Recurrent_contaminants.md): identifies and masks 
recurrent contaminants within the viral database, such as rRNA seqences, vector 
contamination and sequencing adapters present in the viral database.

Users can run each of these two functions individually if only human sequences or 
recurrent contaminants need to be removed, or alternatively (and highly recommended), 
they can run VirMask as a whole:

* [`VirMask`](./Virmask.md): spots and masks human sequences and 
recurrent contaminants from the viral database.

Some help can be obtained with `VirMask.sh -h`.

**Enjoy using VirMask!**

![](assets/logo.svg){: style="display:block; margin-left:auto; margin-right:auto; width:50%;" }