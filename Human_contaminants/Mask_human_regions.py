#!/usr/bin/env python3
"""
This script detects and mask viral sequences contaminated with human reads.
1. Simulates human reads per chromosome.
2. Aligns them to a viral database.
3. Generates coverage BED files.
4. Sequentially masks the viral DB using those BED files.
"""

# Import standard libraries
import os

# Parse command-line arguments
VIRAL_DB = os.path.join(config["virmet_db"], "viral_nuccore", "viral_database.fasta")
HUMAN_CHR_DB = os.path.join(config["virmet_db"], "human", "fasta", "GRCh38.fasta.gz")

DIR = [
    'chr1','chr2','chr3','chr4','chr5', 'chr6', 'chr7','chr8', 'chr9',
    'chr10', 'chr11','chr12', 'chr13','chr14', 'chr15','chr16', 'chr17',
    'chr18', 'chr19','chr20', 'chr21', 'chr22', 'chrX','chrY', 'chrM'
    ]
        
bwa_idx_ext = [".bwt", ".amb", ".ann", ".pac", ".sa"]

# Define functions

def multiext(base, *extensions):
    """Return a list of filenames by appending each extension."""
    return [base + ext for ext in extensions]

def prev_chr(ch):
    """Return the previous chromosome name (for sequential masking)."""
    idx = DIR.index(ch)
    return DIR[idx - 1] if idx > 0 else None

# Define rules

rule all:
    input:
        expand('{dir}/{dir}_ref_count_mapped_newDB.txt', dir=DIR),
        expand('{dir}/{dir}_coverage.bed', dir=DIR),
        "viral_db_masked.fasta"


rule idx_human:
    input:
        human_ref = HUMAN_CHR_DB
    output:
        idx_hum = multiext(HUMAN_CHR_DB, ".fai", ".gzi")
    shell:
        """
        samtools faidx {input.human_ref}
        """

rule extract_chr:
    input:
        idx_hum = rules.idx_human.output.idx_hum,
        human_ref = HUMAN_CHR_DB
    output:
        fasta_file = "{dir}/{dir}.fasta"
    params:
        folder = "{dir}"
    shell:
        """
        samtools faidx {input.human_ref} {params.folder} > {output.fasta_file}
        """

rule art_sim:
    input:
        "{dir}/{dir}.fasta"
    output:
        "{dir}/single_read_f3_errFree.sam"
    params:
        template = "{dir}/single_read_f3"
    shell:
        """
        art_illumina -ss MSv3 -i {input} -l 150 -f 3 --errfree --noALN -o {params.template}
        """

rule create_errfree_bam:
    input:
        "{dir}/single_read_f3_errFree.sam"
    output:
        "{dir}/single_read_f3_errFree.bam"
    shell:
        "samtools view -S -b {input} > {output}"


rule extract_fq:
    input:
        "{dir}/single_read_f3_errFree.bam"
    output:
        "{dir}/single_read_f3_errFree.fastq"
    shell:
        "bedtools bamtofastq -i {input} -fq {output}"

rule ref_bwa_index:
    input:
        VIRAL_DB,
    output:
        multiext(VIRAL_DB, *bwa_idx_ext)
    shell:
        """
        bwa index {input}
        """

rule align:
    input:
        DATA="{dir}/single_read_f3_errFree.fastq",
        INDEX=multiext(VIRAL_DB, *bwa_idx_ext),
    output:
        "{dir}/aln_f3.sam"
    params:
        VIRAL_DB=VIRAL_DB,
    shell:
        """
        bwa mem {params.VIRAL_DB} {input.DATA} > {output}
        """

rule sort:
    input:
        "{dir}/aln_f3.sam"
    output:
        "{dir}/aln_f3_sorted.bam"
    shell:
        "samtools sort -o {output} {input}"
  

rule bedtools_mapped:
    input:
        "{dir}/aln_f3_sorted.bam" 
    output:
        "{dir}/{dir}_ref_count_mapped_newDB.txt"
    shell:
        r"""
        bedtools bamtobed -i {input} | awk -F "\t" '{{print $1}}' | sort | uniq -c | sort -s -rn -k 1,1 > {output}
        """


rule bedfiles:
    input:
        "{dir}/aln_f3_sorted.bam" 
    output:
        "{dir}/{dir}_coverage.bed"
    shell:
        """
        bedtools genomecov -bg -ibam {input} > {output}
        """

rule mask_one_chr:
    """Mask viral DB sequentially for each chromosome."""
    input:
        cur_DB=lambda wildcards: f"masked_to_{prev_chr(wildcards.dir)}.fasta"
        if prev_chr(wildcards.dir) else VIRAL_DB,
        bed_file=lambda wildcards: f"{wildcards.dir}/{wildcards.dir}_coverage.bed"
    output:
        masked_DB="masked_to_{dir}.fasta"
    shell:
        """
        bedtools maskfasta -fi {input.cur_DB} -bed {input.bed_file} -mc N -fullHeader -fo {output.masked_DB}
        """

rule final_masked_db:
    input:
        f"masked_to_{DIR[-1]}.fasta"
    output:
        "viral_db_masked.fasta"
    shell:
        "cp {input} {output}"
