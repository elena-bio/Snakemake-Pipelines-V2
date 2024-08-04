
## Snakemake Workflow

Snakemake is a workflow management system that aims to make the definition and execution of data analysis workflows easy, reproducible, and scalable. It uses a Python-based language to specify workflows in terms of rules that describe 

how to derive output files from input files.

For a comprehensive introduction and installation of Snakemake, refer to the official tutorial: (https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

## Connecting to Computerome Server 

To connect to the Computerome server, use the following command:
```
$ ssh elehos@ssh.computerome.dk
```

For more information about the Computerome server, refer to this [link](https://github.com/elena-bio/snakemake-pipelines-v2/blob/main/computerome%20server.md).

## Environment setup


## Sample data

1. Create a directory for the task
```
   $ mkdir Snakemak
```

1. Download the new dataset from this link https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9205427/#sec4.3.2.
   
1. Synchronize the input files using `rsync`: To do this between your local machine and the remote server use the following command:

```
   `$ rsync "/path/to/local/directory/" user@transfer.computerome.dk:/path/to/remote/directory/
```

* prepared bam file and refrence genome fasta file 
* HG002_35x_PacBio_14kb-15kb.fastq.gz and refrence genome fasta file



## Snakemake workflows 

1. A Snakemake workflow is defined by specifying rules in a Snakefile. To begin, create a new file called Snakefile in your working directory using any text editor you prefer. Rules decompose the workflow into small steps by specifying
  
1. how to create sets of output files from sets of input files. Snakemake automatically determines the dependencies between the rules by matching file names.

usage command
```
   $ nano Snakefile
```

### Reads alignment

The first Snakemake rule maps reads of a given sample to a given reference genome. For alignment, I used the tool minimap2.

`Minimap2` is a versatile alignment tool designed for mapping DNA or RNA sequence reads to a reference genome or for aligning long reads to each other.

`Minimap2` is publicly available and open source.Minimap2 seamlessly works with gzip’d FASTA and FASTQ formats as input, no need to convert between FASTA and FASTQ or decompress gzip’d files first.

Snakemake accepts rule names as targets if the requested rule does not include wildcards, providing flexibility in specifying targets. To manage all desired results efficiently, a rule named **all** can be defined at the top of the 

workflow, including all typically desired target files as input files.The current version used is `Minimap2 2.28-r1209`.

In the Snakefile, I defined the following rule:

usage command
```
rule all:
    input:
        "outputs/output.sam"

rule minimap2:
    input:
        HGOO2="HG002_35x_PacBio_14kb-15kb.fastq.gz",
        ref="GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta"
    output:
        "outputs/output.sam"
    shell:
        "minimap2 -ax map-hifi {input.ref} {input.HGOO2} > {output}"

```

To run jobs in the HPC environment, you must submit them to a queuing system. This system allocates compute resources fairly and efficiently among users and schedules job execution. Below is a generic job submission script for Snakemake rules, designed to avoid the need to create a separate script for each rule. This script is essential for effectively executing tasks within the HPC environment.

```
#!/bin/bash
### Account information
#PBS -W group_list=cu_10160 -A cu_10160
### Job name
#PBS -N snakemake
###Output files
#PBS -e snakemake.err
#PBS -o snakemake.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40
### Memory
#PBS -l mem=150gb
###Requesting time
#PBS -l walltime=24:00:00

echo "Working direcory is $PBS_O_WORKDIR"
cd $PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes

### This ensures that any environment variables, aliases, or functions defined$
source ~/.bashrc

#echo pwd: $(pwd)
#echo 'Testing Queue Submission!' > ./testing.out

# Prepare the environment
#conda activate pbhoney-env
conda activate snakemake-env

# Execute your command
snakemake --cores 40 

```
To execute the `snakemake.sh` script, I used the following command in terminal:
```
$ qsub -W group_list=cu_10160 -A cu_10160 -l nodes=1:ppn=8,mem=40gb,walltime=12:00:00 ./snakemake.sh
```
To check the status of the job that has been submitted 
```
qstat -a; ll -h
```

### Converting sam file to bam file

The second rule converts a SAM file to a BAM file . To speeds up processing, saves space, and allows for quick access to specific data regions.

usage command
```
rule all:
    input:
        "outputs/output.bam"

rule minimap2:
    input:
        HGOO2="HG002_35x_PacBio_14kb-15kb.fastq.gz",
        ref="GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta"
    output:
        "outputs/output.sam"
    shell:
        "minimap2 -ax map-hifi {input.ref} {input.HGOO2} > {output}"
rule sam_to_bam:
    input:
        "outputs/output.sam"
    output:
        "outputs/output.bam"
    shell:
        "samtools view -bS {input} > {output}"

```
To get the desired output, I executed the `snakemake.sh` script again using the command above.

### Sorting read alignments

For later steps, we need the read alignments in the BAM files to be sorted. This can be achieved with the **samtools** `sort` command.

```
rule all:
    input:
        "outputs/sorted.bam"

rule minimap2:
    input:
        HGOO2="HG002_35x_PacBio_14kb-15kb.fastq.gz",
        ref="GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta"
    output:
        "outputs/output.sam"
    shell:
        "minimap2 -ax map-hifi {input.ref} {input.HGOO2} > {output}"
rule sam_to_bam:
    input:
        "outputs/output.sam"
    output:
        "outputs/output.bam"
    shell:
        "samtools view -bS {input} > {output}"
rule samtools_sort:
    input:
        "outputs/output.bam"
    output:
        "outputs/sorted.bam"
    shell:
        "samtools sort -T outputs/sorted -O bam {input} > {output}"

```

### Indexing read alignment 

Next, we need to use `samtools` again to index the sorted read alignments so that we can quickly access reads by the genomic location they were mapped to. This can be done with the following rule:

```
rule all:
    input:
        "outputs/sorted.bam.bai"

rule minimap2:
    input:
        HGOO2="HG002_35x_PacBio_14kb-15kb.fastq.gz",
        ref="GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta"
    output:
        "outputs/output.sam"
    shell:
        "minimap2 -ax map-hifi {input.ref} {input.HGOO2} > {output}"
rule sam_to_bam:
    input:
        "outputs/output.sam"
    output:
        "outputs/output.bam"
    shell:
        "samtools view -bS {input} > {output}"
rule samtools_sort:
    input:
        "outputs/output.bam"
    output:
        "outputs/sorted.bam"
    shell:
        "samtools sort -T outputs/sorted -O bam {input} > {output}"
rule samtools_index:
    input:
        "outputs/sorted.bam"
    output:
        "outputs/sorted.bam.bai"
    shell:
        "samtools index {input}"

```

### Extracting chromosome 6 

```
rule all:
    input:
        "outputs/chromosome6.bam"

rule minimap2:
    input:
        HGOO2="HG002_35x_PacBio_14kb-15kb.fastq.gz",
        ref="GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta"
    output:
        "outputs/output.sam"
    shell:
        "minimap2 -ax map-hifi {input.ref} {input.HGOO2} > {output}"
rule sam_to_bam:
    input:
        "outputs/output.sam"
    output:
        "outputs/output.bam"
    shell:
        "samtools view -bS {input} > {output}"
rule samtools_sort:
    input:
        "outputs/output.bam"
    output:
        "outputs/sorted.bam"
    shell:
        "samtools sort -T outputs/sorted -O bam {input} > {output}"
rule samtools_index:
    input:
        "outputs/sorted.bam"
    output:
        "outputs/sorted.bam.bai"
    shell:
        "samtools index {input}"
rule samtools_extract:
    input:
        bam="outputs/sorted.bam",
        bai="outputs/sorted.bam.bai"
    output:
        "outputs/chromosome6.bam"
    shell:
        "samtools view -b {input.bam} chr6 > {output}"

```


### Runing pbhoney for variant calling 

### Runing sniffles for variant calling

