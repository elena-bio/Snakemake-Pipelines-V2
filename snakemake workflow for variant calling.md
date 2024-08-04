
## Snakemake Workflow

Snakemake is a workflow management system that aims to make the definition and execution of data analysis workflows easy, reproducible, and scalable. It uses a Python-based language to specify workflows in terms of rules that describe how to derive output files from input files.
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

   $ mkdir Snakemak

1. Download the new dataset from this link https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9205427/#sec4.3.2.
1. Synchronize the input files using `rsync`: To do this between your local machine and the remote server use the following command:

   `$ rsync "/path/to/local/directory/" user@transfer.computerome.dk:/path/to/remote/directory/

* prepared bam file and refrence genome fasta file 
* HG002_35x_PacBio_14kb-15kb.fastq.gz and refrence genome fasta file



## Snakemake workflows 

1. A Snakemake workflow is defined by specifying rules in a Snakefile. To begin, create a new file called Snakefile in your working directory using any text editor you prefer. Rules decompose the workflow into small steps by specifying how to create sets of output files from sets of input files. Snakemake automatically determines the dependencies between the rules by matching file names.

usage command
```
   $ nano Snakemake
```

### Reads alignment

The first Snakemake rule maps reads of a given sample to a given reference genome. For alignment, I used the tool minimap2. In the Snakefile, I defined the following rule:

Snakemake accepts rule names as targets if the requested rule does not include wildcards, providing flexibility in specifying targets. To manage all desired results efficiently, a rule named **all** can be defined at the top of the workflow, including all typically desired target files as input files.



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

### Converting sam file to bam file

The second rule converts a SAM file to a BAM file . To speeds up processing, saves space, and allows for quick access to specific data regions



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


### Sorting read alignments

### Indexing read alignment 

### Extracting chromosome 6 

### Runing pbhoney for variant calling 

### Runing snifles for variant calling

