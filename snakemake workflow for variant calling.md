
## Snakemake Workflow

Snakemake is a workflow management system that aims to make the definition and execution of data analysis workflows easy, reproducible, and scalable. It uses a Python-based language to specify workflows in terms of rules that describe how to derive output files from input files.
For a comprehensive introduction and installation of Snakemake, refer to the official tutorial: (https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).


## Environment setup

## Sample data

1. Create a directory for the task

   $ mkdir Snakemak

1. Download the new dataset from this link https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9205427/#sec4.3.2.
1. Synchronize the input files using `rsync`: To do this between your local machine and the remote server use the following command:



## Snakemake workflows 

1. A Snakemake workflow is defined by specifying rules in a Snakefile. To begin, create a new file called Snakefile in your working directory using any text editor you prefer. Rules decompose the workflow into small steps by specifying how to create sets of output files from sets of input files. Snakemake automatically determines the dependencies between the rules by matching file names.

usage command
```
   $ nano Snakemake
```

### Reads alignment

The first Snakemake rule maps reads of a given sample to a given reference genome. For alignment, I used the tool minimap2. In the Snakefile, I defined the following rule:

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

### Sorting read alignments

### Indexing read alignment 

### Extracting chromosome 6 

### Runing pbhoney for variant calling 

### Runing snifles for variant calling

