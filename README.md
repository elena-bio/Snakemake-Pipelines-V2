
# Structural Variant Calling with PBHoney

_PBHoney_ is a tool used for detecting structural variants in DNA, helping researchers identify large-scale changes like deletions or insertions in genomes accurately and efficiently. It analyzes long-read sequencing data to pinpoint genetic differences that can impact health and disease studies.
This guideline utilizes publicly available long-read sequencing data. Initially, I tried to use a prepared BAM file but encountered errors, as explained in the following steps. Instead, I used a new dataset for example PacBio HiFi long-read sequencing data from son of Ashkenazi Jewish ancestry (HG002). All steps and tools used have been included.


## Connecting to Computerome Server 

There are 3 ways to login to Computerome's HPC environment: 
1. SSH via terminal
1. SSH via client
1. Virtual desktop
   
The following explanation is based on SSH via Terminal as a login service. 

1. On a Mac or Linux machine open a new Terminal and execute the following command (replace `username` with the one you received from Computerome). 
   ```
   $ ssh username@ssh.computerome.dk
   ```
1. Next, enter **the temporary password** you received via SMS. 
1. Now you will receive an Entrust Identity mobile notification to confirm. 
1. After you have confirmed the notification you will be logged into the HPC environment.

The disk quota of my personal home directory (`/home/people/user`) was 10 GB. But luckily Computerome provided me with another directory with much larger disk quota. Therefore, I changed the current working directory using this command:
```
$cd /home/projects/cu_10160/people/elehos
```


## Required tools

1. **Mambaforge:** ~~Information how to install `conda` and add the `bioconda` channel is available on https://bioconda.github.io/.~~ [!TODO] add a description for Mambaforge
   - [!TODO] add the command here

1. **PBHoney** [!TODO] Did I download & install it as a standalone package, or did I include it as a Conda dependency?!

1. [!TODO] add other tools, if any


## Environment setup

1. **Set up a dedicated Conda environment for PBHoney**: Creating a separate Conda environment ensures that the dependencies and packages required for PBHoney do not conflict with other tools or projects on the system

   ```
   $ conda create -n pbhoney-env python=2.7.15
   ```

1. **Activate the environment**

   `$ conda activate pbhoney-env`

1. **Create the environment.yaml file**: Specify the configuration of a Conda environment, including packages and dependencies.
   ```
   channels:
   - conda-forge
   - bioconda
   dependencies:
   - python =2.7
   - biopython
   - samtools =0.1.*
   - blasr
   - minimap2  
   - pysam =0.8
   - networkx
   - h5py =2.* 
   - numpy =1.*
   - pbdagcon
   - picard =2.7.*
   - intervaltree_bio
   ```

1. **Update the conda environment**: Update the current environment based on environment file.

   `$ conda env update -f environmentyaml`


## Sample data

1. Create a directory for the task 

   `$ mkdir dna_lr_seq`

1. ~~Checking directory creation with `ls` command~~

1. Download the new dataset from this link https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9205427/#sec4.3.2.

1. Synchronize the input files using `rsync`: To do this between your local machine and the remote server use the following command:

   `$ rsync "/path/to/local/directory/" user@transfer.computerome.dk:/path/to/remote/directory/`

The task will run with the following inputs

* prepared bam file and refrence genome fasta file 
* HG002_35x_PacBio_14kb-15kb.fastq.gz and refrence genome fasta file


## Figure out if PBHoney can generate a VCF file from the prepared bam file

1. Trying to run pie stage of PBHoney on my first sample BAM file

   `$ Honey pie <Bamfile> <refrence genome fastafile>`

   * I encountered this error: `IndexError: list index out of range`

1. Upon review, I have decided to address and fix the issue directly in the source code and I have added a `continue` statement to the loop to bypass instances of bad data with an empty CIGAR string.

   ![alt text](<CIGAR-error.jpeg>)


## Generating a VCF File from a New Input Dataset using PBHoney

1. Alligning PacBio reads with `Minimap2` tool 

   `Minimap2` is a versatile alignment tool designed for mapping DNA or RNA sequence reads to a reference genome or for aligning long reads to each other.
   
   `Minimap2` is publicly available and open source.Minimap2 seamlessly works with gzip’d FASTA and FASTQ formats as input, no need to convert between FASTA and FASTQ or decompress gzip’d files first.

   `$ minimap2 -ax map-pb reference.fa pacbio_reads.fq > pb_alignment.sam`

   In this command example, `-ax` specifies the preset for different types of reads, and `-x `specifies the preset for the overlap layout. The input and output files are specified accordingly, with the output being directed to SAM or PAF files.

   #### Usagge of the command

   `$ minimap2 -ax map-hifi GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta HG002_35x_PacBio_14kb-15kb.fastq.gz > output.sam`

   The current version used is `Minimap2 2.28-r1209`.

   To run jobs in the HPC environment, it is necessary to submit jobs to a queueing system.

1. Extract chromosome number 6

1. Runing PBHoney 

1. Runing PBHoney Tool for Variant Calling 
