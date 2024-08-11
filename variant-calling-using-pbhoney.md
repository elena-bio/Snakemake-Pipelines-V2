
# Structural Variant Calling with PBHoney and Sniffles

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

1. **Mambaforge:** ~~Information how to install `conda` and add the `bioconda` channel is available on https://bioconda.github.io/.~~ :bookmark: [TODO] add a description for Mambaforge
   - `#green` `#ff5500` :bookmark: [TODO] add the command here

1. **PBHoney** :bookmark: [TODO] Did I download & install it as a standalone package, or did I include it as a Conda dependency?!

1. :bookmark: [TODO] add other tools, if any


## Environment setup

1. **Set up  dedicated Conda environments for PBHoney and Sniffles**: Creating a separate Conda environments ensures that the dependencies and packages required for PBHoney and sniffles do not conflict with other tools or projects on the system

   ```
   $ conda create -n pbhoney-env python=2.7.15
   ```
   ```
   $ conda create -n sniffles_env
   ```
 
1. **Activate the environment**
   
usage command:

   `$ conda activate pbhoney-env`

1. 1. **Create the environment.yaml file**: Specify the configuration of a Conda environment, including packages and dependencies.

- In the working directory, create new files called `sniffles_env.yaml` and pbhoney_env.yaml with an editor of your choice.

```
gedit sniffles_env.yaml
```
```
gedit pbhoney_env.yaml
```


   ```yaml
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
   
1. I ran the code again after fixing the issue, but I got this error : Can not convert non-pacbio reads to pbbam record.

   ![Image 2024-06-16 at 23 18 48](https://github.com/elena-bio/snakemake-pipelines-v2/assets/166811346/7d722ca1-b924-4639-98f5-1fbcfe1f60ee)


## Generating a VCF File from a New Input Dataset using PBHoney

1. **Reads alignment with `Minimap2` tool**

   `Minimap2` is a versatile alignment tool designed for mapping DNA or RNA sequence reads to a reference genome or for aligning long reads to each other.
   
   `Minimap2` is publicly available and open source.Minimap2 seamlessly works with gzip’d FASTA and FASTQ formats as input, no need to convert between FASTA and FASTQ or decompress gzip’d files first.

   `$ minimap2 -ax map-pb reference.fa pacbio_reads.fq > pb_alignment.sam`

   In this command example, `-ax` specifies the preset for different types of reads, and `-x `specifies the preset for the overlap layout. The input and output files are specified accordingly, with the output being directed to SAM or PAF files.

   #### Usagge of the command

   `$ minimap2 -ax map-hifi GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta HG002_35x_PacBio_14kb-15kb.fastq.gz > output.sam`

   The current version used is `Minimap2 2.28-r1209`.

   To run jobs in the HPC environment, it is necessary to submit jobs to a queueing system. The HPC environment queueing system shares compute resources fairly and efficiently among users and schedules when jobs should run. I have included a job submission script below. This script is crucial for the effective execution of tasks within the HPC environment.

   ```
   #!/bin/bash
   ### Account information
   #PBS -W group_list=cu_10160 -A cu_10160
   ### Job name
   #PBS -N minimap2
   ###Output files
   #PBS -e minimap2.err
   #PBS -o minimap2.log
   ### Only send mail when job is aborted or terminates abnormally
   #PBS -m n
   ### Number of nodes
   #PBS -l nodes=1:ppn=8
   ### Memory
   #PBS -l mem=40gb
   ###Requesting time
   #PBS -l walltime=12:00:00

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
   conda activate pbhoney-env

   # Execute your command (e.g., minimap2 command)
   minimap2 -ax map-hifi GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta HG002_35x_PacBio_14kb-15kb.fastq.gz > output.sam
   
   ```

   $ qsub -W group_list=<group_name_found_in_step#1> -A <group_name_found_in_step#1> -l nodes=1:ppn=1,mem=4gb, walltime=01:00:00 <your script>`
   
   usage command:
   ```
   $ qsub -W group_list=cu_10160 -A cu_10160 -l nodes=1:ppn=8,mem=40gb,walltime=12:00:00 ./run_minimap2.sh

   ```
   to check the status :
   ```
   $ qstat -a; ll -h

   ```
  


1. **Converting sam file to bam file**
 The second rule converts a SAM file to a BAM file . To speeds up processing, saves space, and allows for quick access to specific data regions.

usage command in bash script:

```
# Execute your command (e.g., converting sam to bam command)
samtools view -Sb output.sam > output.bam
```
I also chenged these comands in the bash script: 
```
###Output files
#PBS -e sam_bam.err
#PBS -o sam_bam.log
```
1. **Sorting read alignments**
   
For later steps, we need the read alignments in the BAM files to be sorted. This can be achieved with the **samtools** `sort` command.

usage command in bash script:
```
# Execute your command (e.g., sorting sam to bam command)
samtools sort output.bam -o sorted.bam
```
i also changed these commands in the bash script:
```
###Output files
#PBS -e sort_bam.err
#PBS -o sort_bam.log
```
### Indexing read alignment 

Next, we need to use `samtools` again to index the sorted read alignments so that we can quickly access reads by the genomic location they were mapped to. This can be done with the following command:

```
# Execute your command (e.g., indexing sorted bam file command)
samtools index sorted.bam
```
1. Extract chromosome number 6
   
 usage command in bash script:   
```
# Execute your command (e.g., extracting chromosome6)
samtools view -b sorted.bam chr6 > chromosome6.bam

```
I also changed these command in the bash script:

```
#PBS -e extract_chr_bam.err
#PBS -o extract_chr_bam.log
```

### Runing Sniffles Tool for Variant Calling
   
   __Sniffles2__: A rapid and accurate structural variant caller designed for long-read sequencing. Sniffles2 efficiently detects structural variants (SVs) across germline, somatic, and population-level studies using data from PacBio and Oxford Nanopore technologies.
    
    **Installation**
    
- ● Use `pip` with the command: `pip install sniffles`
- ● Use `conda` with the command:`conda install sniffles=2.4`

If Sniffles1 is already installed via `conda`, it can be upgraded to Sniffles2 with the following command:

- `conda update sniffles=2.4`
    
    **Requirements**
    
   ```
name: sniffles_env
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.10.12
  - pysam=0.21.0
  - edlib>=1.3.9
  - psutil>=5.9.4
  - pip  # To ensure pip is available for any additional installations
  - pip:
    - sniffles2-plot
 <!-- minimap2 for allignment
  - minimap2
  - samtools=1*
  ``` 
### Runing PBHoney Tool for Variant Calling 
