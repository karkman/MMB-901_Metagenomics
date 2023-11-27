# MMB-901 Microbial metagenomics – Practicals

__Table of Contents:__

1. [Introduction](#introduction)
2. [Setup](#setup)
    1. [Connnecting to Puhti](#connecting-to-puhti-with-visual-studio-code)
    2. [Setting up](#setting-up-the-course-folders)
    3. [Interactive use of Puhti](#interactive-use-of-puhti)
3. [Data](#data)
4. [Quality control](#quality-control)
5. [Metagenome assembly](#metagenome-assembly)
6. [Read-based taxonomy](#read-based-taxonomy)
7. [Assembly QC](#assembly-qc)
8. [Genome-resolved metagenomics](#genome-resolved-metagenomics)
9. [MAG QC and taxonomy](#mag-qc-and-taxonomy)
10. [MAG annotation](#mag-annotation)
11. [Strain engraftment](#strain-engraftment)

## Introduction

On this course we will use a subset from a publicly available human gut metagenomic sequencing data from fecal microbiota transplantation (FMT) study: [Wilson, B.C., Vatanen, T., Jayasinghe, T.N. _et al._ Strain engraftment competition and functional augmentation in a multi-donor fecal microbiota transplantation trial for obesity. Microbiome 9, 107 (2021)](https://doi.org/10.1186/s40168-021-01060-7).  

Using read-based and assembly-based approaches, we will study the human gut microbiome. With read-based approach we will determine the taxonomic composition of one donor microbiome (female super-donor DF16) and compare to those the taxonomic composition of the recipient microbiomes taken at baseline and at 6, 12 and 26 weeks after the FMT. With assembly based approach, we will use [anvi'o](anvio.org) to construct metagenome-assembled genomes (MAGs) from the donor samples. Then we will assees the strain engrafment of selected strains in a selected set of recipients.  

## Setup

### Connecting to Puhti with Visual Studio Code

* Launch Visual Studio Code
* Down left corner you will have a (green) button with "><" (hoover over it and it says "Open a Remote Window"), click it  
* Choose "Connect Current Window to Host..."
* Type in the __YOUR_USERNAME@puhti.csc.fi__ and hit "Enter" (change "user" for your own CSC username)  
* Type your password and hit "Enter"
* In the following dialogue, type __yes__ and hit "Enter"

When the down left corner says `SSH:puhti.csc.fi`, you're connected.

* From the menu select `Terminal > New Terminal` and you should see a new panel. This is the __command line__.

* When you need to logout just type __exit__ in the terminal/command line and hit "Enter"  
(or you can click the down left corner and choose "Close Remote Connection")

### Setting up the course folders

The main course directory is located at `/scratch/project_2009008`.  
There you will set up your own directory where you will perform all the tasks for this course.  
Some of the tools needed on this course can be found from the course project applications folder `/projappl/project_2009008`.  

First list all projects you're affiliated with in CSC.

```bash
csc-workspaces
```

You should see the course project `MMB-901_metagenomics`.
So let's create a folder for you inside the scratch folder, you can find the path in the output from the previous command.

```bash
cd PATH/TO/COURSE/SCRATCH
mkdir $USER
```

Check with `ls`; which folder did `mkdir $USER` create?

This directory (`/scratch/project_2009008/your-user-name`) is your own working directory.  
Every time you log into Puhti, you should use `cd` to navigate to this directory.

Go to your own folder and clone the course repository there.  

```bash
git clone https://github.com/karkman/MMB-901_Metagenomics.git
```

Check what was downloaded and go to that folder. Then again check what is inside.  
**All the scripts are to be run in this folder** (unless instructed otherwise).  

### Interactive use of Puhti

Puhti uses a scheduling system called SLURM. Most jobs are sent to the queue, but smaller jobs can be run interactively.

Interactive session is launched with `sinteractive`.  
You can specify the resources you need for you interactive work interactively with `sinteractive -i`. Or you can give them as options to `sinteractive`.  
You always need to specify the accounting project (`-A`, `--account`). Otherwise for small jobs you can use the default resources (see below).

| Option | Function | Default | Max |  
| --     | --       | --      | --  |  
| -i, --interactive | Set resources interactively |  |  |  
| -t,  --time | Reservation in minutes or in format d-hh:mm:ss | 24:00:00 | 7-00:00:00 |
| -m, --mem | Memory in Mb       | 2000     | 76000  |  
| -j, --jobname |Job name       | interactive     |   |  
| -c, --cores     | Number of cores       | 1      | 8  |  
| -A, --account     | Accounting project       |       |  |  
| -d, --tmp     | $TMPDIR size (in GiB)      |  32     | 760  |  
| -g, --gpu     | Number of GPUs       | 0     | 0 |  

[__Read more about interactive use of Puhti.__](https://docs.csc.fi/computing/running/interactive-usage/#sinteractive-in-puhti)  

Screen is a handy way to run things in the background without losing them when you logout or have connection problems. However, you have to be careful whn using screen, you can easily get lost.  
And to make things even more complicated, the screen sessions are specific to each login node. And Puhti has at least 4 login nodes.  

**Remember to always first open a screen session and only after that run `sinteractive`.**  

Mini manual for screen:  

* `screen -S NAME` - open a screen and give it a session name NAME
* `screen` - open new screen without specifying any name
* `screen -ls` - list all open sessions
* `ctrl + a + d` - to detach from a session (from inside the screen)
* `screen -r NAME` - re-attach to a detached session using the name
* `screen -rD` - re-attach to a attached session
* `exit` - close the screen and kill all processes running inside the screen (from inside the screen)

## Data

First we will retrieve all the metagenomic sequence files for the female super-donor, DF16.  
Go to the publication and find the sequencing project repository accession (BioProject accession).  
Then go to [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) and find all the reads for the project.  Open them in Run selector and find the run accessions of the samples from the super-donor (9 sequencing runs).  

Finallly download only the accession numbers and make a file (`DF16_accessions.txt`) with only the accession, one line per accesion to the data folder in Puhti (`01_DATA/`).  
We will use that file in the next step to download the data to Puhti.  

We will use [Kingfisher](https://wwood.github.io/kingfisher-download/) to download the read files. It has been installed to Puhti under our course project application folder.  
But before downloading, have a look at the documentation and find out how to use it with a list of accessions. We will use the `ena-ftp` method for downloading the reads.  

Apply for resources. This doesn't take that long, so no real need for screen session.  

```bash
sinteractive -A project_2009008 -m 10G -c 4 
```

Download the data to the data folder.  

```bash
cd 01_DATA

/projappl/project_2009008/Kingfisher/bin/kingfisher get \
    ## fill in the needed options to download the read files
    -m ena-ftp \
    --download-threads $SLURM_CPUS_PER_TASK \
    --extraction-threads $SLURM_CPUS_PER_TASK
```

After the job has finished, check what you got? Make sure you have all 9 sequenicng experiments downloaded.  
If not, check which ones are missing and download them individually with `kingfisher` using the run accession(s).  

## Quality control

Before running any real analyses, we should do quality control (QC) for the sequencing data.  

We use two widely used programs that are pre-installed in Puhti:

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for running the QC fpor each sequence file.
* [MultiqC](https://multiqc.info/) to combine the individual reports from FastQC.  

```bash
cd /scratch/project_2009008/$USER/MMB-901_Metagenomics
mkdir 01_DATA/FASTQC

module load biokit
fastqc --threads $SLURM_CPUS_PER_TASK --outdir 01_DATA/FASTQC 01_DATA/*.fastq.gz

module purge
module load multiqc

multiqc --interactive --outdir 01_DATA/FASTQC 01_DATA/FASTQC/
```

After running the QC steps, download the MultiQC report (`.html`) file to your local computer and open it with browser.  
We'll go thru the report together.  

## Metagenome assembly

As the donor samples are from the same individual, we can do a co-assembly with all the data. For the co-assembly, we'll combine all R1 reads to one file and R2 files to another.  

```bash
cat 01_DATA/SRR*_1.fastq.gz > 01_DATA/DF16_1.fastq.gz
cat 01_DATA/SRR*_2.fastq.gz > 01_DATA/DF16_2.fastq.gz
```

Then the actual assembly will be done with [spoades](https://github.com/ablab/spades) using the `--meta` option meant for metagenomic data. As this will take longer, we run it as batch job.  
Have a look at the batch job file. Using [Puhti](https://docs.csc.fi/computing/running/creating-job-scripts-puhti/) and [spades](https://github.com/ablab/spades#sec3.2) manuals, find out which options do we need to define for the SLURM system and to spades.

```bash
less src/spades.sh
```

Then when you understand what we are about to do, submit the job with `sbatch`.  

```bash
sbatch src/spades.sh
```

The assembly will probably queue for a while and run overnight. You can monitor the queue with:

```bash
squeue -l -u $USER
```

And when it has started running, look at the output log file in `00_LOGS`.  

## Read-based taxonomy

While we wait for the assembly to finish, we can run the read-based taxonomic annotation for the donor samples. And later combine some ready-made output files to compare the recipients to the donor.  
We'll use [metaphlan4](https://github.com/biobakery/MetaPhlAn) for the read-based taxonomic annotation. Metaphlan uses marker genes to profile taxonomic copmposition in metagenomic data.  

To make things run a bit faster, we will run metaphlan as an [array job](https://docs.csc.fi/computing/running/array-jobs/). In a nutshell, each job will be run in parallel as individual jobs. This is a handy way to do the same thing for several files that are independent.
Have a look at the array job file and find out how array jobs are defined by comparing it to the spades batch job we ran earlier.  

```bash
less src/metaphlan.sh
```

After that with the help of [metaphlan4 documentation](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4#basic-usage), figure out the different options we need to define.  

And when you have a basic understanding what we are about to do, submit the job(s).  

```bash
sbatch src/metaphlan.sh
```

And again you can monitor the jobs with the same way as spades batch job.  

Taxonomic profiling doesn't take that long, approx 30 min per sample, but as they run in parallel, it will not take n x 30 min, but a lot less depending on the queue.  
While you wait, you can log in to Puhti web interface at [www.puhti.csc.fi](http://www.puhti.csc.fi) and open a Rstudio session (use the default options).  
We will analyse the results in R using few packages for microbiome data analysis.  

But before we can read in the data to R, we need to combine the individual metaphlan outputs and extract the species level annotations from there.  

```bash
module load metaphlan/4.0.6
merge_metaphlan_tables.py 05_TAXONOMY/SRR*.txt > 05_TAXONOMY/metaphlan.txt
awk '$1 ~ "clade_name" || $1 ~ "s__" {print $0}' 05_TAXONOMY/metaphlan.txt |grep -v "t__" > 05_TAXONOMY/metaphlan_species.txt
```

Then you can follow the R instruction in the file `src/taxonomic_profiling.r` and run the analysis in browser interface of Rstudio running at Puhti.  

After we have analysed the taxonomic profiles of the donor, we can combine the rest of the samples to our merged metaphlan table and run the analysis again.  
First copy the taxonomic profiles of additional 192 samples to the metaphlan output folder and re-run the merge command above.  

```bash
cp /scratch/project_2009008/Data/metaphlan/*.txt 05_TAXONOMY/

merge_metaphlan_tables.py 05_TAXONOMY/SRR*.txt > 05_TAXONOMY/metaphlan.txt
awk '$1 ~ "clade_name" || $1 ~ "s__" {print $0}' 05_TAXONOMY/metaphlan.txt |grep -v "t__" > 05_TAXONOMY/metaphlan_species.txt
```

Then re-run the R part.  

## Assembly QC

## Genome-resolved metagenomics

## MAG QC and taxonomy

```bash
export GTDBTK_DATA_PATH="/scratch/project_2009008/DB/release214/"
```

## MAG annotation

## Strain engraftment

## Automatic binning
