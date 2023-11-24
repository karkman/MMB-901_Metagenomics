# MMB-901 Microbial metagenomics – Practicals

## Introduction

On this course we will use a subset from a publicly available human gut metagenomic sequencing data from fecal microbiota transplantation (FMT) study: [Wilson, B.C., Vatanen, T., Jayasinghe, T.N. _et al._ Strain engraftment competition and functional augmentation in a multi-donor fecal microbiota transplantation trial for obesity. Microbiome 9, 107 (2021)](https://doi.org/10.1186/s40168-021-01060-7).  

Using read-based and assembly-based approaches, we will study the human gut microbiome. With read-based approach we will determine the taxonomic composition of one donor microbiome (female super-donor DF16) and compare to those the taxonomic composition of the recipient microbiomes taken at baseline and at 6, 12 and 26 weeks after the FMT. With assembly based approach, we will use [anvi'o](anvio.org) to construct metagenome-assembled genomes (MAGs) from the donor samples. Then we will assees the strain engrafment of selected strains in a selected set of recipients.  

## Setup

### Connecting to Puhti with Visual Studio Code

* Launch Visual Studio Code
* Down left corner you will have a (green) button with "><" (hoover over it and it says "Open a Remote Window"), click it 
* Choose "Connect Current Window to Host..."
* Type in the **user<span>@puhti.csc.fi** and hit "Enter" (change "user" for your own CSC username) 
* Type your password and hit "Enter"
* In the following dialogue, type **yes** and hit "Enter"

When the down left corner says `SSH:puhti.csc.fi`, you're connected.

* From the menu select `Terminal > New Terminal` and you should see a new panel. This is the __command line__.

* When you need to logout just type **exit** in the terminal/command line and hit "Enter"  
(or you can click the down left corner and choose "Close Remote Connection")

## Setting up the course folders

The main course directory is located at `/scratch/project_2009008`.  
There you will set up your own directory where you will perform all the tasks for this course.  
Some of the tools needed on this course can be found from the course project applications folder `/projappl/project_2009008`.  

First list all projects you're affiliated with in CSC.

```
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
**All the scripts are to be run in this folder**.  

## Interactive use of Puhti

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

[**Read more about interactive use of Puhti.**](https://docs.csc.fi/computing/running/interactive-usage/#sinteractive-in-puhti)  

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

First we will retrieve all the metagenomic sequence files for the female super donor, DF16. Go to the publication and find the sequencing project repository accession (BioProject accession).  
Then go to [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) and find all the reads for the project.  Open them in Run selector and find the run accessions of the samples from the super-donor (9 sequencing runs).  
Finallly download only the accession numbers. We will use that file in the next step to download the data to Puhti.  

We will use [Kingfisher](https://wwood.github.io/kingfisher-download/) to download the read files. It has been installed to Puhti under our course project application folder.  
But before downloading, have a look at the documentation and find out how to use it with a list of accessions. We will use the `ena-ftp` method for downloading the reads.  

Apply for resources. Using a screen session or without. This doesn't take that long, so no real need for screen session.  

```bash
sinteractive -A project_2009008 -m 10G -c 4 
```

Download the data to the data folder.  

```bash
cd 01_DATA

~/projappl/Kingfisher/bin/kingfisher get \
    ## fill in the needed options to download the read files
    -m ena-ftp \
    --download-threads $SLURM_CPUS_PER_TASK \
    --extraction-threads $SLURM_CPUS_PER_TASK
```

After the job has finished, check what you got? Make sure you have all 9 sequenicng experiments downloaded.  
If not, check which ones are missing and download them individually with `kingfisher` using the run accession(s).  

## Quality control

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

## Metagenome assembly

```bash
cd 01_DATA
cat SRR*_1.fastq.gz > DF16_1.fastq.gz
cat SRR*_2.fastq.gz > DF16_2.fastq.gz
```

```bash
sbatch src/spades.sh
```

## Read-based taxonomy

## Assembly QC

## Genome-resolved metagenomics

## MAG QC

## Strain engraftment

