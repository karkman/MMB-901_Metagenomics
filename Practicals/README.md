# MMB-901 Microbial metagenomics – Practicals

__Table of Contents:__

1. [Introduction](#introduction)
2. [Setup](#setup)
    1. [Connnecting to Roihu](#connecting-to-Roihu-with-visual-studio-code)
    2. [Setting up](#setting-up-the-course-folders)
    3. [Interactive use of Roihu](#interactive-use-of-Roihu)
3. [Data](#data)
4. [Quality control](#quality-control)
5. [Metagenome assembly](#metagenome-assembly)
6. [Read-based taxonomy](#read-based-taxonomy)
7. [Assembly QC](#assembly-qc)
8. [Genome-resolved metagenomics](#genome-resolved-metagenomics)
9. [MAG QC and taxonomy](#mag-qc-and-taxonomy)
10. [MAG annotation](#mag-annotation)
11. [Strain engraftment](#strain-engraftment)
12. [Automatic binning](#automatic-binning)

## Introduction

On this course we will use a subset from a publicly available human gut metagenomic sequencing data from fecal microbiota transplantation (FMT) study: [Wilson, B.C., Vatanen, T., Jayasinghe, T.N. _et al._ Strain engraftment competition and functional augmentation in a multi-donor fecal microbiota transplantation trial for obesity. Microbiome 9, 107 (2021)](https://doi.org/10.1186/s40168-021-01060-7).  
There is also a 4-year follow-up study, if you're interested: [Long-term health outcomes in adolescents with obesity treated with faecal microbiota transplantation: 4-year follow-up](https://www.nature.com/articles/s41467-025-62752-4).  

Using read-based and assembly-based approaches, we will study the human gut microbiome. With read-based approach we will determine the taxonomic composition of one donor microbiome (female super-donor DF16) and compare to those the taxonomic composition of the recipient microbiomes taken at baseline and at 6, 12 and 26 weeks after the FMT. With assembly based approach, we will use [anvi'o](anvio.org) to construct metagenome-assembled genomes (MAGs) from the donor samples. Then we will assees the strain engrafment of selected strains in a selected set of recipients.  

## Connecting to Roihu

If you not experienced with UNIX and command line, the easiest way is to use the Roihu web interface, [www.roihu.csc.fi](www.roihu.csc.fi). You can log in with Haka or CSC crecentials.  
If you wish to use your own computer, you can connect to Roihu using SSH. For that follow the instructions below.  

### SSH key

You need to set up SSH key before you can connect to CSC from your local machine using SSH.  
If you're not sure whether you have done it, you can check it from [my.csc.fi](my.csc.fi) under your profile and SSH public keys.  

In case you haven't set up SSH key, follow the instructions at [CSC](https://docs.csc.fi/support/tutorials/roihu/#ssh-client).  

### SSH certificate

In addition to the SSH key, you need a SSH certificate that is valid for 24 h. You can follow the instruction at CSC [Roihu tutorial.](https://docs.csc.fi/support/tutorials/roihu/)  

### SSH connection

* Launch Terminal
(e.g. open the Launchpad and type __terminal__)
* Type `ssh USER@roihu-cpu.csc.fi` and hit "Enter" (change `USER` to your own CSC username)
* In the following dialogue, type __yes__ and hit "Enter"
* Type your password and hit "Enter"
* To logout first type __exit__, hit "Enter", and close the Terminal window

## Setup

### Setting up the course folders

The main course directory is located at `/scratch/project_2022405`.  
There you will set up your own directory where you will perform all the tasks for this course.  
Some of the tools needed on this course can be found from the course project applications folder `/projappl/project_2022405`.  

First list all projects you're affiliated with in CSC.

```bash
csc-workspaces
```

You should see the course project `MMB-901_Metagenomics`.
So let's create a folder for you inside the scratch folder, you can find the path in the output from the previous command.

```bash
cd PATH/TO/COURSE/SCRATCH
mkdir $USER
```

Check with `ls`; which folder did `mkdir $USER` create?

This directory (`/scratch/project_2022405/your-user-name`) is your own working directory.  
Every time you log into Roihu, you should use `cd` to navigate to this directory.

Go to your own folder and clone the course repository there.  

```bash
git clone https://github.com/karkman/MMB-901_Metagenomics.git
```

Check what was downloaded and go to that folder. Then again check what is inside.  
__All the scripts are to be run in this folder__ (unless instructed otherwise).  

### Interactive use of Roihu

Roihu uses a scheduling system called SLURM. Most jobs are sent to the queue, but smaller jobs can be run interactively.

Interactive session is launched with `sinteractive`.  
You can specify the resources you need for you interactive work interactively with `sinteractive -i`. Or you can give them as options to `sinteractive`.  
You always need to specify the accounting project (`-A`, `--account`). Otherwise for small jobs you can use the default resources (see below).

| Option | Function | Default | Max |  
| --     | --       | --      | --  |  
| -i, --interactive | Set resources interactively |  |  |  
| -t,  --time | Reservation in minutes or in format d-hh:mm:ss | 08:00:00 | 36:00:00 |
| -m, --mem | Memory in Mb       | 2G     | 64G |  
| -j, --jobname |Job name       | interactive     |   |  
| -c, --cores     | Number of cores       | 2      | 32  |  
| -A, --account     | Accounting project       |       |  |  
| -d, --tmp     | $TMPDIR size      |  0     | 2T  |  
| -g, --gpu     | Number of GPUs       | 0     | 0 |  

[__Read more about interactive use of Roihu.__](https://docs.csc.fi/computing/running/interactive-usage/#the-sinteractive-command)  

[Tmux](https://docs.csc.fi/support/tutorials/tmux/index.html) is a handy way to run things in the background without losing them when you logout or have connection problems. However, you have to be careful whn using tmux, you can easily get lost.  
And to make things even more complicated, the sessions are specific to each login node. And Roihu has at least 4 login nodes.  

__Remember to always first open a screen session and only after that run `sinteractive`.__  

Mini manual for Tmux:  

* `tmux new -s NAME` - open a new session and give it name NAME
* `tmux ls` - list all open sessions
* `ctrl + b + d` - to detach from a session (from inside the session)
* `tmux attach -t NAME` - re-attach to a detached session using the name
* `ctrl + d` or `exit` - close the session and kill all processes running inside it (from inside the session)
* `ctrl + b + [` - scroll mode in the session (use arrow keys to scroll), otherwise you cannot scroll in the session  

## Data

First we will retrieve all the metagenomic sequence files for the female super-donor, __DF16__.  
Go to the publication and find the sequencing project repository accession (BioProject accession).  
Then go to [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) and find all the reads for the project.  Open them in Run selector and find the run accessions of the samples from the super-donor (9 sequencing runs).  

Finallly download only the accession numbers and make a file (`DF16_accessions.txt`) with only the accession, one line per accesion, to the data folder in Roihu (`01_DATA/`).  
We will use that file in the next step to download the data to Roihu.  

We will use a tool called [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) from NCBI to download the sequencing read files. It is installed on Roihu by default under the biokit module, so we need to load that first.  
But before downloading, have a look at the documentation to understand the different options we are using.  

Apply for resources. This doesn't take that long, so no real need for screen session.  

```bash
sinteractive -A project_2022405 -m 10G -c 8 -t 03:00:00
```

Then navigate to your course folder and run the following commands to download the data.

As we have a list of accessions, we could either download each one by one, write them separately on the command line, or we can make a for loop that reads the accession file and downloads them one by one.  
Or if you know how [parallel](https://www.gnu.org/software/parallel/parallel_tutorial.html) works or want to ask AI, you can make a parallel version of the for-loop. Just load also parallel module.  

```bash
export TMPDIR=.

module load biokit

mkdir -p 01_DATA

for ACC in `cat 01_DATA/DF16_accessions.txt`; do
    fasterq-dump \
        --split-3 \
        --skip-technical \
        --outdir 01_DATA \
        --threads $SLURM_CPUS_PER_TASK \
        --progress \
        $ACC
done
```

After the job has finished, check what you got? Make sure you have all 9 sequencing experiments downloaded.  
If not, check which ones are missing and download them individually with `fasterq-dump` using the run accession(s).  

To save some space, we should compress the sequence files (`fasterq-dump` downloads uncompressed fastq files).  
This is rather slow, so consider running it with tmux. Before opening a tmux session, check the commands from above.  
Also remember to allocate some resources. More CPUs will make the compression faster. We will use [pigz](https://zlib.net/pigz/) for the compression, which is a parallel implementation of gzip.

```bash
pigz -p $SLURM_CPUS_PER_TASK 01_DATA/*.fastq
```

After this is done, check that each sequence file in the folder `01_DATA` has the ending `.gz`.  

In case `fasterq-dump` gave errors and some of the files are missing, you can copy them from the course data folder (`/scratch/project_2022405/GutBugsData/Rawdata`).  

## Quality control

Before running any real analyses, we should do quality control (QC) for the sequencing data.  

We use two widely used programs that are pre-installed in Roihu:

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for running the QC fpor each sequence file.
* [MultiQC](https://multiqc.info/) to combine the individual reports from FastQC.  

```bash
cd /scratch/project_2022405/$USER/MMB-901_Metagenomics
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

Then the actual assembly will be done with [spades](https://github.com/ablab/spades) using the `--meta` option meant for metagenomic data. As this will take longer, we'll run it as batch job.  

Using AI, prepare a batch job script for running metagenomic assembly based on this information. Ask for a detailed description of the file content and the different options, so you understand what is going on.  
You will need to request 12 CPUs, 150G of memory and 500G of local scratch space. The maximum time for the job is 16 hours.
Roihu uses a scheduling system called SLURM.
Make sure the files names, folders, options and project number are correct.  
Save the script as `src/spades.sh`.  

You can also use [Roihu](https://docs.csc.fi/computing/running/creating-job-scripts-roihu/) and [spades](https://github.com/ablab/spades#sec3.2) manuals, to learn more.  

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

Next we run the read-based taxonomic annotation for the donor samples. And later combine some ready-made output files to compare the recipients to the donor.  
We'll use [metaphlan4](https://github.com/biobakery/MetaPhlAn) for the read-based taxonomic annotation. Metaphlan uses marker genes to profile taxonomic compposition in metagenomic data.  
You can read about the basic usage from here: [metaphlan4 documentation](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4.2#basic-usage).  

To make things run a bit faster, we will run metaphlan as an [array job](https://docs.csc.fi/computing/running/array-jobs/). In a nutshell, all jobs will be run in parallel as individual jobs. This is a handy way to do the same thing for several files that are independent.  

Use CSC documentation and metaphlan manual (and/or AI) to prepare a batch job script for running metaphlan4 as an array job based on this information. If using AI, ask for a detailed description of the file content and the different options, so you understand what is going on.  
The You will need to request 8 CPUs, 50G of memory and 100G of local scratch space. The maximum time for the job is 1 hour.  
The version of metaphlan we'll use is 4.2.4 and it is installed on Roihu under the module `py-metaphlan/4.2.4` (needs also `bio-apps` module).

Make sure the files names, folders, options and project number are correct and save the script as `src/metaphlan.sh`.  

And when you have a basic understanding what we are about to do, submit the job(s).  

```bash
sbatch src/metaphlan.sh
```

And again you can monitor the jobs with the same way as spades batch job.  

Taxonomic profiling doesn't take that long, approx 30 min per sample, but as they run in parallel, it will not take n x 30 min, but a lot less depending on the queue.  
While you wait, you can log in to Roihu web interface at [www.Roihu.csc.fi](http://www.Roihu.csc.fi) and open a Rstudio session (use the default options).  
We will analyse the results in R using few packages for microbiome data analysis.  

When all the metaphlan jobs are finished, we can go on with the data analysis of the microbial community.  
But before we can read the data into R, we need to combine the individual metaphlan outputs and extract the species level annotations from there.  

```bash
module load metaphlan/4.2.4
merge_metaphlan_tables.py 03_TAXONOMY/SRR*.txt > 03_TAXONOMY/metaphlan.txt
```

Then follow the R instruction in the file `src/taxonomic_profiling.r` and run the analysis in Roihu using the browser interface of Rstudio.  

After we have analysed the taxonomic profiles of the donor, we combine the rest of the samples to our merged metaphlan table and run the analysis again.  
First copy the taxonomic profiles of additional ~192 samples to the metaphlan output folder and re-run the merge command above.  

```bash
cp /scratch/project_2022405/GutBugsData/metaphlan/*.txt 03_TAXONOMY/

merge_metaphlan_tables.py 03_TAXONOMY/SRR*.txt > 03_TAXONOMY/metaphlan.txt
```

Then re-run the R part.  

## Assembly QC

Although it is not straightforward to assess the quality of a metagenomic assembly, we can still run a QC ananlysis and at least see what we got.  
We will use the metagenomic version of [QUAST](https://quast.sourceforge.net/docs/manual.html) for the job.  

QUAST can be found from Roihu, so just need to load the quast module.  
But first allocate some resources.  

```bash
sinteractive --account project_2022405
```

While you wait for the resources, have a look at the quast manual and read about the options were using.  
Then when you're connected to a computing node and have read about the options from the manual, run quast.

```bash
module load quast/5.2.0 
metaquast.py 02_ASSEMBLY/contigs.fasta --max-ref-num 0 --threads $SLURM_CPUS_PER_TASK -o 02_ASSEMBLY/QUAST --fast
```

When the job is finished, free the resources with:  

```bash
exit
```

And have a look at the output folder and inspect the results. QUAST manual will help you with the different output files.  

## Genome-resolved metagenomics

Instead of manual binning, we will use automatic binning software to reconstruct some of the original genomes (metagenome-assembled genomes, MAGs) in our donor samples.  
There are several different tools for automated binning and they all perform very differently (as you probably read in the blog post).  
The programs we will use are [Metabat](LINK) and [Semibin2](https://github.com/BigDataBiology/SemiBin). Both use the differential coverage information and the sequence composition of the contigs as the basis to cluster them into bins.  
There's also a [´omics vocabulary](https://anvio.org/vocabulary/), that can help you to understand some key terms and concepts in microbial 'omics.  

### Mapping

To obtain the differential coverage information for each contig in our database, we will map the original reads from each sample individually to the assembly.  
We use [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for the mapping and the first step is to create an index of the fasta file.  

```bash
mkdir -p 04_MAPPING
bowtie2-build 02_ASSEMBLY/contigs.fasta 04_MAPPING/contigs --threads $SLURM_CPUS_PER_TASK
```

As we will be running the same job for all nine samples, the quickest way is to create an array job for the task.  
Copy the `metaphlan.sh` script to a new file called `mapping.sh` in the same `src` folder and modify it accordingly:

* Change the name of the job and the names of the log files
* Change the time to 1 hour
* Change the memory to 12G
* Instead of metaphlan, load bowtie2 and samtools modules

Remember that each of the array jobs reads the `DF16_accessions.txt` file and picks the right sample to analyse.  
In more detail; each array has its own ID from 1–9 stored in the environmental variable `SLURM_ARRAY_TASK_ID`. Each arrays reads the corresponding line from the file and stores it to the variable `SAMPLE_ACC`.  
Make sure the path to the accession file is correct in the following line: `SAMPLE_ACC=$(sed -n ${SLURM_ARRAY_TASK_ID}p 01_DATA/DF16_accessions.txt)`  

__Do not run__ the following, but these are the mapping and profiling commands that you need to change in the file:  

```bash
bowtie2 \
    -x 04_MAPPING/contigs \
    -1 01_DATA/${SAMPLE_ACC}_1.fastq.gz \
    -2 01_DATA/${SAMPLE_ACC}_2.fastq.gz \
    --no-unal \
    -p $SLURM_CPUS_PER_TASK \
    -S $TMPDIR/${SAMPLE_ACC}.sam

samtools view -bS $TMPDIR/${SAMPLE_ACC}.sam > $TMPDIR/${SAMPLE_ACC}.bam
samtools sort -o 04_MAPPING/${SAMPLE_ACC}.sorted.bam $TMPDIR/${SAMPLE_ACC}.bam
samtools index 04_MAPPING/${SAMPLE_ACC}.sorted.bam
```

After the script is ready, you can submit the job(s).  

```bash
sbatch src/mapping.sh
```

### Metabat

INSTRUCTIONS HERE

```bash
05_BINNING/METABAT2/
```

### Semibin2


Semibin2 uses deep learning in metagenomic binning and has pre-trained models for several different environments, including human gut. However, with multiple samples the pre-trained models cannot be used and we would need to train the model with our data. The model training takes some time, although it might be more accurate. Running Semibin2 with pre-trained model is a lot faster and we will prefer speed over accuracy on this course. In case we have time, we can try the multi-sample mode for our data.  

We will use the same files we created and used in the manual binning with anvi'o. The job can be run interactively or as a batch job. You will need at least 6 CPUs, 50G of memory, 100G of local storage and many hours (to be sure this time).  
Allocate the resources or write a batch job script and use the following command to run Semibin2.  

Also make sure to run this from our course main folder (MMB-901_Metagenomics) so that all the paths are right.  

```bash
module load bio-apps/v202603
module load py-semibin/2.6.0

SemiBin2 single_easy_bin \
    --input-fasta 02_ASSEMBLY/contigs.fasta \
    --input-bam 04_MAPPING/*.bam \
    --tmpdir $TMPDIR \
    --threads $SLURM_CPUS_PER_TASK \
    --output 05_BINNING/SEMIBIN2/
```

The outputs from Semibin2 will be written to `05_BINNING/SEMIBIN2/` and the folder will contain information about the bins, some log files and one folder with all the genome bins as compressed fasta files.  

NOT SURE ABOUT THIS, COULD SOMEONE CHECK?  
Semibin2 produces a file (`contig_bins.tsv`) that has all the needed information, but the format is not exactly as anvi'o would wish.
So we need to do some command line magic. Or in other words use `awk` to remove the first line and add a prefix to each bin name (as anvi'o does not accept only numbers as bin names).

```bash
awk 'NR!=1 {print $1"\tSemibin2_"$2}' 08_AUTOMATED_BINNING/contig_bins.tsv > 08_AUTOMATED_BINNING/semibin_collection.txt
```

## MAG QC and taxonomy

When the binning is done, it's time to assess the quality of the bins and give a taxonomic annotation to each bin.  
We will use [CheckM2](https://github.com/chklovski/CheckM2) for assessing the quality of our MAGs and [GTDB-Tk](https://github.com/ecogenomics/gtdbtk) for the taxonomic annotation.  

The MAG quality control and taxonomic annotation will be run as a batch job to be more efficient and just for practice. Run them as separate batch jobs and start with the one that probably takes longer.  
Below you can find the scripts. Resources for the batch jobs are as follows:

__CheckM2:__

* 20G of memory
* 1 hour of time
* 200G of local storage
* 10 CPUs

__GTDB-Tk:__

* 120G of memory
* 2 hours of time
* 200G of local storage
* 10 CPUs

You can use the spades script as a template.

Remember to analyse the bins from each binning method separately and make sure you have the extensions correct for each method (for example wth Semibin2 it's `.fa.gz`).  

### MAG QC with CheckM2

```bash
module load  bio-apps/v202603
module load checkm2/1.1.0

checkm2 predict \
    --input BIN_FOLDER \
    --output-directory BIN_FOLDER/checkm2 \
    --extension fa \
    --tmpdir $TMPDIR \
    --threads $SLURM_CPUS_PER_TASK
```

### MAG taxonomy with GTDB-Tk

```bash
module load  bio-apps/v202603
module load py-gtdbtk/2.7.2

gtdbtk classify_wf \
    --genome_dir BIN_FOLDER \
    --out_dir BIN_FOLDER/gtdbtk \
    --extension fa \
    --skip_ani_screen \
    --cpus $SLURM_CPUS_PER_TASK \
    --pplacer_cpus 1 \
    --scratch_dir $TMPDIR \
    --tmpdir $TMPDIR
```

## DASTOOL or COMPARISON OF BINNING METHODS??

If DASTOOL, then outputs from it will be used in the next step. 
If we just compare methods, we could count the number of high-quality bins from each method and choose the better. 
Can be done with simple `awk| wc -l` one-liner.

THEN WE NEED MAG DEREPLICATION, MAPPING OF RECIPIENT READS TO MAGs, AND STRAIN ENGRAFTMENT ANALYSIS.

## Phage part

INSTRUCTIONS HERE

```bash
06_PHAGES/
```

## Strain engraftment with inStrain

The next step is to determine the strain engraftment of the donor MAGs to few selected recipients.  
Before we can run the workflow, we need to fetch the recipient data and map it against the donor MAGs.  

First make a whole new folder for all this. And put all output files and files we make in this folder.

```bash
mkdir -p 07_RECIPIENTS
cd 07_RECIPIENTS
```

### Fetch recipient data

We will download all the data from three recipients. Two from the FMT and one from the placebo group.  
The subject identifiers for these are:

```bash
TF13
TF29
TF45
```

Find all accession numbers for these subjects from SRA and write them down.  
You should have 12 read accessions to download.  

Put the recipient data inside a `Data` folder in the `07_RECIPIENTS` folder. Make sure you store them as compressed files (`.gz`). And it might be a good idea to make an array batch job for this task.  

```bash
mkdir Data
```

Make an array job that downloads the read files for each accession with `fasterq-dump` and also compresses them with `pigz`. Use the metaphlan script as an example and the same commands we used to download the donor metagenomes.  

```bash
sbatch src/YOUR_ARRAY_SCRIPT.sh
```

And of there were errors or problems with the download, you can copy the files from the course data folder (`/scratch/project_2022405/GutBugsData/Rawdata`).

### Mapping of recipient reads to donor MAGs

```bash
# we can use the same bowtie index
# and a modified version of the mapping.sh script
```

We might need this??  
__samples.txt:__

```bash
sample  r1  r2
TF13_12wk_FMT   Data/SRR11941425_1.fastq.gz    Data/SRR11941425_2.fastq.gz
TF13_BL_FMT     Data/SRR11941661_1.fastq.gz    Data/SRR11941661_2.fastq.gz
TF13_6wk_FMT    Data/SRR11941662_1.fastq.gz    Data/SRR11941662_2.fastq.gz
TF13_26wk_FMT   Data/SRR11941663_1.fastq.gz    Data/SRR11941663_2.fastq.gz
TF29_BL_FMT     Data/SRR11941593_1.fastq.gz    Data/SRR11941593_2.fastq.gz
TF29_6wk_FMT    Data/SRR11941594_1.fastq.gz    Data/SRR11941594_2.fastq.gz
TF29_26wk_FMT   Data/SRR11941595_1.fastq.gz    Data/SRR11941595_2.fastq.gz
TF29_12wk_FMT   Data/SRR11941596_1.fastq.gz    Data/SRR11941596_2.fastq.gz
TF45_BL_placebo Data/SRR11941426_1.fastq.gz    Data/SRR11941426_2.fastq.gz
TF45_6wk_placebo        Data/SRR11941485_1.fastq.gz    Data/SRR11941485_2.fastq.gz
TF45_26wk_placebo       Data/SRR11941486_1.fastq.gz    Data/SRR11941486_2.fastq.gz
TF45_12wk_placebo       Data/SRR11941487_1.fastq.gz    Data/SRR11941487_2.fastq.gz
```



Then we can import the resulting file `semibin_collection.txt` to anvi'o.  

```bash
module load anvio/8

anvi-import-collection \
    -C Semibin2 \
    -c 03_ANVIO/CONTIGS.db \
    -p 04_MAPPING/MERGED/PROFILE.db \
    --contigs-mode \
    08_AUTOMATED_BINNING/semibin_collection.txt
```

Then open the interactive interface and from "Bins" tab, click "Load bin collection" and select the correct collection. This will take some time, so be ___p-a-t-i-e-n-t___.  

You can also try to make similar figures from few bins as in the blog post [Visualizing the fate of contigs across metagenomic binning algorithms](https://merenlab.org/2020/01/02/visualizing-metagenomic-bins/). You need to export your own collection and semibin collection from anvi'o and then copy the script that was used to make those figures. And you also need to modify the script to read in the correct input files having the binning results. Read the blog post carefully and inspect the script and its comments.
