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
12. [Automatic binning](#automatic-binning)

## Introduction

On this course we will use a subset from a publicly available human gut metagenomic sequencing data from fecal microbiota transplantation (FMT) study: [Wilson, B.C., Vatanen, T., Jayasinghe, T.N. _et al._ Strain engraftment competition and functional augmentation in a multi-donor fecal microbiota transplantation trial for obesity. Microbiome 9, 107 (2021)](https://doi.org/10.1186/s40168-021-01060-7).  
There is also a 4-year follow-up study, if you're interested: [Long-term health outcomes in adolescents with obesity treated with faecal microbiota transplantation: 4-year follow-up](https://www.nature.com/articles/s41467-025-62752-4).  

Using read-based and assembly-based approaches, we will study the human gut microbiome. With read-based approach we will determine the taxonomic composition of one donor microbiome (female super-donor DF16) and compare to those the taxonomic composition of the recipient microbiomes taken at baseline and at 6, 12 and 26 weeks after the FMT. With assembly based approach, we will use [anvi'o](anvio.org) to construct metagenome-assembled genomes (MAGs) from the donor samples. Then we will assees the strain engrafment of selected strains in a selected set of recipients.  

## Setup

### SSH key

You need to set up SSH key before you can connect to CSC from your local machine using SSH.  
If you're not sure whether you have done it, you can check it from [my.csc.fi](my.csc.fi) under your profile and SSH public keys.  

In case you haven't set up SSH key, follow the instructions at [CSC webpages](https://docs.csc.fi/computing/connecting/ssh-keys/).  

__Save the SSH key to the default location.__  

### Connecting to Puhti with Visual Studio Code

* Launch Visual Studio Code
* Down left corner you will have a (green) button with "><" (hoover over it and it says "Open a Remote Window"), click it  
* Choose "Connect Current Window to Host..."
* Type in the __YOUR_USERNAME@puhti.csc.fi__ and hit "Enter" (change "user" for your own CSC username)  
* Type your password and hit "Enter"
* In the following dialogue, type __yes__ and hit "Enter"

### Notes for Windows user with VS Code

_If you have problems in connecting, follow the following instructions:_

* _When conneting to remote host open the config file that VS Cocde suggests_
* _Add a newline after your username with the following:_ `MACs hmac-sha2-512`

When the down left corner says `SSH:puhti.csc.fi`, you're connected.

* From the menu select `Terminal > New Terminal` and you should see a new panel. This is the __command line__.

* When you need to logout just type __exit__ in the terminal/command line and hit "Enter"  
(or you can click the down left corner and choose "Close Remote Connection")

### Setting up the course folders

The main course directory is located at `/scratch/project_2016640`.  
There you will set up your own directory where you will perform all the tasks for this course.  
Some of the tools needed on this course can be found from the course project applications folder `/projappl/project_2016640`.  

First list all projects you're affiliated with in CSC.

```bash
csc-workspaces
```

You should see the course project `MMB-901`.
So let's create a folder for you inside the scratch folder, you can find the path in the output from the previous command.

```bash
cd PATH/TO/COURSE/SCRATCH
mkdir $USER
```

Check with `ls`; which folder did `mkdir $USER` create?

This directory (`/scratch/project_2016640/your-user-name`) is your own working directory.  
Every time you log into Puhti, you should use `cd` to navigate to this directory.

Go to your own folder and clone the course repository there.  

```bash
git clone https://github.com/karkman/MMB-901_Metagenomics.git
```

Check what was downloaded and go to that folder. Then again check what is inside.  
__All the scripts are to be run in this folder__ (unless instructed otherwise).  

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

__Remember to always first open a screen session and only after that run `sinteractive`.__  

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

Finallly download only the accession numbers and make a file (`DF16_accessions.txt`) with only the accession, one line per accesion, to the data folder in Puhti (`01_DATA/`).  
We will use that file in the next step to download the data to Puhti.  

We will use a tool called [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) from NCBI to download the sequencing read files. It is installed on Puhti by default under the biokit module, so we need to load that first.  
But before downloading, have a look at the documentation to understand the different options we are using.  

Apply for resources. This doesn't take that long, so no real need for screen session.  

```bash
sinteractive -A project_2016640 -m 10G -c 8 -t 03:00:00
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
This is rather slow, so consider running it in screen. Before opening a screen, check how screen works from above.  

```bash
pigz -p $SLURM_CPUS_PER_TASK 01_DATA/*.fastq
```

After this is done, check that each sequence file in the folder `01_DATA` has the ending `.gz`.  

## Quality control

Before running any real analyses, we should do quality control (QC) for the sequencing data.  

We use two widely used programs that are pre-installed in Puhti:

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for running the QC fpor each sequence file.
* [MultiQC](https://multiqc.info/) to combine the individual reports from FastQC.  

```bash
cd /scratch/project_2016640/$USER/MMB-901_Metagenomics
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
Puhti uses a scheduling system called SLURM.
Make sure the files names, folders, options and project number are correct.  
Save the script as `src/spades.sh`.  

You can also use [Puhti](https://docs.csc.fi/computing/running/creating-job-scripts-puhti/) and [spades](https://github.com/ablab/spades#sec3.2) manuals, to learn more.  

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

To make things run a bit faster, we will run metaphlan as an [array job](https://docs.csc.fi/computing/running/array-jobs/). In a nutshell, all jobs will be run in parallel as individual jobs. This is a handy way to do the same thing for several files that are independent.  
CSC uses SLURM job scheduling system for array jobs.  

Use AI to prepare a batch job script for running metaphlan4 as an array job based on this information. Ask for a detailed description of the file content and the different options, so you understand what is going on.  
The You will need to request 8 CPUs, 50G of memory and 100G of local scratch space. The maximum time for the job is 1 hour.  
The version of metaphlan we'll use is 4.2.4 and the database is located at `/scratch/project_2016640/DBs/metaphlan`.  

Make sure the files names, folders, options and project number are correct and save the script as `src/metaphlan.sh`.  
You can also use [metaphlan4 documentation](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4.2#basic-usage) to learn more.  

And when you have a basic understanding what we are about to do, submit the job(s).  

```bash
sbatch src/metaphlan.sh
```

And again you can monitor the jobs with the same way as spades batch job.  

Taxonomic profiling doesn't take that long, approx 30 min per sample, but as they run in parallel, it will not take n x 30 min, but a lot less depending on the queue.  
While you wait, you can log in to Puhti web interface at [www.puhti.csc.fi](http://www.puhti.csc.fi) and open a Rstudio session (use the default options).  
We will analyse the results in R using few packages for microbiome data analysis.  

When all the metaphlan jobs are finished, we can go on with the data analysis of the microbial community.  
But before we can read the data into R, we need to combine the individual metaphlan outputs and extract the species level annotations from there.  

```bash
module load metaphlan/4.2.4
merge_metaphlan_tables.py 05_TAXONOMY/SRR*.txt > 05_TAXONOMY/metaphlan.txt
```

Then follow the R instruction in the file `src/taxonomic_profiling.r` and run the analysis in Puhti using the browser interface of Rstudio.  

After we have analysed the taxonomic profiles of the donor, we combine the rest of the samples to our merged metaphlan table and run the analysis again.  
First copy the taxonomic profiles of additional ~192 samples to the metaphlan output folder and re-run the merge command above.  

```bash
cp /scratch/project_2016640/GutBugsData/metaphlan/*.txt 05_TAXONOMY/

merge_metaphlan_tables.py 05_TAXONOMY/SRR*.txt > 05_TAXONOMY/metaphlan.txt
```

Then re-run the R part.  

## Assembly QC

Although it is not straightforward to assess the quality of a metagenomic assembly, we can still run a QC ananlysis and at least see what we got.  
We will use the metagenomic version of [QUAST](https://quast.sourceforge.net/docs/manual.html) for the job.  

QUAST can be found from Puhti, so just need to load the quast module.  
But first allocate some resources.  

```bash
sinteractive --account project_2016640
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

We will use [anvi'o](http://www.anvio.org) for genome-resolved metagenomics. Anvi'o is a multi-omics analysis and visualization software. The website has multiple tutorials and comprehensive manual on how to use it with different data sets and for various analyses. We will cover only a small part of it.  
There's also a [´omics vocabulary](https://anvio.org/vocabulary/), that can help you to understand some key terms and concepts in microbial 'omics.  

Allocate computing resources. We'll need 40G of memory, 6 CPUs and 100(G) of temporary storage for at least 4 hours.  

```bash
# sinteractive
```

### Contigs database

The first task is to create the contigs database from our assembled contigs. To make things run a bit smoother, we'll first remove all contigs < 2500nt.  
During the contigs database creation, anvi'o does gene calling with prodigal, calculates the tetranucleotide frequencies for each contigs and splits longer contigs into ~20 000 nt chunks called splits.  

```bash
mkdir -p 03_ANVIO
module load anvio/8

anvi-script-reformat-fasta 02_ASSEMBLY/contigs.fasta --min-len 2500 -o 03_ANVIO/contigs2500.fasta
anvi-gen-contigs-database -f 03_ANVIO/contigs2500.fasta -T 6 -o 03_ANVIO/CONTIGS.db
```

### Annotation of contigs database

After the contigs database has been created, we'll add only few annotations to the database that aid in the binning process. There are many other commands to add other annotations to the database.  
`anvi-run-hmms` annotates single-copy core genes and ribosomal RNA's that are used to estimate the completeness and contamination of a bin.  
`anvi-run-scg-taxonomy` annotates single-copy core genes with taxonomic information.  

```bash
anvi-run-hmms -c 03_ANVIO/CONTIGS.db -T $SLURM_CPUS_PER_TASK
anvi-run-scg-taxonomy -c 03_ANVIO/CONTIGS.db -T $SLURM_CPUS_PER_TASK
```

### Mapping and profiling

To obtain the differential coverage information for each contig in our database, we will map the original reads from each sample individually to the formatted fasta file (contigs < 2500 nt removed) and make a single profile database from each formatted mappign output file.  
We used [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for the mapping and the first step is to create an index of the fasta file.  

```bash
mkdir -p 04_MAPPING
bowtie2-build 03_ANVIO/contigs2500.fasta 04_MAPPING/contigs --threads $SLURM_CPUS_PER_TASK
```

As we will be running the same job for all nine samples, the quickest way is to create an array job for the task.  
Copy the `metaphlan.sh` script to a new file called `mapping.sh` in the same `src` folder and modify it accordingly:

* Change the name of the job and the names of the log files
* Change the time to 1 hour
* Change the memory to 12G
* Instead of metaphlan, load `anvio/8` module

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
    -S $LOCAL_SCRATCH/${SAMPLE_ACC}.sam

samtools view -bS $LOCAL_SCRATCH/${SAMPLE_ACC}.sam > $LOCAL_SCRATCH/${SAMPLE_ACC}.bam

anvi-init-bam $LOCAL_SCRATCH/${SAMPLE_ACC}.bam -T $SLURM_CPUS_PER_TASK -o 04_MAPPING/${SAMPLE_ACC}.bam

anvi-profile \
    -i 04_MAPPING/${SAMPLE_ACC}.bam \
    -c 03_ANVIO/CONTIGS.db \
    --output-dir 04_MAPPING/${SAMPLE_ACC} \
    -S ${SAMPLE_ACC} \
    -T $SLURM_CPUS_PER_TASK \
    --min-contig-length 5000
```

After the script is ready, you can submit the job(s).  

```bash
sbatch src/mapping.sh
```

### Merging the profiles

When all the mapping jobs have been finished, we can merge all the single profiles into a merged profile database.  
For this task you will need 40G of memory and about 30 min.  

```bash
anvi-merge 04_MAPPING/SRR*/PROFILE.db -o 04_MAPPING/MERGED -c 03_ANVIO/CONTIGS.db --enforce-hierarchical-clustering
```

### Interactive use and binning

Now we should have all that we need to visualize our metagenome in anvi'o and start binning. This part is easiest to do with VS Code, but other options are also possible. Except not in Puhti web interface (as far as I know).  
To be able to tunnel the interactive interface from Puhti to your own computer, everyone will need a separate port.  
We will allocate the ports in class.  

When you know your port number (`XXXX`), assign it to the environemntal variable `$ANVIO_PORT`.  

```bash
module load anvio/8
ANVIO_PORT=XXXX
```

For the interactive work you will need around 10G of memory and only 1 CPU.  

Then we need to set up the ports in VS Code. We'll go this thru together. If you're using other software for tunneling, figure out how to set the port forwarding there.  

Forward the port in VS Code under PORTS tab as follows: `NODEID.bullx:PORT`  

After the port has been set, you can run the `anvi-interactive` to launch the interactive interface.  

```bash
anvi-interactive -c 03_ANVIO/CONTIGS.db -p 04_MAPPING/MERGED/PROFILE.db -P $ANVIO_PORT
```

Other usefull commands we will need during binnning.  

```bash
anvi-refine -c 03_ANVIO/CONTIGS.db -p 04_MAPPING/MERGED/PROFILE.db -P $ANVIO_PORT -C PreCluster -b Bin_1

anvi-rename-bins \
    -c 03_ANVIO/CONTIGS.db \
    -p 04_MAPPING/MERGED/PROFILE.db \
    --collection-to-read PreCluster \
    --collection-to-write Bins \
    --prefix DF16 \
    --report-file 03_ANVIO/Bins_report.txt

anvi-summarize -c 03_ANVIO/CONTIGS.db -p 04_MAPPING/MERGED/PROFILE.db -C Bins -o 03_ANVIO/SUMMARY_Bins
```

## MAG QC and taxonomy

When the binning is done, it's time to assess the quality of them and give a taxonomic annotation to each bin.  
We will use [CheckM2](https://github.com/chklovski/CheckM2) for assessing the quality of our MAGs and [GTDB-Tk](https://github.com/ecogenomics/gtdbtk) for the taxonomic annotation.  

But since there are probably a lot of low quality bins (low completeness), let's run `anvi-rename-bins` and `anvi-summarize` once more and ask anvi'o to mark all bins it thinks are good quality. This way we can reduce the amount of bins we need to analyse and everything will be a bit faster. We will use a bit more strict cut-off for the completion to concentrate only on the high-quality MAGs.  

Allocate resources for the job. You will need 10G of memory and 30 min.  

```bash
sinteractive -A project_2016640 -t 00:30:00 -m 10G
```

```bash
anvi-rename-bins \
    -c 03_ANVIO/CONTIGS.db \
    -p 04_MAPPING/MERGED/PROFILE.db \
    --collection-to-read \ # The latest collection you have goes here, before the "\"
    --collection-to-write MAGs \
    --call-MAGs \
    --min-completion-for-MAG 90 \
    --max-redundancy-for-MAG 10 \
    --prefix DF16 \
    --report-file 03_ANVIO/MAG_report.txt

anvi-summarize \
    -c 03_ANVIO/CONTIGS.db \
    -p 04_MAPPING/MERGED/PROFILE.db \
    -o 03_ANVIO/SUMMARY_MAGs \
    -C MAGs
```

Now you can release the resources.  

```bash
exit
```

### MAG folder

Now each bin that had completeness and redundancy (according to anvi'o) over 70 % and under 10 %, respectively, will have "MAG" in its name. For example `DF16_MAG_00001`. So we can only pick those for further analyses. We will make a folder with softlinks to all of our MAGs and then analyze only these with CheckM2 and GTDB-Tk.  

```bash
mkdir -p 06_GENOMES
cd 06_GENOMES
ln -s ../03_ANVIO/SUMMARY_MAGs/bin_by_bin/*MAG*/*MAG*-contigs.fa ./
cd ..
```

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

### MAG QC with CheckM2

```bash
export CHECKM2DB="/scratch/project_2016640/DBs/CheckM2_database/uniref100.KO.1.dmnd"

/projappl/project_2016640/tax_tools/bin/checkm2 predict \
    --input 06_GENOMES \
    --output-directory 06_GENOMES/checkm2 \
    --extension fa \
    --tmpdir $LOCAL_SCRATCH \
    --threads $SLURM_CPUS_PER_TASK
```

### MAG taxonomy with GTDB-Tk

```bash
export GTDBTK_DATA_PATH="/scratch/project_2016640/DBs/release226/"

/projappl/project_2016640/tax_tools/bin/gtdbtk classify_wf \
    --genome_dir 06_GENOMES \
    --out_dir 06_GENOMES/gtdbtk \
    --extension fa \
    --skip_ani_screen \
    --cpus $SLURM_CPUS_PER_TASK \
    --pplacer_cpus 1 \
    --scratch_dir $LOCAL_SCRATCH \
    --tmpdir $LOCAL_SCRATCH
```

## MAG annotation

After quality control and taxonomic annotation of all MAGs, we will choose two for strain engraftment analysis. The Fig. 2 in the original publication can help you choose MAGs that could be of interest. We want to determine whether the donor strains have colonised the recipients. We can talk together which ones you could choose, but keep in mind that we might not have good quality genomes from all of those, so you need to use the CheckM2 and GTDB-Tk results to verify what you have.  

When you have picked two, annotate them both with Bakta using the following command. Make sure the path to the genome is right (`GENOME_BIN`) and use the genus level annotation of the MAG as `GENOME_NAME`.  

And of course allocate some resources: 4 CPUs, 20Gb of memory and 1 hour. It takes around 15-20 min per genome.  

```bash
sinteractive -A project_2016640 ...
```

Set the Bakta database path.  

```bash
export BAKTA_DB=/scratch/project_2016640/DBs/bakta/db
```

```bash
/projappl/project_2016640/bakta/bin/bakta \
    06_GENOMES/GENOME_BIN.fa  \
    --skip-pseudo \
    --skip-sorf \
    --prefix GENOME_NAME \
    --locus GENOME_NAME \
    --threads $SLURM_CPUS_PER_TASK \
    --output 06_GENOMES/GENOME_NAME 
```

## Strain engraftment

The next step is to determine the strain engraftment of the selected MAGs from the donor to few selected recipients.  
We will use anvio workflows, which is a snakemake wrapper for anvi'o. You can learn more about anvi'o workflows here: [https://merenlab.org/2018/07/09/anvio-snakemake-workflows](https://merenlab.org/2018/07/09/anvio-snakemake-workflows) and about Snakemake, the workflow manager used with anvi'o workflows, from here: [https://snakemake.readthedocs.io/en/stable/](https://snakemake.readthedocs.io/en/stable/).  

Before we can run the workflow, we need to fetch the recipient data, pre-process the annotated genome files (Genbank files from Bakta) for anvi'o and prepare few files for the workflow.  

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

Put the recipient data inside a `Data` folder in the `07_RECIPINTS` folder. Make sure you store them as compressed files (`.gz`). And it might be a good idea to make an array batch job for this task.  

```bash
mkdir Data
```

Make an array job that downloads the read files for each accession with `fasterq-dump` and also compresses them with `pigz`. Use the metaphlan script as an example and the same commands we used to download the donor metagenomes.  

```bash
sbatch src/YOUR_ARRAY_SCRIPT.sh
```

### Process Genbank files

Process both selected MAGs. The input is the genbank file (`.gbff`) in the bakta output folder.  
Write the output files to a new folder called `Genomes` in our `07_RECIPIENTS` folder. Add the genome name as the prefix (option `-O`).  

```bash
mkdir Genomes

module load anvio/8

anvi-script-process-genbank \
    -i PATH/TO/GENOME_NAME.gbff \
    -O PATH/TO/GENOME_NAME \
    --annotation-source bakta \
    --annotation-version 1.5.1
```

### Anvi'o workflow

We need to create three files for the workflow. One with the location of the genome files called `fasta.txt`. The second files has the paths to the read files called `samples.txt`. And the last file is a configuration file that has the instructions for what do run in the workflow. This file will be called `config.json`. Put all files in the `07_RECIPINTS` folder.  

Examples of each file below. The first two have to be tab-separated. The third one is a JSON file and has a very different syntax. But this file you do not need to change.  

__fasta.txt:__

```bash
name    path    external_gene_calls gene_functional_annotation
GENOME_NAME Genomes/GENOME_NAME-contigs.fa   Genomes/GENOME_NAME-external-gene-calls.txt  Genomes/GENOME_NAME-external-functions.txt
GENOME_NAME Genomes/GENOME_NAME-contigs.fa   Genomes/GENOME_NAME-external-gene-calls.txt  Genomes/GENOME_NAME-external-functions.txt
```

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

__config.json:__  

```bash
{
    "fasta_txt": "fasta.txt",
    "anvi_gen_contigs_database": {
        "--project-name": "{group}",
        "--description": "",
        "--skip-gene-calling": "",
        "--ignore-internal-stop-codons": "",
        "--skip-mindful-splitting": "",
        "--contigs-fasta": "",
        "--split-length": "",
        "--kmer-size": "",
        "--skip-predict-frame": "",
        "--prodigal-translation-table": "",
        "threads": 6
    },
    "anvi_run_hmms": {
        "run": true,
        "threads": 6,
        "--also-scan-trnas": true,
        "--installed-hmm-profile": "",
        "--hmm-profile-dir": ""
    },
    "samples_txt": "samples.txt",
    "bowtie": {
        "additional_params": "--no-unal",
        "threads": 4
    },
    "samtools_view": {
        "additional_params": "-F 4",
        "threads": 2
    },
    "anvi_profile": {
        "threads": 4,
        "--sample-name": "{sample}",
        "--overwrite-output-destinations": true,
        "--report-variability-full": "",
        "--skip-SNV-profiling": "",
        "--profile-SCVs": "",
        "--description": "",
        "--skip-hierarchical-clustering": "",
        "--distance": "",
        "--linkage": "",
        "--min-contig-length": "",
        "--min-mean-coverage": "",
        "--min-coverage-for-variability": "",
        "--cluster-contigs": "",
        "--contigs-of-interest": "",
        "--queue-size": "",
        "--write-buffer-size-per-thread": "",
        "--max-contig-length": ""
    },
    "anvi_merge": {
        "--sample-name": "{group}",
        "--overwrite-output-destinations": true,
        "--description": "",
        "--skip-hierarchical-clustering": "",
        "--enforce-hierarchical-clustering": true,
        "--distance": "",
        "--linkage": "",
        "threads": "6"
    },
    "import_percent_of_reads_mapped": {
        "run": true,
        "threads": ""
    },
    "annotate_contigs_database": {
        "threads": ""
    },
    "bowtie_build": {
        "additional_params": "",
        "threads": 6
    },
    "anvi_init_bam": {
        "threads": 4
    },
    "references_mode": "true",
    "all_against_all": "true",
    "kraken_txt": "",
    "collections_txt": "",
    "output_dirs": {
        "FASTA_DIR": "02_FASTA",
        "CONTIGS_DIR": "03_CONTIGS",
        "QC_DIR": "01_QC",
        "MAPPING_DIR": "04_MAPPING",
        "PROFILE_DIR": "05_ANVIO_PROFILE",
        "MERGE_DIR": "06_MERGED",
        "TAXONOMY_DIR": "07_TAXONOMY",
        "SUMMARY_DIR": "08_SUMMARY",
        "SPLIT_PROFILES_DIR": "09_SPLIT_PROFILES",
        "LOGS_DIR": "00_LOGS"
    },
    "max_threads": 12,
    "config_version": "3",
    "workflow_name": "metagenomics"
}
```

When all the files have been created, check that everything is formatted correctly by doing a dry run.  

```bash
module load anvio/8
anvi-run-workflow --workflow metagenomics --config-file config.json --dry-run
```

Then when all the files are ready and correctly formatted, you can either create a batch file to run the workflow or do it interactively.  
No matter which way you choose, you will need 12 CPUs, 50G of memory and 2 hours.  

And the commands to run the workflow are below. Make sure you run it inside the `07_RECIPIENTS` folder.  

```bash
module load anvio/8
anvi-run-workflow --workflow metagenomics --config-file config.json
```

## Automatic binning

There are several different tools for automated binning and they all perform very differently (as you probably read in the blog post). We will use [SemiBin2](https://github.com/BigDataBiology/SemiBin) for the automated binning.  

Semibin2 uses deep learning in metagenomic binning and has pre-trained models for several different environments, including human gut. However, with multiple samples the pre-trained models cannot be used and we would need to train the model with our data. The model training takes some time, although it might be more accurate. Running Semibin2 with pre-trained model is a lot faster and we will prefer speed over accuracy on this course. In case we have time, we can try the multi-sample mode for our data.  

We will use the same files we created and used in the manual binning with anvi'o. The job can be run interactively or as a batch job. You will need at least 6 CPUs, 50G of memory, 100G of local storage and many hours (to be sure this time). My test run took 35 min.  
Allocate the resources or write a batch job script and use the following command to run Semibin2.  

Also make sure to run this from our course main folder (MMB-901_Metagenomics) so that all the paths are right.  

```bash
/projappl/project_2016640/Semibin2/bin/SemiBin2 single_easy_bin \
    --input-fasta 03_ANVIO/contigs2500.fasta \
    --input-bam 04_MAPPING/SRR11941565.bam \
    --environment human_gut \
    --tmpdir $LOCAL_SCRATCH \
    --threads $SLURM_CPUS_PER_TASK \
    --output 08_AUTOMATED_BINNING
```

The outputs from Semibin2 will be written to `08_AUTOMATED_BINNING` and the folder will contain information about the bins, some log files and one folder with all the genome bins as compressed fasta files.  

If you want to try multi-sample binning, re-run the previous command with all the BAM-files and remove the option `--environment human_gut`. Also remember to change the name of the output folder. This will take significantly longer, so make a batch job with more time and more cores (10-16 CPUs).  

After the binning is ready, you can run CheckM2 and GTDB-Tk for the resulting bins and compare the results to your own binning. Both tools accept the genome fasta files in compressed format (`.fa.gz`). Just specify the extension correctly.  

We can also import the binning results to anvi'o and visually inspect whether we agree with Semibin2 or not.  
Semibin2 produces a file (`contig_bins.tsv`) that has all the needed information, but the format is not exactly as anvi'o would wish.
So we need to do some command line magic. Or in other words use `awk` to remove the first line and add a prefix to each bin name (as anvi'o does not accept only numbers as bin names).

```bash
awk 'NR!=1 {print $1"\tSemibin2_"$2}' 08_AUTOMATED_BINNING/contig_bins.tsv > 08_AUTOMATED_BINNING/semibin_collection.txt
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
