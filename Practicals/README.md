# MMB-901 Microbial metagenomics – Practicals

__Table of Contents:__

1. [Introduction](#introduction)
2. [Connecting to Roihu](#connecting-to-roihu)
3. [Setup](#setup)
4. [Data](#data)
5. [Quality control](#quality-control)
6. [Metagenome assembly](#metagenome-assembly)
7. [Read-based taxonomy](#read-based-taxonomy)
8. [Assembly QC](#assembly-qc)
9. [Genome-resolved metagenomics](#genome-resolved-metagenomics)
10. [MAG taxonomy](#mag-taxonomy)
11. [Viral metagenomics](#viral-metagenomics)
12. [Strain engraftment](#strain-engraftment-with-instrain)

## Introduction

On this course we will use a subset from a publicly available human gut metagenomic sequencing data from fecal microbiota transplantation (FMT) study: [Wilson, B.C., Vatanen, T., Jayasinghe, T.N. _et al._ Strain engraftment competition and functional augmentation in a multi-donor fecal microbiota transplantation trial for obesity. Microbiome 9, 107 (2021)](https://doi.org/10.1186/s40168-021-01060-7).  
There is also a 4-year follow-up study, if you're interested: [Long-term health outcomes in adolescents with obesity treated with faecal microbiota transplantation: 4-year follow-up](https://www.nature.com/articles/s41467-025-62752-4).  

Using read-based and assembly-based approaches, we will study the human gut microbiome. With read-based approach we will determine the taxonomic composition of one donor microbiome (female super-donor DF16) and compare to those the taxonomic composition of the recipient microbiomes taken at baseline and at 6, 12 and 26 weeks after the FMT. With assembly based approach, we will construct metagenome-assembled genomes (MAGs) from the donor samples. Then we will assess the strain engrafment of selected strains in a selected set of recipients.  

## Connecting to Roihu

If you not experienced with UNIX and command line, the easiest way is to use the Roihu web interface, [roihu.csc.fi](https://roihu.csc.fi/). You can log in with Haka or CSC crecentials.  
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

[Tmux](https://docs.csc.fi/support/tutorials/tmux/index.html) is a handy way to run things in the background without losing them when you logout or have connection problems. However, you have to be careful when using tmux, you can easily get lost.  
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

Finally download only the accession numbers and make a file (`DF16_accessions.txt`) with only the accession, one line per accession, to the data folder in Roihu (`01_DATA/`).  
We will use that file in the next step to download the data to Roihu.  

We will use a tool called [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) from NCBI to download the sequencing read files. It is installed on Roihu by default under the bio-apps module, so we need to load that first.  
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

module load bio-apps/v202603
module load sra-tools/3.3.0

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

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for running the QC for each sequence file.
* [MultiQC](https://multiqc.info/) to combine the individual reports from FastQC.  

```bash
cd /scratch/project_2022405/$USER/MMB-901_Metagenomics
mkdir 01_DATA/FASTQC

module load bio-apps/v202603
module load fastqc/0.12.1

fastqc --threads $SLURM_CPUS_PER_TASK --outdir 01_DATA/FASTQC 01_DATA/*.fastq.gz

module load py-multiqc/1.35

multiqc --interactive --outdir 01_DATA/FASTQC 01_DATA/FASTQC/
```

After running the QC steps, download the MultiQC report (`.html`) file to your local computer and open it with browser.  
We'll go through the report together.  

## Metagenome assembly

As the donor samples are from the same individual, we can do a co-assembly with all the data. For the co-assembly, we'll combine all R1 reads to one file and R2 files to another.  

```bash
cat 01_DATA/SRR*_1.fastq.gz > 01_DATA/DF16_1.fastq.gz
cat 01_DATA/SRR*_2.fastq.gz > 01_DATA/DF16_2.fastq.gz
```

Then the actual assembly will be done with [spades](https://github.com/ablab/spades) using the `--meta` option meant for metagenomic data and `--only-assembler` to skip the read correction part.  Add  `--memory 400` to ensure spades limits memory to what is requested in the SBATCH job.

As this will take longer, we'll run it as batch job. The scheduling system in Roihu is called SLURM.  

Using [Roihu](https://docs.csc.fi/computing/running/creating-job-scripts-roihu/) and [spades](https://github.com/ablab/spades#sec3.2) manuals, prepare a batch job script for running metagenomic assembly based on this information.  
If using AI, ask for a detailed description of the file content and the different options, so you understand what is going on.  

You will need to request 12 CPUs, 400G of memory. The maximum time for the job is 16 hours.

Use --tmp-dir $TMPDIR option to tell which location to use for temporary memory. Add `module load bio-apps` and `module load spades/4.3.0` lines before you run the spades.

Use the paired-end FASTQ files DF16_1.fastq.gz and DF16_2.fastq.gz from 01_DATA/FASTQ/.

Make sure the files names, folders, options and project number are correct.  
Save the script as `src/spades.sh`.  

You can also use  and  manuals, to learn more.  

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
The You will need to request 8 CPUs and 50G of memory. The maximum time for the job is 1 hour.

The version of metaphlan we'll use is 4.2.4 and it is installed on Roihu under the module `py-metaphlan/4.2.4` (needs also `bio-apps` module).

Hint for the assay job: create a text file that lists the accession numbers as shown below, and include it in the Slurm job. Then the metaphlan script can call the $SAMPLE_ACC variable.

```
SAMPLE_ACC=$(sed -n ${SLURM_ARRAY_TASK_ID}p 01_DATA/DF16_accessions.txt)
```
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
module load bio-apps/v202603
module load py-metaphlan/4.2.4
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

Although it is not straightforward to assess the quality of a metagenomic assembly, we can still run a QC analysis and at least see what we got.  
We will use the metagenomic version of [QUAST](https://quast.sourceforge.net/docs/manual.html) for the job.  

QUAST can be found from Roihu, so just need to load the quast module.  
But first allocate some resources.  

```bash
sinteractive --account project_2022405
```

While you wait for the resources, have a look at the quast manual and read about the options were using.  
Then when you're connected to a computing node and have read about the options from the manual, run quast.

```bash
module load quast/5.3.0 
metaquast.py 02_ASSEMBLY/contigs.fasta --max-ref-num 0 --threads $SLURM_CPUS_PER_TASK -o 02_ASSEMBLY/QUAST --fast
```

When the job is finished, free the resources with:  

```bash
exit
```

And have a look at the output folder and inspect the results. QUAST manual will help you with the different output files.  


We will filter the assemblies using seqkit. This only takes a few seconds to run and little memory.

```
module load bio-apps/v202603
module load seqkit/2.13.0
seqkit seq -m 1000 contigs.fasta > filtered_contigs.fasta
```
## Genome-resolved metagenomics

Instead of manual binning, we will use automatic binning software to reconstruct some of the original genomes (metagenome-assembled genomes, MAGs) in our donor samples.  
There are several different tools for automated binning and they all perform very differently (as you probably read in the blog post).  
The programs we will use are [Metabat](LINK) and [Semibin2](https://github.com/BigDataBiology/SemiBin). Both use the differential coverage information and the sequence composition of the contigs as the basis to cluster them into bins.  
There's also a [´omics vocabulary](https://anvio.org/vocabulary/), that can help you to understand some key terms and concepts in microbial 'omics.  

### Mapping

To obtain the differential coverage information for each contig in our database, we will map the original reads from each sample individually to the assembly.  
We use [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for the mapping and the first step is to create an index of the fasta file. Indexing can be done interactively or as a batch job, using 20 min of time and 20 GB of memory.

```bash
mkdir -p 04_MAPPING
module load bio-apps/v202603
module load bowtie2/2.5.4

bowtie2-build 02_ASSEMBLY/filtered_contigs.fasta 04_MAPPING/contigs --threads $SLURM_CPUS_PER_TASK
```

As we will be running the same job for all nine samples, the quickest way is to create an array job for the task.  
Copy the `metaphlan.sh` script to a new file called `mapping.sh` in the same `src` folder and modify it accordingly:

* Change the name of the job and the names of the log files
* Change the time to 50 minutes
* Change the memory to 12G
* Instead of metaphlan, load bowtie2 and samtools modules

Remember that each of the array jobs reads the `DF16_accessions.txt` file and picks the right sample to analyse.  
In more detail; each array has its own ID from 1–9 stored in the environmental variable `SLURM_ARRAY_TASK_ID`. Each arrays reads the corresponding line from the file and stores it to the variable `SAMPLE_ACC`.  
Make sure the path to the accession file is correct in the following line: `SAMPLE_ACC=$(sed -n ${SLURM_ARRAY_TASK_ID}p 01_DATA/DF16_accessions.txt)`  

__Do not run__ the following, but these are the mapping and profiling commands that you need to change in the file:  

```bash
mkdir -p 04_MAPPING
module load bio-apps/v202603
module load bowtie2/2.5.4
module load samtools/1.24

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

### MetaBAT2

MetaBAT2 is an automated metagenomic binning tool that groups assembled contigs into metagenome-assembled genomes (MAGs). It uses sequence composition (tetranucleotide frequencies) and contig abundance (coverage) to identify groups of contigs that likely originate from the same microorganism. Unlike SemiBin2, MetaBAT2 does not use deep learning or require model training, which makes it relatively straightforward and computationally efficient.

We will use the assembly files and mapped reads created in the previous step. MetaBAT2 requires the assembled contigs in FASTA format and a contig-depth file generated from sorted BAM files. The depth file can be generated using jgi_summarize_bam_contig_depths, which is included with MetaBAT2.

Since we did a co-assembly, we will use all nine BAM files together in one MetaBAT2 run for the binning. The resulting genome bins can then be compared with the bins obtained using SemiBin2 (below).The job will use 8 CPUs and 32 GB of memory, with an initial time allocation of 4 hours.  

Create a folder for MetaBAT2 binning:  

```bash
05_BINNING/METABAT2/
```

To run MetaBAT2, we need three files:

1. Assembled (filtered) contigs (.fa / .fasta) to be grouped into MAGs
2. Sorted bam files (.bam)
3. Depth file (.txt) which will be created in the script

 Create a bash script in your src directory and run it using sbatch, or run interactively using 5G of memory, 30 min of time, and 8 CPUs. If you use the script below, make sure you change all paths and filenames accordingly.  

```bash
#!/bin/bash
#SBATCH --job-name=metabat2
#SBATCH --account=YOUR_PROJECT
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --output=metabat2_%j.out
#SBATCH --error=metabat2_%j.err

# Load MetaBAT2
module load bio-apps/v202603
module load metabat/2.15

# Input and output directories
CONTIGS="PATH_TO_FILTERED_CONTIG_FILE"
BAM_DIR="PATH_TO_BAM_DIR"
OUT_DIR="05_BINNING/METABAT2_OUT"

mkdir -p "$OUT_DIR"

# List all sorted BAM files using accession numbers
BAMS=()
for SAMPLE_ACC in $(cat 01_DATA/DF16_accessions.txt); do
    BAMS+=("$BAM_DIR/${SAMPLE_ACC}.sorted.bam")
done

# Generate contig depth file
jgi_summarize_bam_contig_depths \
    --outputDepth "$OUT_DIR/depth.txt" \
    "${BAMS[@]}"

# Run MetaBAT2
metabat2 \
    -i "$CONTIGS" \
    -a "$OUT_DIR/depth.txt" \
    -o "$OUT_DIR/bin" \
    -t "$SLURM_CPUS_PER_TASK"
```

Run the script:

```bash
sbatch src/metabat2.sh
```

### Semibin2

Semibin2 uses deep learning in metagenomic binning and has pre-trained models for several different environments, including human gut. However, with multiple samples the pre-trained models cannot be used and we would need to train the model with our data. The model training takes some time, although it might be more accurate. Running Semibin2 with pre-trained model is a lot faster and but will prefer accuracy over speed on this course and use the differential coverage information from the nine samples to train the model.  

We will use the same files as with Metabat. The job can be run interactively or as a batch job. You will need at least 8 CPUs, 60 GB of memory, and 3 h.
Allocate the resources or write a batch job script and use the following command to run Semibin2.  

Also make sure to run this from our course main folder (MMB-901_Metagenomics) so that all the paths are right.  

```bash
module load bio-apps/v202603
module load py-semibin/2.6.0

SemiBin2 single_easy_bin \
    --input-fasta 02_ASSEMBLY/filtered_contigs.fasta \
    --input-bam 04_MAPPING/*.bam \
    --tmpdir $TMPDIR \
    --threads $SLURM_CPUS_PER_TASK \
    --output 05_BINNING/SEMIBIN2/
```

The outputs from Semibin2 will be written to `05_BINNING/SEMIBIN2/` and the folder will contain information about the bins, some log files and one folder with all the genome bins as compressed fasta files.  

Semibin2 produces a file (`contig_bins.tsv`) that has all the needed information.
So we need to do some command line magic. Or in other words use `awk` to remove the first line and add a prefix to each bin name.

```bash
awk 'NR>1 {print $1"\tSemiBin2_"$2}' 05_BINNING/SemiBin2/contig_bins.tsv > 05_BINNING/SemiBin2/semibin_collection.txt
```

## DASTool

[DASTool](https://github.com/cmks/DAS_Tool) is an automated method that integrates the results of binning algorithms to calculate an optimized, non-redundant set of bins.

Familiarize yourself with the GitHub pages of DASTool and prepare a bash script to modify the input using the Fasta_to_Contigs2Bin.sh helper script. For the SemiBin2 contigs, you will need to gunzip them first. 

Use the code below to save the uncompressed files in $TMPDIR.

Remember to use the $TMPDIR/SemiBin2 as input for your Fasta_to_Contigs2Bin.sh helper script. Use the interactive node with 5GB memory, 8 CPUs, and 1 h time.
Direct the output to a new directory /scratch/project_2022405/$USER/MMB-901_Metagenomics/06_DASTOOL/DF16.

```bash
# Load DAS Tool
module load bio-apps/v202603
export PATH="/projappl/project_2022405/das_tool-1.1.7/bin:$PATH"

# Convert SemiBin2 bins (.fa.gz)
mkdir -p "$TMPDIR/SemiBin2"

for f in "$SEMIBIN"/*.fa.gz; do
    [[ -e "$f" ]] || continue
    gunzip -c "$f" > "$TMPDIR/SemiBin2/$(basename "${f%.gz}")"
done
```


Once the files are prepared, we will run DASTool with the options --write_bin_evals to write evaluations of input bin sets and --write_bins to export bins as fasta files. Use the same interactive session to run DASTool.

```bash
# Load DAS Tool
module load bio-apps/v202603
export PATH="/projappl/project_2022405/das_tool-1.1.7/bin:$PATH"

# Output directory
OUTDIR="/scratch/project_2022405/$USER/MMB-901_Metagenomics/06_DASTOOL/DF16"
CONTIGS="/scratch/project_2022405/$USER/MMB-901_Metagenomics/02_ASSEMBLY/DF16/contigs.fasta"

# Run DAS Tool
DAS_Tool \
    -i "$OUTDIR/SemiBin2_contigs2bin.tsv,$OUTDIR/MetaBAT2_contigs2bin.tsv" \
    -l SemiBin2,MetaBAT2 \
    -c "$CONTIGS" \
    -o "$OUTDIR/DF16_DASTool" \
    --write_bins \
    --write_bin_evals \
    --threads "$SLURM_CPUS_PER_TASK"
```

After DASTool has run open the  DF16_DASTool_DASTool_summary.tsv on your computer to view the results. Then create a list of all high-quality MAGs (HQ bins) with SCG_completeness 90 or more and SCG_redundancy 5 or less in a file named DF16_high_quality_bins.txt.

Check how many bins you had and how many of them were high quality before proceeding.

Then write a script to move the high-quality bins to a new directory, HQ_BINS, using DF16_high_quality_bins.txt to select the correct bins to move.

## dRep

Usually, we need to  dereplicate the bins for further analyses. This can be done using dRep. We will also prepare some files for or strain sharing analysis using dRep. Check the [dRep](https://github.com/MrOlm/drep) manual to see what is does. By chance, the DasTool genomes are already ANI 98% or more divergent, so we don't need to do further dereplication.

However, we will need to run [parse_stb.py](https://github.com/MrOlm/drep/blob/master/helper_scripts/parse_stb.py) to determine which scaffold belongs to which genome before we do strain analyses.  
Export the path and run the Python program on the fasta files of the HQ bins. We will need the representative_genomes.stb file later for InStrain analysis.

```bash
export PATH="/projappl/project_2022405/drep-4-0-3/bin:$PATH"
parse_stb.py --reverse -f HQ_BINS/*fa -o representitve_genomes.stb
```

## MAG taxonomy

When the binning and quality assessment are done, we need to give a taxonomic annotation to each bin.

We will use [GTDB-Tk](https://github.com/ecogenomics/gtdbtk) for the taxonomic annotation.  

The MAG taxonomic annotation can be run either interactively or as a batch job.

Below you can find the scripts. Resources for the batch/interactive jobs are as follows:

__GTDB-Tk:__

* 50G of memory
* 30 minutes of time
* 10 CPUs

You can use the spades script as a template.

### MAG taxonomy with GTDB-Tk

```bash
module load  bio-apps/v202603
module load py-gtdbtk/2.7.2

gtdbtk classify_wf \
    --genome_dir HQ_BINS \
    --out_dir HQ_BINS/gtdbtk \
    --extension fa \
    --skip_ani_screen \
    --cpus $SLURM_CPUS_PER_TASK \
    --pplacer_cpus 1 \
    --scratch_dir $TMPDIR \
    --tmpdir $TMPDIR
```

## Viral metagenomics

Prophages replicate with their host, so if a donor strain engrafts, its prophage genes arrive in the recipient gut with it. Some phage-encoded genes benefit the host bacterium, and in doing so the prophage itself. Anti-phage immunity genes are a clear case: carried within a prophage, they shield the host from other phages, and therefore the prophage from rivals.

__Research question:__ which donor prophages carry putative anti-phage system genes, and can we find them in the engrafting bacterial strains?

In this exercise you will predict prophages in donor MAGs, annotate their genes, and shortlist candidates. Keep in mind that annotation is incomplete, so absence means little, and that presence is a hypothesis about host fitness rather than evidence.

Make a directory for all virus analyses in your own directory:
```bash
mkdir 07_PHAGES
```

### Viral contig identification with geNomad

There are many different tools for predicting viral contigs from metagenomes. In this course, we will use [geNomad](https://github.com/apcamargo/genomad), which can identify viral and proviral sequences, and is a commonly used tool in viral metagenomics. 

__How does geNomad identify viral sequences?__
Check from the publication [Camargo, A.P., Roux, S., Schulz, F. et al. Identification of mobile genetic elements with geNomad. Nat Biotechnol 42, 1303–1312 (2024)](https://www.nature.com/articles/s41587-023-01953-y) 

Good documentation also [here](https://portal.nersc.gov/genomad/pipeline.html)

We will run geNomad on the high-quality bins to identify putative proviral sequences. Since geNomad generates a separate output directory for each input FASTA file, we will first concatenate all bin FASTA files into a single file to obtain a unified output. To preserve information about the original bin assignments, we will add a bin-specific prefix to each contig header before combining the files. (This can be run in the login node.)

```bash
for fasta in /path-to-high-quality-bins/*.fa; do
  bin=$(basename "$fasta" .fa)
  sed "s/^>/>${bin}__/" "$fasta"
done > 07_PHAGES/all_hq_bins.fna
```

We now have a single FASTA file as input for geNomad, while retaining the original bin assignment information in the contig headers. You can check the names of the first 10 contigs using:

```bash
grep '^>' 07_PHAGES/all_hq_bins.fna 
```

Write a batch job script to run geNomad following the script below.

Note that geNomad is not available in Roihu modules so we have installed the software for you. Therefore you need to __call the tool from the project applications__ with:

```bash
export PATH="/projappl/project_2022405/genomad-1.13.0/bin:$PATH"
```

geNomad creates a lot intermediate files, so we need 100G [flash storage](https://docs.csc.fi/supercomputing/disk/#temporary-local-disk-areas) and a cleanup step in the command. geNomad searches a big database of protein profiles which takes up a lot of memory. We can split the search into chunks using the `--splits` parameter to prevent the execution form failing due to insufficient memory.

```bash
#!/bin/bash
#SBATCH --job-name=genomad
#SBATCH --output=00_LOGS/%x-%j.out
#SBATCH --error=00_LOGS/%x-%j.err
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --account=project_project_2022405
#SBATCH --bb="#BB_LUA SBF storagesize=100G path=/sbf/data"

export PATH="/projappl/project_2022405/genomad-1.13.0/bin:$PATH"
DB=/scratch/project_2022405/DATABASES/genomad_db-1.9

# Add MAGs fasta here
MAGS="07_PHAGES/hq_bins.fna"

OUTDIR="07_PHAGES/GENOMAD"
mkdir -p "$OUTDIR"

genomad end-to-end \
    --cleanup \
    --splits 16 \
    --threads "$SLURM_CPUS_PER_TASK" \
    "$MAGS" \
    "$OUTDIR" \
    "$DB"
```
Check the used options (and others) by `genomad -h` and `genomad end-to-end -h`. (Run `export PATH="/projappl/project_2022405/genomad-1.13.0/bin:$PATH"` before running the command)

* What does end-to-end mode do?
* What kind of information can you expect as an output?

#### geNomad output 

Explore the output. Have a look at the log files first:

* What steps did geNomad run?
* How many viral contigs were identified in the high-quality MAGs? What about plasmids?

You can look at output with:

```bash
less -S file.tsv 	#-S disables line wrapping
```
Press `q` to quit

OR with `nano`.

Count predictions in *virus_summary.tsv output:
```bash
tail -n +2 file.tsv | wc -l
```

Find summary table, where viral contigs are listed:
* What viral taxa were predicted?
* What length do viral contigs have?
* Are there proviruses predicted?

### Quality control with CheckV

We will assess the quality and completeness of viral contigs identified by geNomad with [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/). 

* How does CheckV predict the quality? 
Check from the publication [Nayfach, S., Camargo, A.P., Schulz, F. et al. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nat Biotechnol 39, 578–585 (2021)](https://www.nature.com/articles/s41587-020-00774-7)

Note that CheckV needs its database, which is already downloaded to `/scratch/project_2022405/DATABASES/`

Run CheckV interactively or using a batch script following the script below.

You need: 8G memory and 8 CPU. Running CheckV will only take a few minutes.

```bash
#!/bin/bash
#SBATCH --job-name=checkv
#SBATCH --output=00_LOGS/%x-%j.out
#SBATCH --error=00_LOGS/%x-%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --account=project_project_2022405

export PATH="/projappl/project_2022405/checkVv1.1.1/bin:$PATH"
DB=/scratch/project_2022405/DATABASES/checkv-db-v1.5

# Add geNomad virus.fna here
VIRAL_CONTIGS="07_PHAGES/GENOMAD..."

OUTDIR="07_PHAGES/CHECKV"
mkdir -p "$OUTDIR"

checkv end_to_end \
    "$VIRAL_CONTIGS" \
    "$OUTDIR" \
    -d "$DB" \
    -t "$SLURM_CPUS_PER_TASK"
```

Check the options you used with `checkv -h` and `checkv end_to_end -h`.

#### CheckV output

Explore the output, especially the summary file quality._summary.tsv:

* Are there any proviruses predicted? Are those contigs that were flagged as proviruses by geNomad now also flagged as proviruses by CheckV? Any new proviruses compared to geNomad predictions?*
* Are most contigs of low, medium or high quality? See how the quality corresponds to completeness.
* Are there any 100% complete viral genomes listed?
* Any warnings?

Check for proportions of quality predictions (you can ask AI to explain command):

```bash
awk -F'\t' 'NR>1 && $8!="" {count[$8]++; total++}
END {
    print "Quality\tCount\tPercentage"
    for (tier in count)
        printf "%s\t%d\t%.2f%%\n", tier, count[tier], 100*count[tier]/total
}' path-to/quality_summary.tsv
```

*We used the geNomad `.fna` output as input for CheckV. However, geNomad has already excised the identified proviral sequences from their host-derived flanking regions and renamed the contigs to include provirus coordinates (e.g. `contig_1480|provirus_18530_57545`). Since the host-derived regions have already been removed, CheckV may not recognize some sequences as proviruses. In other cases, it may recommend further trimming or identify additional proviruses that geNomad did not detect. CheckV produces a proviruses.fna file containing sequences trimmed according to its recommendations, with updated contig headers and coordinates (e.g. `NODE_1301_length_66349_cov_37.879129_1 3260-66349/66349`). For more accurate proviral length, gene count, viral gene and host gene estimates, this file should be used as input for a second CheckV run. Although CheckV may recommend further trimming, running it more than twice is generally not standard practice.


In a real project, __CheckV output is typically used for filtering__ some predictions out. Generally, high- (>90% completeness) and medium-confidence (50-90% completeness) estimates are quite accurate. Common thresholds for metagenomic viral contigs include:
* at least 1 virus gene identified by CheckV;
* host to virus gene count ratio no more than 1:1;
* length minimum of 5 kbp or 10 kbp, unless a genome is >=50% complete (but not shorter than 1 kbp anyway).


In this course, we won't filter any viral contigs to keep things simple enough.


## Strain engraftment with inStrain

The next step is to determine the strain engraftment of the donor MAGs to few selected recipients.  
Before we can run the workflow, we need to fetch the recipient data and map it against the donor MAGs.  

First, make a whole new folder for the strain engraftment part. And put all output files and files we make in this folder.

```bash
mkdir -p 08_RECIPIENTS
cd 08_RECIPIENTS
```

### Fetch recipient data

We will download all the data from three recipients. Two from the FMT and one from the placebo group.  
The subject identifiers for these are:

```bash
TF38
TF29
TF45
```

Find all accession numbers for these subjects from SRA and write them down.  
You should have 15 read accessions to download.  

Put the recipient data inside a `Data` folder in the `08_RECIPIENTS` folder. Make sure you store them as compressed files (`.gz`). And it might be a good idea to make an array batch job for this task.  

```bash
mkdir Data
```

Make an array job that downloads the read files for each accession with `fasterq-dump` and also compresses them with `pigz`. Use the metaphlan script as an example and the same commands we used to download the donor metagenomes.  

```bash
sbatch src/YOUR_ARRAY_SCRIPT.sh
```

And if there were errors or problems with the download, you can copy the files from the course data folder (`/scratch/project_2022405/GutBugsData/Rawdata`).

### Mapping of recipient reads to donor assemblies

We will use the donor assembly as the reference for mapping the recipient reads. We can use the same Bowtie index we created for the donor reads in the binning step, located in /scratch/project_2022405/$USER/MMB-901_Metagenomics/04_MAPPING/.

To map the recipient reads to the donor assembly, we will use a modified version of the mapping script we used for the donor reads. For this you will need to a file with all the recipient accessions, similar to the one we made for the donor.  

Create the files and put them in the `07_RECIPIENTS` folder. Then modify the mapping script to read the recipient accessions and map them to the donor assembly.
After that, you can submit the job(s) and wait for them to finish.

### inStrain

[inStrain](https://instrain.readthedocs.io/en/latest/index.html) is a tool that uses a reference set of genomes or our assemblies and Bowtie2 mapping results to identify unique strains in a sample. The previous mapping results of the recipients against the donor assemblies will be used to identify possible identical strains in the participants.

We will first need to run progical, which predicts genes for us in the assembly. Modify the script below and run interactively/use a batch job with 20G of memory for 20 min with 8 CPUs.

```bash
export PATH="/projappl/project_2022405/instrain-v1.12.0/bin:$PATH"
prodigal -i assembly.fasta -d genes.fna -a genes.faa
```

inStrain has been installed for you in projappl.

Check the [tutorial](https://instrain.readthedocs.io/en/latest/tutorial.html) and create an inStrain profile command using the scaffold file representative_genomes.stb we created using the dRep helper script.

Modify the script below to include the correct BAM files for each recipient mapping, each donor mapping, and the DF16 assembly contigs. This can be done in an array script with 50 GB of memory, 5 h of time, and 6 CPUs.

```bash
export PATH="/projappl/project_2022405/instrain-v1.12.0/bin:$PATH"
inStrain profile mappingfile.bam DF16_contigs.fasta -o outputlocation.IS -p 6 -g genes.faa -s scaffoldtobin.stb --database_mode
```
