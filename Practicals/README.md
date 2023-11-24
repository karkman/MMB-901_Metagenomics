# MMB-901 Microbial metagenomics – Practicals

## Introduction

On this course we will use a subset from a publicly available human gut metagenomic sequencing data from fecal microbiota transplantation (FMT) study: [Wilson, B.C., Vatanen, T., Jayasinghe, T.N. _et al._ Strain engraftment competition and functional augmentation in a multi-donor fecal microbiota transplantation trial for obesity. Microbiome 9, 107 (2021)](https://doi.org/10.1186/s40168-021-01060-7).  

Using read-based and assembly-based approaches, we will study the human gut microbiome. With read-based approach we will determine the taxonomic composition of one donor microbiome (female super-donor DF16) and compare to those the taxonomic composition of the recipient microbiomes taken at baseline and at 6, 12 and 26 weeks after the FMT. With assembly based approach, we will use [anvi'o](anvio.org) to construct metagenome-assembled genomes (MAGs) from the donor samples. Then we will assees the strain engrafment of selected strains in a selected set of recipients.  

## Setup



## Data 

First we will retrieve all the metagenomic sequence files for the female super donor, DF16. Go to the publication and find the sequencing project repository accession (BioProject accession).  
Then go to [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) and find all the reads for the project.  Open them in Run selector and find the run accessions of the samples from the super-donor (9 sequencing runs).  
Finallly download only the accession numbers. We will use that file in the next step to download the data to Puhti.  

We will use [Kingfisher](https://wwood.github.io/kingfisher-download/) to download the read files. It has been installed to Puhti under our course project application folder.  
But before downloading, have a look at the documentation and find out how to use it. We will use the `ena-ftp` method for downloading the reads.  


```bash

~/projappl/Kingfisher/bin/kingfisher get \
    --run-identifiers-list Recipient_accession.txt \
    -m ena-ftp \
    --download-threads $SLURM_CPUS_PER_TASK \
    --extraction-threads $SLURM_CPUS_PER_TASK
```

# Quality control

## Metagenome assembly

## Read-based taxonomy

## Assembly QC

## Genome-resolved metagenomics

## MAG QC

## Strain engraftment

