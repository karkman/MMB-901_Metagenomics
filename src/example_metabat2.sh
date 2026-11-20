#!/bin/bash -l
#SBATCH --job-name=metabat2
#SBATCH --account=project_XXXXXXX
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=00_LOGS/metabat2_%j.out

module load bio-apps/v202603
module load metabat/2.15

mkdir -p 00_LOGS
mkdir -p 05_BINNING/MetaBAT2

# Generate contig coverage file
jgi_summarize_bam_contig_depths \
    --outputDepth 05_BINNING/MetaBAT2/depth.txt \
    04_MAPPING/*.sorted.bam

# Run MetaBAT2
metabat2 \
    -i 02_ASSEMBLY/DF16/contigs.fasta \
    -a 05_BINNING/MetaBAT2/depth.txt \
    -o 05_BINNING/MetaBAT2/bin \
    -t $SLURM_CPUS_PER_TASK
