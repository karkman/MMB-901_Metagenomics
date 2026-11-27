#!/bin/bash -l
#SBATCH --job-name bowtie2_index
#SBATCH --output 00_LOGS/bowtie2_index_%j.out
#SBATCH --error 00_LOGS/bowtie2_index_%j.err
#SBATCH --time 02:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 50G
#SBATCH --account project_XXXXXXX

module load bio-apps/v202603
module load bowtie2/2.5.4

mkdir -p 04_MAPPING

bowtie2-build \
    --threads $SLURM_CPUS_PER_TASK \
    02_ASSEMBLY/contigs.fasta \
    04_MAPPING/contigs
