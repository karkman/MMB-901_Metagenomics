#!/bin/bash -l
#SBATCH --job-name spades
#SBATCH --output 00_LOGS/spades_%j.out
#SBATCH --error 00_LOGS/spades_%j.err
#SBATCH --time 16:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 12
#SBATCH --mem 150G
#SBATCH --account project_XXXXXXX
#SBATCH --gres=nvme:500

cd /scratch/project_XXXXXXX/$USER/MMB-901_Metagenomics
mkdir -p 00_LOGS

module load spades/3.15.5

spades.py \
        -1 01_DATA/DF16_1.fastq.gz \
        -2 01_DATA/DF16_2.fastq.gz \
        --threads $SLURM_CPUS_PER_TASK \
        --meta \
        --only-assembler \
        --memory 160 \
        --tmp-dir $LOCAL_SCRATCH \
        -o 02_ASSEMBLY
