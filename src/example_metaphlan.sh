#!/bin/bash -l
#SBATCH --job-name metaphlan
#SBATCH --output 00_LOGS/metaphlan_%A_%a.out
#SBATCH --error 00_LOGS/metaphlan_%A_%a.err
#SBATCH --time 00:60:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 50G
#SBATCH --account project_XXXXXXX
#SBATCH --array=1-9
#SBATCH --gres=nvme:100

module load metaphlan/4.2.4
mkdir -p 00_LOGS
mkdir -p 05_TAXONOMY

SAMPLE_ACC=$(sed -n ${SLURM_ARRAY_TASK_ID}p 01_DATA/DF16_accessions.txt)

metaphlan \
    01_DATA/${SAMPLE_ACC}_1.fastq.gz,01_DATA/${SAMPLE_ACC}_2.fastq.gz \
    --nproc $SLURM_CPUS_PER_TASK \
    --sample_id ${SAMPLE_ACC} \
    --input_type fastq \
    --db_dir /scratch/project_XXXXXXX/DBs/metaphlan \
    --output 05_TAXONOMY/${SAMPLE_ACC}.txt \
    --mapout 05_TAXONOMY/${SAMPLE_ACC}.bowtie2.bz2 
