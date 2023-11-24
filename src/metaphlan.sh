#!/bin/bash -l
#SBATCH --job-name metaphlan
#SBATCH --output 00_LOGS/metaphlan_%A_%a.out
#SBATCH --error 00_LOGS/metaphlan_%A_%a.err
#SBATCH --time 00:30:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 12
#SBATCH --mem 20G
#SBATCH --account project_2009008
#SBATCH --array=1-9
#SBATCH --gres=nvme:100

module load metaphlan/4.0.6

SAMPLE_ACC=$(sed -n ${SLURM_ARRAY_TASK_ID}p 01_DATA/DF16_accessions.txt)

metaphlan \
    01_DATA/${SAMPLE_ACC}_1.fastq.gz,01_DATA/${SAMPLE_ACC}_2.fastq.gz \
    --nproc $SLURM_CPUS_PER_TASK \
    --unclassified_estimation \
    --sample_id ${SAMPLE_ACC} \
    --input_type fastq \
    --bowtie2db /scratch/project_2009008/DB/metaphlan \
    --output 05_TAXONOMY/${SAMPLE_ACC}.txt \
    --bowtie2out 05_TAXONOMY/${SAMPLE_ACC}.bowtie2.bz2 
