#!/bin/sh
#SBATCH -J rnaseq
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH -a 1-6

module load gcc
module load STAR
module load fastqc
module load bedtools2
module load samtools
module load deepTools/2.0-gcb01
module load HTSeq/0.6.1-fasrc01

names_file=/data/westlab/jes157/data_input/fastq/rnaseq_names.txt
READ=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' $names_file)

/data/westlab/jes157/scripts/rnaseq/align-rnaseq.sh $READ /data/westlab/jes157/data_input/fastq/rnaseq