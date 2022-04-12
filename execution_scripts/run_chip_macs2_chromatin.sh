#!/bin/sh
#SBATCH -J ChIP
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH -a 1-17%4

module load gcc
module load STAR
module load fastqc
module load bedtools2
module load samtools
module load deepTools/2.0-gcb01
module load macs2/2.1.0.20140616-fasrc01

names_file=/data/westlab/jes157/data_input/fastq/chip_chromatin/chip_chromatin_ids.txt
READ1=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' $names_file)

/data/westlab/jes157/scripts/chip/align-single-chip.sh $READ1 /data/westlab/jes157/data_input/fastq/chip_chromatin