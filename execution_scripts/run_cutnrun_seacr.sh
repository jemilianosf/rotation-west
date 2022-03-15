#!/bin/sh
#SBATCH -J CutAndRun
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH -a 1-12%4

module load gcc
module load STAR
module load fastqc
module load bedtools2
module load samtools
module load deepTools/2.0-gcb01
module load macs2/2.1.0.20140616-fasrc01

names_file=/data/westlab/jes157/data_input/fastq/cutnrun_zic_data.txt
READ1=$(awk -v i=$SLURM_ARRAY_TASK_ID 'BEGIN{FS="\t"} NR==i {print $1}' $names_file)
READ2=$(awk -v i=$SLURM_ARRAY_TASK_ID 'BEGIN{FS="\t"} NR==i {print $2}' $names_file)

/data/westlab/jes157/scripts/cutnrun/cutnrun-seacr.sh $READ1 $READ2 /data/westlab/jes157/data_input/fastq/cutnrun_zic