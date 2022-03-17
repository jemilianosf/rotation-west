#!/bin/sh
#SBATCH -J macs
#SBATCH --mem=8G
#SBATCH -a 1-12
module load bedtools2

module load macs2/2.1.0.20140616-fasrc01

names_file=/data/westlab/jes157/data_input/fastq/cutnrun_zic_dirs.txt
peaks_dir=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' $names_file)
peaks_name=${peaks_dir##*/}
cd $peaks_dir

/data/westlab/jes157/scripts/bam2bed.sh $peaks_name

#8. Peak Calling Using SEACR
echo 'Calling peaks with MACS...'
macs2 callpeak -t sequence.final.bed -f BED -g mm -n ${peaks_dir##*/}.macs_n.fdr01 -q 0.01 --bdg || { echo ERROR:macs2 does not work, exiting... ; exit 1; }
