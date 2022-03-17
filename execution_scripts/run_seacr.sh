#!/bin/sh
#SBATCH -J seacr
#SBATCH --mem=8G
#SBATCH -a 1-8


names_file=/data/westlab/jes157/data_input/fastq/chip_zic_dirs.txt
peaks_dir=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' $names_file)


cd $peaks_dir
module load Anaconda3/2019.10-gcb02
source activate /data/westlab/conda

module load bedtools2

chromSize='/data/westlab/reference/gencode_mouse_vM21/STAR_mm10_2.7.2b/chrNameLength.txt'
seacr='/data/westlab/genomeTools/SEACR/SEACR_1.3.sh'

bedtools genomecov -bg -scale 1.0 -i sequence.final.bed -g $chromSize > sequence.normalized.bedgraph 
bash $seacr sequence.normalized.bedgraph 0.01 non stringent ${peaks_dir##*/}.seacr_peaks.bed
