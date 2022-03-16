# !/bin/bash


peaks_dir=$1


echo "................Calling peaks with SEACR..."
cd $peaks_dir
module load bedtools2
module load Anaconda3/2019.10-gcb02
source activate /data/westlab/conda

chromSize='/data/westlab/reference/gencode_mouse_vM21/STAR_mm10_2.7.2b/chrNameLength.txt'
seacr='/data/westlab/genomeTools/SEACR/SEACR_1.3.sh'

bedtools genomecov -bg -scale 1.0 -i sequence.final.bed -g $chromSize > sequence.normalized.bedgraph 
bash $seacr sequence.normalized.bedgraph 0.01 non stringent ${peaks_dir##*/}.seacr_peaks.bed
