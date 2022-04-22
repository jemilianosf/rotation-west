#!/bin/sh
#SBATCH -J download_chip
#SBATCH -c 4
#SBATCH --mem-per-cpu=4G


cd /data/westlab/jes157/data_input/fastq/chip_chromatin/

for fq in $(cat chip_chromatin_pe_ftp.txt.txt);
do 
  wget $fq
done
