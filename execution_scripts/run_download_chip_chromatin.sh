#!/bin/sh
#SBATCH -J download_chip
#SBATCH -c 4
#SBATCH --mem-per-cpu=4G

module load sratoolkit/2.9.6-gcb01

cd /data/westlab/jes157/data_input/fastq/chip_chromatin/

for id in $(cat chip_chromatin_ids.txt);
do 
 fastq-dump --gzip $id
done
