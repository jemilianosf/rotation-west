#!/bin/sh
#SBATCH -J download_chip

cd /data/westlab/jes157/data_input/fastq/chip_zic/

for url in $(cat chip_zic_cerebellum_sample_urls.txt);
do 
 wget $url
done

