#!/bin/sh
#SBATCH -J download_rnaseq

cd /data/westlab/jes157/data_input/fastq/rnaseq/

for url in $(cat rnaseq_zic_kd_urls.txt);
do 
 wget $url
done

