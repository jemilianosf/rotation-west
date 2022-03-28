# Gets raw coverage, using only uniquely mapping reads
for bam in $(ls *bam);
do 
  bamCoverage --extendReads --minMappingQuality 255 --normalizeUsing None --binSize 1 -b ${bam} -o ${bam%.*}_coverage.bw
done