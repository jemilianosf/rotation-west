#!/bin/bash

if [ -z "$1" ]; then
 echo usage: $0 '<name>'
 exit
fi

NAME=$1

#8.BAM to Bed Conversion
echo ".............Convert sortbyCoord.bam to Bed..."

bedtools bamtobed -i "$NAME"Aligned.sortedByCoord.out.bam > sequence.bed || { echo 'ERROR: could not produce sequence.bed from .bam alignment, exiting...' ; exit 1; }

#Sorting feature file by Chromosome
echo 'Starting Filtering...'
CRMREFIX='chr'
    for (( i = 1 ; i <= 25; i++ ))
    do
    CRM=$CRMREFIX$i
    grep -w $CRM sequence.bed > "$CRM".bed  #output chromosome specific bed file
    bedtools sort -i "$CRM".bed > sorted."$CRM".bed
    done
grep -w "chrX\|chrY" sequence.bed > chrXY.bed
bedtools sort -i chrXY.bed > sorted.chrXY.bed
cat sorted.chr*bed > sorted.all.bed #combining all sorted feature files
grep -w "chrM" sequence.bed | wc -l > chrM.txt

# Filtering duplicates
echo "Filtering Duplicates"
macs2 filterdup -i sorted.all.bed --keep-dup=1 -o sequence.final.bed

rm *chr*


