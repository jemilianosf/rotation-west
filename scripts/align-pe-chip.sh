#!/bin/bash
# Author: VJ Ramesh, Adapted from Melyssa Minto Script 
# Date: 7/21/2021
# Updated: 7/21/2021
# This script will align ChIP seq paired-end fastqs

# usage: ./align.sh sample /absolute/path/to/fastq
# programs needed to run this: fastqc, STAR, UCSC genome tools aka Kent Utils, samtools, macs2, htseq, trimmomatic

# Check if the user input arguments <sample name> and <data directory>
if [ -z "$1" ]; then
 echo usage: $0 '<name> <datadir>'
 exit
elif [ -z "$2" ]; then
 echo 'no data dir'
 exit
fi

# set up all of my variables
NAME=$1 #<name>
RD=$2 #<datadir>
NAME1=${NAME}"_1"
NAME2=${NAME}"_2"
WD='/data/westlab/jes157/data_output/chip/' # where the alignments will be saved
QC='/data/westlab/jes157/data_output/chip/QC/'                 # where the fastQC results will be saved
TRIMMED='/data/westlab/jes157/trimmed/' # where the trimmed fastqa will be saved
trim='/data/westlab/genomeTools/Trimmomatic-0.38/trimmomatic-0.38.jar'
adapters='/data/westlab/genomeTools/Trimmomatic-0.38/adapters/'
REF='/data/westlab/reference/gencode_mouse_vM21/STAR_mm10_2.7.2b/'
GTF='/data/westlab/reference/gencode_mouse_vM21/gencode.vM21.annotation.gtf'

#begin pipeline
echo ".............$NAME.............."

echo "making sample directory..."
mkdir $WD$NAME


#1. unzipping files
echo "unzipping fastq file..."
cd $RD
# unziping read1
if [ -e $NAME1'.fastq.gz' ]
then
 gunzip $NAME1'.fastq.gz'
fi
# unzipping read2
if [ -e $NAME2'.fastq.gz' ]
then
 gunzip $NAME2'.fastq.gz'
fi

#2.  Get fastQC
echo "..........2. FastQC................"

fastqc $NAME1'.fastq' --outdir=$QC
fastqc $NAME2'.fastq' --outdir=$QC

#3. Trim adaptors 
echo "..........3. Trimming Adaptors......"
java -jar $trim PE \
-threads 2 \
-phred33 \
$RD'/'$NAME1'.fastq' $RD'/'$NAME2'.fastq' $TRIMMED$NAME1'.paired.fastq' $TRIMMED$NAME1'.unpaired.fastq' $TRIMMED$NAME2'.paired.fastq' $TRIMMED$NAME2'.unpaired.fastq' \
ILLUMINACLIP:$adapters'/'TruSeq3-PE-NEB.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36 2>&1 | tee $WD$NAME'/trimmonatic.log'

#4. sorting trimmed reads 
echo "...........4. Sorting Trimmed Fastqs..."
cat $TRIMMED$NAME1'.paired.fastq' | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" >$TRIMMED$NAME1'.sorted.fastq'
cat $TRIMMED$NAME2'.paired.fastq' | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" >$TRIMMED$NAME2'.sorted.fastq'

#5. Align to mm10#

echo ".............6.  Aligning Sorted Fastq to MM10 genome"

cd $WD$NAME
STAR \
--genomeDir $REF \
--runThreadN 8 \
--readFilesIn $TRIMMED$NAME1'.sorted.fastq' $TRIMMED$NAME2'.sorted.fastq' \
--outFileNamePrefix $NAME \
--limitBAMsortRAM 20000000000 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard

#7.creating normalized non-stranded bigwig with deeptools

echo ".............7.  Beginning Bamcoverage ..."


#index bamfile for bamcoverage 
samtools index -b "$NAME"Aligned.sortedByCoord.out.bam "$NAME"Aligned.sortedByCoord.out.bam.bai



bamCoverage \
-b "$NAME"Aligned.sortedByCoord.out.bam \
-p 8  \
-o $NAME'_norm.bw' \
--normalizeUsing BPM \
--effectiveGenomeSize 2654621837 \
--ignoreForNormalization chrX 2>&1 | tee $WD$NAME'/bamcoverage.log'


#8.BAM to Bed Conversion
echo ".............8.  Convert sortbyCoord.bam to Bed..."

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
rm sequence.bed sorted.all.bed *chr*


# peak calling
#Flags:
#-t: treatment filename
#-g: genome size ( hs: 2.7e9, mm: 1.87e9, ce: 9e7, dm:1.2e8)
#-n: name of experiment
#-q: The qvalue (minimum FDR) cutoff to call significant regions. Default is 0.01
#--nomodel:While on, MACS will bypass building the shifting model
#--shift --100:
#--ext:extend reads in 5’->3’ direction to fix-sized fragments.
    #1. To find enriched cutting sites such as some DNAse-Seq datasets. In this case, all 5’ ends of sequenced reads should be extended in both direction to smooth the pileup signals. If the wanted smoothing window is 200bps, then use ‘–nomodel –shift -100 –extsize 200’.
    #2. For certain nucleosome-seq data, we need to pileup the centers of nucleosomes using a half-nucleosome size for wavelet analysis (e.g. NPS algorithm). Since the DNA wrapped on nucleosome is about 147bps, this option can be used: ‘–nomodel –shift 37 –extsize 73’.
#--bdg:MACS will store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files

# For ChIP seq, no extra arguments needed
# For DNAse seq use these arguments: -nomodel --shift -100 --extsize 200
# For ATAC seq use these arguments: -nomodel --shift 37 --extsize 

echo Calling ChIP peaks with MACS...
macs2 callpeak -t sequence.final.bed -f BED -g mm -n ./${NAME}_macs_n.fdr01 -q 0.01 --bdg --nomodel --ext 147 || { echo ERROR:macs2 does not work, exiting... ; exit 1; }

rm sequence.final.bed 

echo "................10.  Finished Peak Calling, Now Tidying Up"


rm *Aligned.sortedByCoord.out.bam
rm *chr*
rm ${TRIMMED}${NAME}*

# now that we are done with everything we can zip up the fastq file to save space
cd $RD
gzip -5 $NAME1'.fastq'
gzip -5 $NAME2'.fastq'

# FIN!
echo "pipeline is done!"
