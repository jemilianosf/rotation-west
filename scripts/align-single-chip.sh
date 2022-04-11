#!/bin/bash
# Author: VJ Ramesh, Adapted from Melyssa Minto Script 
# Date: 7/21/2021
# Updated: 7/21/2021
# This script will align ChIP seq single-end fastqs

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
WD='/data/westlab/jes157/data_output/chip/' # where the alignments will be saved
QC='/data/westlab/jes157/data_output/chip/QC/'                 # where the fastQC results will be saved
TRIMMED='/data/westlab/jes157/trimmed/' # where the trimmed fastqa will be saved
trim='/data/westlab/genomeTools/Trimmomatic-0.38/trimmomatic-0.38.jar'
adapters='/data/westlab/genomeTools/Trimmomatic-0.38/adapters/'
REF='/data/westlab/reference/gencode_mouse_vM21/STAR_mm10_2.7.2b/'
GTF='/data/westlab/reference/gencode_mouse_vM21/gencode.vM21.annotation.gtf'

###1. start aligment
echo .............$NAME..............

## get fastqc
cd $RD
# unzip file
echo unzipping fastq
FILE=$NAME'.fastq.gz'
if [ -e $NAME'.fastq.gz' ]
then
 gunzip $NAME'.fastq.gz'
fi


##2. FastQCc
echo fastQC...
fastqc $NAME'.fastq' --outdir=$QC

##3. Trim adaptors and output to a specified folder. The name of the trimmed file is $NAME + _trimmed.fastq.gz
cd $TRIMMED
echo trimming adaptors...
java -jar $trim SE -threads 8 -phred33 $RD'/'$NAME'.fastq' $NAME'.trimmed.fastq' ILLUMINACLIP:${adapters}TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


#sorting trimmed reads
echo sorting trimmed fastqs
cat $TRIMMED$NAME'.trimmed.fastq' | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" >$TRIMMED$NAME'.sorted.fastq'


## Run alignment
cd $WD     # changing directory to where the alignments should be saved
mkdir $NAME # making a folder to hold the samples alignment files
cd $NAME        # changing directory into the folder created above
echo run alignment...
STAR --genomeDir $REF --runThreadN 8 --readFilesIn $TRIMMED$NAME'.sorted.fastq' --outFileNamePrefix $NAME --limitBAMsortRAM 20000000000 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard


# Filter out reads that have aligned
#Flags:
#-F #: Do not output alignments with any bits set in # present in the FLAG field
#-b:Output in BAM format
#-h:Include the header in the output
#-o: ouput to file
echo getting accepted hits...
samtools view -F 4 -b -h -o "$NAME"accepted_hits.bam *Aligned.sortedByCoord.out.bam || { echo 'ERROR: could not remove unmapped reads from .bam alignment, exiting...' ; exit 1; }

echo beginning bamcoverage normalization...
# index bamfile for bamcoverage 
samtools index -b "$NAME"accepted_hits.bam "$NAME"accepted_hits.bam.bai

##### This next section is for ATAC, ChIP, or even DNAse data ############
## Calling peaks

#Big to Bigwig
bamCoverage \
-b *accepted_hits.bam \
-p 8  \
-o $NAME'_norm.bw' \
--normalizeUsing BPM \
--effectiveGenomeSize 2730871774 \ 
--ignoreForNormalization chrX 2>&1 | tee $WD$NAME'/bamcoverage.log'


#Converting BAM to BED
#Flags:
#-i: input
echo convert bam to bed...
bedtools bamtobed -i "$NAME"accepted_hits.bam > sequence.bed || { echo 'ERROR: could not produce sequence.bed from .bam alignment, exiting...' ; exit 1; }

#Sorting feature file by Chromosome
echo 'Starting filtering...'
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
echo fitering duplicates
macs2 filterdup -i sorted.all.bed --keep-dup=1 -o sequence.final.bed

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
macs2 callpeak -t sequence.final.bed -f BED -g mm -n ./macs_n.fdr01 -q 0.01 --bdg --nomodel --ext 147 || { echo ERROR:macs2 does not work, exiting... ; exit 1; }


echo "................10.  Finished Peak Calling, Now Tidying Up"


rm *sortedByName.bam
rm *chr*
rm ${TRIMMED}${NAME}*

# now that we are done with everything we can zip up the fastq file to save space
cd $RD
gzip -5 $NAME'.fastq'

# FIN!
echo "pipeline is done!"
