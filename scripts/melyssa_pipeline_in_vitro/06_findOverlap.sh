#!/bin/bash
# Author: Melyssa Minto
# This script will use bedtools to find the overlap in peaks (Zic, DNase, H3K27ac) in the anchors of the 
# predicted chroamtin loops in adult cerebellum

cd ../data_output/diffbind_cutnrun_zic/melyssa_pipeline

P56ANCHOR=combined_MAPS_peaks.txt
P4ANCHOR=P4_loops.bed

cat $P56ANCHOR | awk -v OFS='\t' '{ print $1, $2, $3 }' > anchors.bed
cat $P56ANCHOR | awk -v OFS='\t' '{ print $4, $5, $6 }' >> anchors.bed
cat $P4ANCHOR | awk -v OFS='\t' '{ print $1, $2, $3 }' >> anchors.bed

bedtools intersect -wa -wb -a anchors.bed -b ~/repos/rotation-west/data_output/bed/cutnrun_zic_seacr_peaks_union.bed > zic_anchors.bed
