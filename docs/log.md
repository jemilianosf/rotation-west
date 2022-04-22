# Log

#2022-03-08

- Read background papers

- Ideas:
- looking at enhancer variation associated with autism 
- are fast evolving regulatory sequences related to schizophrenia / autism involved in neuron maturation rather than cell identity?
- How to assign genes to binding sites?

- Look through / work through VJ scripts

#2022-03-09
- Read background papers
- Look through / work through VJ scripts with Zic data
- Activity dependent transcription Seminar

#2022-03-10
- Activity dependent transcription Seminar
- Run cutnrun pipeline 

#2022-03-15
- How to do diff binding between timepoints
- Process ChIP data
- Read grant
- Overlap processed ChIP and CutNRun

#2022-03-16
- Why is ChIP pipeline not working on some samples
It was still due to an earlier version of the pipeline.
- Run ChIP pipeline
- Run SEACR for all samples
- CutNRun DESEq
- Meetings

#2022-03-17
- Why is seacr  not working on bash scripts?
- Run MAC2 for CutNRun after exploring in the browser SEACR, MACS and MACS2 for ChIP, and macs2 really looks better than either seacr or previous call
- Made new analysis: results/zic_cutnrun_chip_reproccessed_overlap.Rmd

TL;DR:
From looking at the browser
- New macs2 chip peaks look better than previous ones on the browser. They are shorter, and more centered on the signal summit.
- seacr chip peaks don't look as good as macs2 peaks. 
- Macs2 seems to work well for cutnrun as well. Calls less peaks where signal looks more noisy compared to seacr, but also misses some peaks that look real that seacr doesn't.

From looking at the overlaps
- When comparing macs2 vs seacr, most macs2 peaks are also called by seacr, but seacr calls many more. These are likely a mix of seacr false positives, and macs2 false negatives.
- If we want to compare cutnrun vs chip with the same caller, I would trust more the macs2 comparison, which works reasonably well on both datasets.


- CutNRun DESEq -  test per each combination - results/zic_cutnrun_deseq.Rmd
- Think  more about the cutnrun data: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1802-4

#2022-03-18
- Finish overlaps analysis
  - Sort by overlap size
  - Change colors
- DESeq with seacr peaks 


#2022-03-21
- Annotation of DESeq2 diff binding seacr peaks in vitro in results/zic_cutnrun_deseq.Rmd
- Started DESeq2 chip analysis

- Left to do: general ChIPseeker of all peaks (seacr and macs2)

#2022-03-22
- Overlaps of DESeq2 diff peaks and chip in vivo peaks
- ChIPseeker of all peaks
- Changed the overlap thing to make symmetric overlaps when comparing by element.
- Put together slides

#2022-03-23
- Add number of peaks to slides
- Get library stats
- Meeting

#2022-03-24
- RNAseq overlap with diff peaks
 -> only a few diff bound genes are affected by 
DESeq2 ChIPseq + overlap of diff peaks
 -> There is overlap between changing peaks in vivo + in culture 
 
#2022-03-25
- Re-organize DESeq script
- Rank promoters by peak lfc and visualize top peaks
- Assigning peaks to genes with list of cerebellum mapped enhancer - gene pairs
- Assigning peaks to genes with HiC
- Overlap with H3K27me3 peaks
- Correlation of tracks
- Density plot around tss, around peaks

#2022-03-28
- Added split markdowns of gene overlaps into overlaps and overlaps kdrna markdowns. 
- PCA with chip dataset results/diffbind_step7_pca_chip_vs_cutnrun.Rmd.
- MACS paper

#2022-03-29
To do:
- Do suggestions on overlaps kdrna analysis:
  - Focus on genes that change, where is the nearest diff peak? where is the nearest peak?
- Do enhancer-promoter peak assignments with P4 hic data.
- Check how to use a more stringent SEACR peak list.

#2022-03-30
- Violin plots for knockdown genes with diff peaks
- Continue with enhancer-promoter peak assignments

#2022-03-31
- Continue with enhancer-promoter peak assignments

#2022-04-01
- Get table and annotated peaks
- General exploration of peaks / enhancers overlaps
- Re-do previous overlaps but now with enhancers

#2022-04-04
- Do the 3 vs 7 day comparison, add that to existing plots/analyses
  - I did that comparison, it looks similar to 7 vs 1 results, although not entirely
  - I realize the other comparison was not 7vs3 but KD vs 7, so I guess the plots I already have with 7 vs 1 are fine
- Do example plots from DESeq2 counts

#2022-04-05
- Do the peak - gene assignments with P22 loops:
  - Do assignments
  - General exploration - just add the new loops to the assignment
  - Re - do previous overlaps but now with P22 enhancers
  
- Check methods of both previous papers

#2022-04-06
- Meeting:
To do:
 - Re-do the enhancer assignments using the whole bin / Melyssa's strategy
 - Dig more into the RNA seq knockdown data
 - Compare HiC p22 vs p4
 
#2022-04-07
- Added correlation Zic vs Chd7 knockdown rna seq analysis

#2022-04-08
- Re-do the enhancer assignments using the whole bin / Melyssa's strategy

#2022-04-11
To do for the week:
- Re-do the RNAseq in vitro + correlation with Chd7 / Chd4 cKOs
  - Got count matrices

- Re-align the ChIP seq data:
  Reddy et al: P4 CTCF and Chd7, both WT and Chd7cKO.
  Goodman et al: the P22 CTCF and Chd4, both WT and Chd4cKO,
  - Download data / went over the scripts w Melyssa 

- Re-do the enhancer assignments using the whole bin / Melyssa's strategy
 
#2022-04-12
- Re-align the ChIP seq data:
  - Run realignments
- Re-do the RNAseq in vitro + correlation with Chd7 / Chd4 cKOs
  - DESEq 
  - Correlations
  
#2022-04-12
- Re-do the RNAseq in vitro + correlation with Chd7 / Chd4 cKOs
  - Correlations

- RNAseq Zic KD exploration of targets
Main points:

- Zic is a negative transcriptional regulator of Chd7: 
- Melyssa showed Zic binding in Chd7 enhancers 
- In vitro CutNRun shows Zic binding in Chd7 promoter
- Chd7 expression goes up in both Zic1 and Zic2 knockdowns
- Zic2 expression does not change in Chd7 
- Zic1 expression does not change in Chd7 ??? 

- Chd7 cKO and Zic1/2 KD show overlap of disregulated genes
- In particular in Zic1 KD, negative correlation in fold change with Chd7 cKO

- Zic and Chd7 might also overlap at their target enhancers:
- Different direction of change for genes that are significant in both 
- Overlap of Zic and Chd7 ChIP???


#2022-04-13

- ChIP seq overlaps
  - Clean the hardac directories / copy all that's left
  - Create session in IGV to show overlaps
  - Do a quick overlap of peaks

- Re-assign enhancer-promoters peaks with Melyssa's pipeline

#2022-04-14
- Re-assign enhancer-promoters peaks with Melyssa's pipeline
- Process the new ChIPs
- Look into the rank rank hypergeometric overlap
- Re-analyze Chd7/Chd4 rnaseq 



- Other things that came up:
  - Look at seacr signal to get more stringent peaks
  - Compare gene-lists by rank rank hypergeometric overlap
  - Compare Chd7/4 knockdowns with time changes
  - Look at the peaks that go up both in vivo and in vitro and what they look like in expression.
  
#2022-04-15
- Going through:  Re-assign enhancer-promoters peaks with Melyssa's pipeline script

#2022-04-18
- Check new chip tracks

- Do overlaps with new chip tracks and also with assigned files
-  how to do rank rank test

- read 1 paper

#2022-04-19
- Do overlaps plot of enhancers assignments
- Do RRHOs
- Explore bowtie2 output:
- Done - So the issue was PE files as SE

#2022-04-20
- Run:
/data/westlab/jes157/data_output/chip/chip_chromatin_pe

- read 2 paper

#2022-04-21

- CHD KO vs Zic KD RNA SEQ comparison
  - Re-analyze rnaseq for chd7
  - 
- Read Chd paper 

#2022-04-22
- Worked on direct target analyses.


