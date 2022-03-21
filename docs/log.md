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


#2022-03-19
- Annotation of DESeq2 diff binding seacr peaks in vitro in results/zic_cutnrun_deseq.Rmd
- Started DESeq2 chip analysis

- Left to do: general ChIPseeker of all peaks (seacr and macs2)


#2022-03-19


