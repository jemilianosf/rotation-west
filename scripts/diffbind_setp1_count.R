## load cutnrun peaks as granges objects
library(rtracklayer)
library(tidyverse)

cutnrun_peaks_dir <- "../data_output/cutnrun_zic/seacr_peaks/"
cutnrun_peaks_files <- list.files(cutnrun_peaks_dir,pattern = ".bed")
cutnrun_peaks_full <- paste0(cutnrun_peaks_dir,cutnrun_peaks_files)


cutnrun_peaks_granges <- lapply(cutnrun_peaks_full, import.bedGraph)

names(cutnrun_peaks_granges) <- str_split(cutnrun_peaks_files,"_",n=2) %>% map_chr(1)

names(cutnrun_peaks_granges) <- paste0("cutnrun-",names(cutnrun_peaks_granges))

## get union
cutnrun_peaks_union_granges <- GenomicRanges::union(unlist(GRangesList(cutnrun_peaks_granges)),
                                                    unlist(GRangesList(cutnrun_peaks_granges)))
## remove peaks that overlap blacklist

mm10_blacklist <- import("../data_input/bed/mm10-blacklist.v2.bed.gz")
cutnrun_peaks_union_granges <- cutnrun_peaks_union_granges[!overlapsAny(cutnrun_peaks_union_granges, mm10_blacklist)]

# Export union
export.bed(cutnrun_peaks_union_granges,"../data_output/bed/cutnrun_zic_seacr_peaks_union.bed")


library(csaw)
bam_files <- list.files("../data_output/cutnrun_zic/bam/",pattern = ".bam$",full.names = T)
counts <- csaw::regionCounts(bam.files = bam_files,
                             regions = cutnrun_peaks_union_granges,
                             param = csaw::readParam(dedup = FALSE, minq = 255, pe = "both"))

# Get counts matrix
colnames(counts) <- str_split(basename(bam_files),"_",2) %>% map_chr(1)

colData(counts)$rep <- str_split(colnames(counts),"-") %>% map_chr(1)
colData(counts)$time <- str_split(colnames(counts),"-") %>% map_chr(2)

# Get DESEq object
library(DESeq2)

deseq_obj <- DESeqDataSet(counts, design = ~ time)


saveRDS(deseq_obj, "../data_output/diffbind_cutnrun_zic/deseq_obj.Rds")

