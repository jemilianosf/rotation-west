library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

deseq_obj_contrasts <- readRDS("../data_output/diffbind_cutnrun_zic/deseq_obj_contrasts.Rds")
res_7v1 <- readRDS("../data_output/diffbind_cutnrun_zic/res_7v1_ashr.Rds")

peaks_7v1 <- granges(rowRanges(deseq_obj_contrasts[tidyr::replace_na(res_7v1$padj,1) < 0.05 &
                                             tidyr::replace_na(res_7v1$log2FoldChange,0) > 0]))
peaks_1v7 <- granges(rowRanges(deseq_obj_contrasts[tidyr::replace_na(res_7v1$padj,1) < 0.05 &
                                             tidyr::replace_na(res_7v1$log2FoldChange,0) < 0]))
peaks_all <- granges(rowRanges(deseq_obj_contrasts))

elementMetadata(peaks_7v1) <- res_7v1[tidyr::replace_na(res_7v1$padj,1) < 0.05 &
                                        tidyr::replace_na(res_7v1$log2FoldChange,0) > 0]
elementMetadata(peaks_1v7) <- res_7v1[tidyr::replace_na(res_7v1$padj,1) < 0.05 &
                                        tidyr::replace_na(res_7v1$log2FoldChange,0) < 0]

elementMetadata(peaks_all) <- res_7v1

peaks_7v1_anno <- annotatePeak(peak = peaks_7v1,TxDb = txdb,
                               level = "gene",
                               assignGenomicAnnotation = T,
                               annoDb = "org.Mm.eg.db")
peaks_1v7_anno <- annotatePeak(peak = peaks_1v7,
                               TxDb = txdb,
                               level = "gene",
                               assignGenomicAnnotation = T,
                               annoDb = "org.Mm.eg.db")
peaks_all_anno <- annotatePeak(peak = peaks_all,
                               TxDb = txdb,
                               level = "gene",
                               assignGenomicAnnotation = T,
                               annoDb = "org.Mm.eg.db")

saveRDS(peaks_7v1_anno, file = "../data_output/diffbind_cutnrun_zic/peaks_7v1_up_annotation.Rds")
saveRDS(peaks_1v7_anno, file = "../data_output/diffbind_cutnrun_zic/peaks_7v1_down_annotation.Rds")
saveRDS(peaks_all_anno, file = "../data_output/diffbind_cutnrun_zic/peaks_7v1_all_annotation.Rds")


readr::write_tsv(as.data.frame(peaks_7v1_anno@anno), file = "../data_output/diffbind_cutnrun_zic/peaks_7v1_up_annotation.tsv")
readr::write_tsv(as.data.frame(peaks_1v7_anno@anno), file = "../data_output/diffbind_cutnrun_zic/peaks_7v1_down_annotation.tsv")
readr::write_tsv(as.data.frame(peaks_all_anno@anno), file = "../data_output/diffbind_cutnrun_zic/peaks_7v1_all_annotation.tsv")



