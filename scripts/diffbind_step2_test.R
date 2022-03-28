library(DESeq2)

# Load deseq obj
deseq_obj <- readRDS("../data_output/diffbind_cutnrun_zic/deseq_obj.Rds")

# Prefiltering
keep <- rowSums(counts(deseq_obj)) >= 20
deseq_obj <- deseq_obj[keep,]

# Build model
deseq_obj_contrasts <- DESeq(deseq_obj)

# Test
res_3v1 <- results(deseq_obj_contrasts, contrast = c("time","D3Z","D1Z"),
                   alpha = 0.05)
res_5v3 <- results(deseq_obj_contrasts, contrast = c("time","D5Z","D3Z"),
                   alpha = 0.05)
res_7v5 <- results(deseq_obj_contrasts, contrast = c("time","D7Z","D5Z"),
                   alpha = 0.05)
res_7v1 <- results(deseq_obj_contrasts, contrast = c("time","D7Z","D1Z"),
                   alpha = 0.05)

# Test with LRT analogous to anova 
deseq_obj_lrt <- DESeq(deseq_obj,
                       test="LRT",
                       reduced=~1)

res_all <- results(deseq_obj_lrt, alpha = 0.05)


# Shrinkage
res_3v1_ashr <- lfcShrink(deseq_obj_contrasts,
                          res  = res_3v1,
                          type="ashr")
res_5v3_ashr <- lfcShrink(deseq_obj_contrasts,
                          res  = res_5v3,
                          type="ashr")
res_7v5_ashr <- lfcShrink(deseq_obj_contrasts,
                          res  = res_7v5,
                          type="ashr")
res_7v1_ashr <- lfcShrink(deseq_obj_contrasts,
                          res  = res_7v1,
                          type="ashr")
res_all_ashr <- lfcShrink(deseq_obj_contrasts,
                          res  = res_all,
                          type="ashr")

# Get diff ranges

get_ranges_up <- function(res, deseq, padj = 0.05,lfc = 0){
  granges(rowRanges(deseq)[tidyr::replace_na(res$padj,1) < padj & tidyr::replace_na(res$log2FoldChange,0) > lfc])
}
get_ranges_down <- function(res, deseq, padj = 0.05,lfc = 0){
  granges(rowRanges(deseq)[tidyr::replace_na(res$padj,1) < padj & tidyr::replace_na(res$log2FoldChange,0) < lfc])
}


res_3v1_up <- get_ranges_up(res_3v1_ashr, deseq_obj_contrasts)
res_3v1_down <-  get_ranges_down(res_3v1_ashr, deseq_obj_contrasts)

res_5v3_up <- get_ranges_up(res_5v3_ashr, deseq_obj_contrasts)
res_5v3_down <-  get_ranges_down(res_5v3_ashr, deseq_obj_contrasts)

res_7v5_up <- get_ranges_up(res_7v5_ashr, deseq_obj_contrasts)
res_7v5_down <-  get_ranges_down(res_7v5_ashr, deseq_obj_contrasts)

res_7v1_up <- get_ranges_up(res_7v1_ashr, deseq_obj_contrasts)
res_7v1_down <-  get_ranges_down(res_7v1_ashr, deseq_obj_contrasts)

res_all_any <- granges(rowRanges(deseq_obj_lrt)[tidyr::replace_na(res_all_ashr$padj,1) < 0.05 ])

res_7v1_up_lfc2 <- get_ranges_up(res_7v1_ashr, deseq_obj_contrasts, lfc = 2)
res_7v1_down_lfc2 <-  get_ranges_down(res_7v1_ashr, deseq_obj_contrasts, lfc = -2)

# Write results objects
saveRDS(res_3v1_ashr, "../data_output/diffbind_cutnrun_zic/res_3v1_ashr.Rds")
saveRDS(res_5v3_ashr, "../data_output/diffbind_cutnrun_zic/res_5v3_ashr.Rds")
saveRDS(res_7v5_ashr, "../data_output/diffbind_cutnrun_zic/res_7v5_ashr.Rds")
saveRDS(res_7v1_ashr, "../data_output/diffbind_cutnrun_zic/res_7v1_ashr.Rds")
saveRDS(res_all_ashr, "../data_output/diffbind_cutnrun_zic/res_all_ashr.Rds")

# Write results tables
deseq_obj_contrasts_ranges <- as.data.frame(granges(rowRanges(deseq_obj_contrasts)))[,1:3]

readr::write_tsv(cbind(deseq_obj_contrasts_ranges,as.data.frame(res_3v1_ashr)),
                 "../data_output/diffbind_cutnrun_zic/res_3v1_ashr.tsv")
readr::write_tsv(cbind(deseq_obj_contrasts_ranges,as.data.frame(res_5v3_ashr)),
                 "../data_output/diffbind_cutnrun_zic/res_5v3_ashr.tsv")
readr::write_tsv(cbind(deseq_obj_contrasts_ranges,as.data.frame(res_7v5_ashr)),
                 "../data_output/diffbind_cutnrun_zic/res_7v5_ashr.tsv")
readr::write_tsv(cbind(deseq_obj_contrasts_ranges,as.data.frame(res_7v1_ashr)),
                 "../data_output/diffbind_cutnrun_zic/res_7v1_ashr.tsv")
readr::write_tsv(cbind(deseq_obj_contrasts_ranges,as.data.frame(res_all_ashr)),
                 "../data_output/diffbind_cutnrun_zic/res_all_ashr.tsv")

# Write deseq fitted objects
saveRDS(deseq_obj_contrasts, "../data_output/diffbind_cutnrun_zic/deseq_obj_contrasts.Rds")
saveRDS(deseq_obj_lrt, "../data_output/diffbind_cutnrun_zic/deseq_obj_lrt.Rds")

# DiffBind bedfiles

rtracklayer::export.bed(res_3v1_up,"../data_output/diffbind_cutnrun_zic/res_3v1_up.bed")
rtracklayer::export.bed(res_3v1_down,"../data_output/diffbind_cutnrun_zic/res_3v1_down.bed")


rtracklayer::export.bed(res_5v3_up,"../data_output/diffbind_cutnrun_zic/res_5v3_up.bed")
rtracklayer::export.bed(res_5v3_down,"../data_output/diffbind_cutnrun_zic/res_5v3_down.bed")


rtracklayer::export.bed(res_7v5_up,"../data_output/diffbind_cutnrun_zic/res_7v5_up.bed")
rtracklayer::export.bed(res_7v5_down,"../data_output/diffbind_cutnrun_zic/res_7v5_down.bed")


rtracklayer::export.bed(res_7v1_up,"../data_output/diffbind_cutnrun_zic/res_7v1_up.bed")
rtracklayer::export.bed(res_7v1_down,"../data_output/diffbind_cutnrun_zic/res_7v1_down.bed")

rtracklayer::export.bed(res_all_any,"../data_output/diffbind_cutnrun_zic/res_all_any.bed")

rtracklayer::export.bed(res_7v1_up_lfc2,"../data_output/diffbind_cutnrun_zic/res_7v1_up_lfc2.bed")
rtracklayer::export.bed(res_7v1_down_lfc2,"../data_output/diffbind_cutnrun_zic/res_7v1_down_lfc2.bed")



