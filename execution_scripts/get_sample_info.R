library(tidyverse)

reports_tsv <- map(list.files("data_output/ena_reports/full_reports", full.names = TRUE), read_tsv)
chip_chromatin_ids <- readLines("execution_scripts/chip_chromatin_ids.txt")

reports_tsv <- do.call(rbind, map(reports_tsv, function(x) {
  x <- x[x$run_accession %in% chip_chromatin_ids,]
  return(x)
}))

reports_tsv %>%
  filter(library_layout == "PAIRED") %>%
  separate(fastq_ftp, sep = ";", into = c("r1","r2")) %>%
  select(r1,r2) %>%
  unlist() %>%
  writeLines(con = "execution_scripts/chip_chromatin_pe_ftp.txt")

reports_tsv %>%
  filter(library_layout == "PAIRED") %>%
  pull(run_accession)

write_tsv(reports_tsv, "execution_scripts/chip_chromatin_sample_info.tsv")
