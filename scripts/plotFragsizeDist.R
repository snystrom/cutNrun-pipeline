library(dplyr)
library(GenomicRanges)
library(ggplot2)

args <- commandArgs(trailingOnly = T)

bedFile <- args[1]

bed <- readr::read_tsv(bedFile, col_names = c("seqnames",
                                              "start",
                                              "end",
                                              "id",
                                              "score"))

bed %>% 
  GRanges %>% 
  data.frame %>% 
  ggplot(aes(width)) +
    geom_histogram(binwidth = 1) +
    ggtitle(bedFile) +
    theme_bw()

outFile <- glue::glue(basename(bedFile), "_fragSize.png")
ggsave(outFile, width = 15, height = 7.5)