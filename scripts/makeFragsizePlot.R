library(dplyr)
library(ggplot2)
library(GenomicRanges)





readNarrowPeak <- function(path){
  # read MACS2 narrowpeak file into data.frame, coerce to GRanges
  # Uses readr to enforce explicit data types.
  # Should not warn unless there's an issue.
  
  bedCols <- c("chr", "start", "end", "id", "score", "strand", "fold.change", "pValue", "qValue", "summitPos")
  bedDf <- readr::read_tsv(path,
                  col_names = bedCols, 
                  col_types = readr::cols(chr = readr::col_character(),
                                          start = readr::col_integer(),
                                          end = readr::col_integer(),
                                          id = readr::col_character(),
                                          score = readr::col_integer(),
                                          strand = readr::col_character(),
                                          fold.change = readr::col_double(),
                                          pValue = readr::col_double(),
                                          qValue = readr::col_double(),
                                          summitPos = readr::col_integer())
  )
  
  bed <- GenomicRanges::GRanges(bedDf)
  return(bed) 
}

readFragmentBed <- function(bed){
  # Imports a fragment bedfile, returns a GRanges object
  
  frags <- readr::read_tsv(bed, 
                           col_names = c("chr", "start", "end", "name", "score"),
                           col_types = readr::cols(chr = readr::col_character(),
                                                   start = readr::col_integer(),
                                                   end = readr::col_integer(),
                                                   name = readr::col_character(),
                                                   score = readr::col_integer())) %>% 
    GenomicRanges::GRanges()
  
  return(frags)
}

annotateFragmentsInPeaks <- function(frags, peaks){
  frags %>% 
    data.frame %>% 
    dplyr::mutate(inPeak = GRanges(.) %over% peaks) %>% 
    dplyr::mutate(inPeak = ifelse(inPeak, "Inside", "Outside"))
}


makeFragsizePlot <- function(frags_annotated, title = waiver()){
  frags_annotated %>% 
    ggplot(aes(width)) +
      geom_histogram(binwidth = 1, aes(fill = inPeak, color = inPeak)) +
      facet_wrap(~inPeak, scales = "free_y", ncol = 1) +
      theme_bw() +
      theme(legend.position = "none", 
            axis.text = element_text(size = 40/.pt),
            strip.text.x = element_text(size = 50/.pt),
            axis.title = element_text(size = 40/.pt)) +
      scale_x_continuous(breaks = seq(0, 1000, by = 50)) +
      labs(x = "Fragment Size",
           y = "Number of Fragments",
           title = title) 
}

writeFragsizePlot <- function(fragSizePlot, path, device = "png", defaults = T, ...){
  # saves a ggplot object w/ default width x height (10x10)
  # can set defaults = F to pass alternate params to ...
  # defaults to png image, can be changed with `device`
  
  if (defaults == T){
    ggsave(path, fragSizePlot, device = device, width = 10, height = 10, units = "in")
  } else {
    ggsave(path, fragSizePlot, device = device, ...)
  }
}

###
# MAIN
###

# By default the title of the plot is the name of the bed file
# TODO: use optArgs to parse cmdline flags instead of going by position
# ?TODO: Autodetect filetype from output extension?

# runs only when script is run by itself
# Thanks, Stack: https://stackoverflow.com/questions/2968220/is-there-an-r-equivalent-of-the-pythonic-if-name-main-main
if (sys.nframe() == 0){
  
  argv <- commandArgs(trailingOnly = TRUE)
  
  if (length(argv) != 3){
     stop("usage: Rscript --vanilla makeFragsizePlot.R <fragments.bed> <peaks.narrowPeak> <output_name>") 
  }
  
  
  fragment_bed_path <- argv[1]
  narrowPeak_path <- argv[2]
  output_name <- argv[3]
  
  peaks <- readNarrowPeak(narrowPeak_path)
  
  readFragmentBed(fragment_bed_path) %>% 
    annotateFragmentsInPeaks(., peaks) %>% 
    makeFragsizePlot(., title = fragment_bed_path) %>% 
    writeFragsizePlot(., output_name)
}

