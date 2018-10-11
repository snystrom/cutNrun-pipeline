library(magrittr)
library(GenomicRanges)

thresholdPeakCall <- function(bigWig, signal_threshold = 5, merge_dist = 5, compute_score = T){
  # Input: bigWig path
  # signal_threshold = require this number of reads (or normalized reads) to be considered
  # merge_dist = max distance for two intervals to be merged. Typically set to 1/2 binsize of bigwig (?)
  # Output:
  # GRanges object w/ score column corresponding to sum of signal within peak
  bw <- rtracklayer::import.bw(bigWig)
  
  peaks <- bw %>% 
    data.frame %>% 
    dplyr::filter(score > signal_threshold) %>% 
    GRanges %>% 
    reduce(min.gapwidth = merge_dist) 
 
  if (compute_score == T) {
    peak_signal <- subsetByOverlaps(bw, peaks)
    mcols(peak_signal)$peak <- findOverlaps(bw, peaks) %>% subjectHits
    peaks$score <- peak_signal %>% 
      data.frame %>% 
      dplyr::group_by(peak) %>% 
      dplyr::summarise(score = sum(score)) %>% 
      .$score
    
    return(peaks)
  } else
  {
    return(peaks)
  }
}


main <- function() {
  args <- commandArgs(trailingOnly = T)
  bigwig_input <- args[1]
  peakFile <- args[2]
  
  peaks <- thresholdPeakCall(bigwig_input)
  
  rtracklayer::export.bed(peaks, peakFile)
}

main()
