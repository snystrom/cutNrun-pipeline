#!/usr/bin/env Rscript

# Author: Spencer Nystrom (snystrom@email.unc.edu)
# Z-normalizes input bigwig file by chromosome arm
# Usage: Rscript --vanilla zNorm.r <input.bw> <output.bw>

library(rtracklayer)

argv <- commandArgs(trailingOnly = TRUE)

if (length(argv) != 2){
   stop("usage: Rscript --vanilla zNorm.r <input.bw> <output.bw>") 
}


bw <- import.bw(argv[1])

bwList <- split(bw, bw@seqnames) # split by chromosome arm

# For filtering out specific regions you don't want. testing if using a blacklist file will overcome need for this.
#chrKeep <- c("chr2L", "chr2R", "chr3R", "chr3L", "chr4", "chrX")
#bwList <- bwList[chrKeep] # only keep chr2,3,4,X , no heterochromatin

write(paste("chr", "Mean", "Sd", sep = " "), stdout())
zScoreChrList <- lapply(bwList, function(x){
	# Calculate z-Score per chromosome arm:
	# z = (x-u)/sd
	# Mean and SD must be calculated as if from a frequency table, as deeptools will merge regions with identical scores
    chrScore <- x@elementMetadata$score

    nfreq <- sum(x@ranges@width)
    zMean <- (sum(x@ranges@width * x@elementMetadata$score)/nfreq)
    zSD <- sqrt((sum(x@ranges@width * (x@elementMetadata$score^2)) - nfreq * (zMean)^2)/(nfreq - 1))

    x@elementMetadata$score <- (chrScore - zMean)/zSD
    # Print ChrStats for report:
    chrName <- unique(x@seqnames)
    write(paste(chrName, zMean,zSD, sep = " "), stdout())
    #return(x)
    if (all(is.na(mcols(x)$score))) {
        # Chromosomes with no reads will have NaN score for whole chromosome
        return(NULL)
    } else {
        return(x)
    }
})

zScoreChrList_noNull <- Filter(Negate(is.null), zScoreChrList) 
zScoreChrList_noNull <- GRangesList(zScoreChrList_noNull)
zScoreBw <- unlist(zScoreChrList_noNull) # recombine into Granges object

export.bw(con = argv[2], object = zScoreBw) # write bigwig

