#!/usr/bin/env Rscript

# Author: Spencer Nystrom (snystrom@email.unc.edu)
# Z-normalizes input bigwig file by chromosome arm

library(rtracklayer)

argv <- commandArgs(trailingOnly = TRUE)

if (length(argv) != 2){
   stop("usage: Rscript --vanilla zNorm.r <input.bw> <output.bw>") 
}

bw <- import.bw(argv[1])

bwList <- split(bw, bw@seqnames) # split by chromosome arm

zScoreChrList <- endoapply(bwList, function(x){
	# Calculate z-Score per chromosome arm:
	# z = (x-u)/sd
    chrScore <- x@elementMetadata$score
    zMean <- mean(chrScore)
    zSD <- sd(chrScore)
    x@elementMetadata$score <- (chrScore - zMean)/zSD
    return(x)
})


zScoreBw <- unlist(zScoreChrList) # recombine into Granges object

export.bw(con = argv[2], object = zScoreBw) # write bigwig


