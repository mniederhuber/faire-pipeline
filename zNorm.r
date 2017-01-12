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

# For filtering out specific regions you don't want. testing if using a blacklist file will overcome need for this.
#chrKeep <- c("chr2L", "chr2R", "chr3R", "chr3L", "chr4", "chrX")
#bwList <- bwList[chrKeep] # only keep chr2,3,4,X , no heterochromatin

print(paste("chr", "Mean", "Sd", sep = " "))
zScoreChrList <- endoapply(bwList, function(x){
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
    return(x)
})


zScoreBw <- unlist(zScoreChrList) # recombine into Granges object

zScoreCSV <- data.frame(zScoreBw)
write.csv(zScoreCSV, "./zScoreCSV.csv", row.names = F)

export.bw(con = argv[2], object = zScoreBw) # write bigwig


