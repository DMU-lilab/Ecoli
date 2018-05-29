#!/usr/bin/R

# NOTE: This script is design for Ecoli genome
# The RLE funciton (RegionAsBed) has a broad usage

RegionAsBed <- function(marker, cf.length=5, tolerance.length=NULL,to_up_len = 2,chrom) {
	# cf.length: At least 5 values(p value <= cutoff)
	# tolerance.lenght: tolerance for False within the hotsport region
	# to_up_len: the minimum value for continuous True within the region
	#library(zoo)
	r <- rle(marker)
	if(is.null(tolerance.length)){
		end <- with(r,cumsum(lengths)[values & lengths>=cf.length])
		start <- end - with(r,lengths[values & lengths>=cf.length]) + 1
		if(length(start) == 0 || length(end) == 0){
			df <- data.frame(chrom=character(), start=integer(), end=integer())
		} else {
			df <- data.frame(chrom, start, end)
		}
	} else {
		end<-cumsum(r$lengths)
		start<- end - r$lengths + 1
		df.tmp <- data.frame(start = start  , end = end , value = r$values)
		df.tmp$len <- df.tmp$end - df.tmp$start + 1 
		df.tolerance <- df.tmp[!df.tmp$value & df.tmp$len <= tolerance.length,] 
		tolerance_rowname <- as.numeric(rownames(df.tolerance))
		tolerance_up <- tolerance_rowname - 1 
		tolerance_down <- tolerance_rowname + 1 
		tolerance_mark <- tolerance_rowname[df.tmp$len[tolerance_up] >= to_up_len & df.tmp$len[tolerance_down] >= to_up_len]
		dt.to <- data.table(df.tmp[tolerance_mark,])
		if (nrow(dt.to)!=0) {
			dt.to$id <- 1:nrow(dt.to)
			dt.to.index <- dt.to[,.(index = start:end), by = id] 		
			marker[dt.to.index$index] <- TRUE
		}
		
		r <- rle(marker)
        	end <- with(r,cumsum(lengths)[values & lengths>=cf.length])
	        start <- end - with(r,lengths[values & lengths>=cf.length]) + 1
        	if(length(start) == 0 || length(end) == 0){
            		df <- data.frame(chrom=character(), start=integer(), end=integer())
        	} else {
            		df <- data.frame(chrom, start, end)
        	}
    	}	
    return(df)
}

library(data.table) 

args <- commandArgs(T)

if (length(args) < 2) {
    message("\nUsage: Rscript GetHotspot.R <input.txt> <output.txt>\n")
    q()
}

input <- args[1] # eg. input.txt
output <- args[2] # eg. output.txt

# read the input file
a <- read.table(input,header=T)

# calculate the Pvalue
a$Pvalue <- 1 - exp(-(1/mean(a[,2])) * a[,2])
a$marker <- a$Pvalue < 0.05

# excute the RLE algorithm
m <- a$marker
hotsport <- RegionAsBed(m,5,1,2,"Ecoli")
hotsport$start <- a$Pos[hotsport$start]
hotsport$end <- a$Pos[hotsport$end]

# write out the result
write.table(hotsport, file=output, row.names=F, quote=F, sep="\t")
