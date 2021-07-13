library(data.table)
options(scipen=999)

# read file containing chromosome lengths
args <- commandArgs(trailingOnly=T)
chrLenFile <- args[1]

chrLen <- fread(chrLenFile, col.names = c("chr", "len")) # physical lengths

# midpoints of 5kb windows for each chromosome 
bed <- chrLen[, seq(2500, max(len), by=2500), by = .(chr,len)]
setnames(bed, "V1", "midpoint")
bed <- bed[bed$midpoint %% 5000 != 0] 

# make start and end of intervals. If last window of chrom, use chr length + 1 as end bc bed file is zero-indexed
bed[, c("start", "end") := .(midpoint - 2500, 
                             ifelse(midpoint == max(midpoint), max(len)+1, midpoint + 2500)), by = chr] 

bed[, "len" := NULL][, "midpoint" := NULL]
setcolorder(bed, c("chr", "start", "end"))

# write out
write.table(bed, file = "", quote = F, sep = "\t", col.names = F, row.names = F) 
