library(data.table)

# Make 1kb interval bed file for X. birchmanni chromosomes =====

# read file containing chromosome lengths
chrLen <- fread("xbir10x_chrlengths.txt", col.names = c("chr", "len")) # physical lengths

# midpoints of 1kb windows for each chromosome 
bed <- chrLen[, seq(500, max(len), by=5e2), by = .(chr,len)]
setnames(bed, "V1", "midpoint")
bed <- bed[bed$midpoint %% 1000 != 0] 

# make start and end of intervals. If last window of chrom, use chr length + 1 as end bc bed file is zero-indexed
bed[, c("start", "end") := .(midpoint - 500, 
                             ifelse(midpoint == max(midpoint), max(len)+1, midpoint + 500)), by = chr] 

bed[, "len" := NULL]
setcolorder(bed, c("chr", "start", "end", "midpoint"))

# write out
write.table(bed, file = "", quote = F, sep = "\t", col.names = F, row.names = F) 