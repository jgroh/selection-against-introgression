library(data.table)
library(tools)
library(magrittr)

# 1. Make 1kb bed file =====

# read files containing chromosome lengths
chrCM <- fread("cM_lengths_birchmanni10x.txt", col.names = c("chr", "cM")) # cM lengths
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
write.table(bed, file = "allChr_1kb.bed", quote = F, sep = "\t", col.names = F, row.names = F)


# 2. Make recombination bed file =====

# This reads all recombination files and constructs a bed file of rec intervals
# extending to the ends of chromosomes using the rec rates for the outermost intervals

recombBed <- rbindlist(
  lapply(list.files("LD_recMap/",full.names=T), function(x){
    rChrom <- fread(x, col.names = c("chr","start","end","mean_2Nerrr","V1","median_2Ner","V3"))
    rChrom <- rChrom[, c("chr", "start", "end", "median_2Ner")]
    
    if(min(rChrom$start) != 0){
      rChrom <- rbind(list(chr=rChrom[1,chr],
                           start=0,
                           end=min(rChrom$start),
                           median_2Ner=rChrom[start==min(start),median_2Ner]),
                      rChrom)
    }
    if(max(rChrom$end) < chrLen[chr==rChrom[1,chr],len]){
      rChrom <- rbind(rChrom,
                      list(chr=rChrom[1,chr],
                           start=max(rChrom$end),
                           end=chrLen[chr==rChrom[1,chr],len+1], # add +1 bc bed is zero-indexed
                           median_2Ner=rChrom[end==max(end),median_2Ner]))
    }
    rChrom
  }
  )
)

# Use regression of total Rho against Morgan lengths to get estimate of 2Ne
chrRho <- recombBed[, max(sum(rep(median_2Ner,end-start))), by = chr]
setnames(chrRho, "V1", "totalRho")

chrGenLen <- merge(chrCM, chrRho, by = "chr")
chrGenLen[, Morgan := cM/100]

z <- lm(totalRho ~ 0 + Morgan, data = chrGenLen) # force intercept through zero
with(chrGenLen, plot(totalRho ~ Morgan))

Ne2 <- z$coefficients[1]

# divide Rho values by 2Ne to get r
recombBed[, r := median_2Ner/Ne2][, median_2Ner := NULL]

# write out
write.table(bed, file = "allChrLDRecMap.bed", quote = F, sep = "\t", col.names = F, row.names = F)

