library(data.table)
options(scipen=999)

# Make bed file of recombination per bp recombination rates for X. birchmanni genome =====
chrLen <- fread("xbir10x_chrlengths.txt", col.names = c("chr", "len")) # physical lengths
chrCM <- fread("cM_lengths_birchmanni10x.txt", col.names = c("chr", "cM")) # cM lengths

# This combines all recombination files,
# extending intervals to the ends of chromosomes using the rec rates for the outermost intervals

recombBed <- rbindlist(
  lapply(list.files("LD_recMap/",full.names=T), function(x){
    rChrom <- fread(x, col.names = c("chr","start","end","mean_2Ner","V1","median_2Ner","V3"))
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

# what portion of the genome is above the threshold? 1.6%
#sum(recombBed[median_2Ner >= 0.005, end-start])/sum(recombBed[, end-start])
#recombBed[median_2Ner >= 0.005, median_2Ner := 0.005]
#hist(log(recombBed[,rep(median_2Ner, end-start)]))


# Use regression of total Rho against Morgan lengths to get estimate of 2Ne
chrRho <- recombBed[, max(sum(rep(median_2Ner,end-start))), by = chr]
setnames(chrRho, "V1", "totalRho")

chrGenLen <- merge(chrCM, chrRho, by = "chr")
chrGenLen[, Morgan := cM/100]

z <- lm(totalRho ~ 0 + Morgan, data = chrGenLen) # force intercept through zero
with(chrGenLen, plot(totalRho ~ Morgan))

summary(z)
Ne2 <- z$coefficients[1]

# divide Rho values by 2Ne to get r
recombBed[, r := median_2Ner/Ne2][, median_2Ner := NULL]

# write out
# since this will be the map file (for bedtools), needs 5 columns where 5th column is score
recombBed[, name := paste0("map-",seq_len(.N))]
setcolorder(recombBed, c("chr","start","end","name","r"))
write.table(recombBed, file = "", quote = F, sep = "\t", col.names = F, row.names = F)
