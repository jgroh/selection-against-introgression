library(tools)
library(magrittr)

# Make bed files for gene density

# given a genome file that contains scaffold lengths and chromosome lengths, 
# make bed file for interpolation intervals

# get max snp position for each chromosome for interpolation intervals
par1 <- fread(paste0("ancestry-probs-par1_allchrs_ACUA_historical_2018.tsv"))
a <- as.data.table(str_split_fixed(colnames(par1)[-1], ":", 2))
setnames(a, c("scaff", "len"))
scaffLengths <- a[, as.numeric(max(len)), by=scaff]
setnames(scaffLengths, "V1", "len")

# make bed 
bed <- scaffLengths[, seq(0, max(len), by=1e3), by = .(scaff,len)]
setnames(bed, "V1", "midpoint")

bed[, c("start", "end") := .(ifelse(midpoint == 0, 0, midpoint - 500), 
                             ifelse(midpoint == max(midpoint), max(len)+1, midpoint + 500)), by = scaff]
bed[, c("len","midpoint") := NULL]

# write out
write.table(bed, file = "allChr_1kb.bed", quote = F, sep = "\t", col.names = F, row.names = F)
