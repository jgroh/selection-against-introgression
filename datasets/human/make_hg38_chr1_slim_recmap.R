library(data.table)
library(ggplot2)
hg38 <- fread("recomb-hg38/genetic_map_GRCh38_merged.tab")
chr1 <- hg38[chrom == "chr1"]

#chr1[, myrate := (pos_cm - shift(pos_cm))/(pos-shift(pos))]
#chr1[, boxplot(myrate)]
#which(chr1$recomb_rate > 100)
#chr1[11977:11980]
#ggplot(chr1[11000:12500], aes(x = pos, y = pos_cm)) + geom_point() + geom_line(aes(x = 13376798))

# clearly a major outlier - is this believable for recombination hotspots?
# try dropping this SNP and recalculating?
chr1 <- chr1[recomb_rate < 200]
chr1[, myrate := (pos_cm/100 - shift(pos_cm)/100)]
#ggplot(chr1[11000:12500], aes(x = pos, y = pos_cm)) + geom_point() + geom_line(aes(x = 13376798))

chr1[, length(unique(pos))] # need to know this number to input number of loci in slim
out <- chr1[!is.na(myrate), .(pos, myrate)]
out[, max(pos)-min(pos)]

fwrite(chr1[!is.na(myrate), .(myrate)], 
       file = "hg38_chr1_slim_recmap.txt", 
       quote=F, sep = "\t", 
       col.names = F)

fwrite(chr1[!is.na(myrate), .(pos, pos_cm, myrate)], 
       file = "hg38_chr1_slim_recmap_verbose.txt", 
       quote=F, sep = "\t", 
       col.names = F)
