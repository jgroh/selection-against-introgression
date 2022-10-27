library(data.table)

dcode <- fread("aau1043_datas3")
dcode <- dcode[Chr== "chr1"]

# cumulative cM
#ggplot(dcode, aes(x = Begin, y = cM)) + geom_point() 

# rates
#ggplot(dcode, aes(x = Begin, y = cMperMb)) + geom_point()

# want evenly space physical positions, so that deleterious loci are evenly spaced on the physical map
xout <- dcode[, seq(min(End), max(End), by = 1000)]

# cm positions of evenly spaced bp positions
dcode_interp <- dcode[, approx(x = End, y = cM, xout = xout)]
setnames(dcode_interp, c("pos", "cM"))
dcode_interp[, Morgan_dist := (cM-shift(cM))/100]
dcode_interp[, Morgan_dist := c(Morgan_dist[2], Morgan_dist[2:nrow(dcode_interp)])]

#ggplot(dcode_interp, aes(x=pos, y = Morgan_dist)) + geom_point()

fwrite(dcode_interp[, Morgan_dist],
       file = "hg38_chr1_slim_recmap.txt",
       quote=F, sep = "\t",
       col.names = F)

fwrite(dcode_interp[, .(pos, Morgan_dist)],
       file = "hg38_chr1_slim_recmap_verbose.txt",
       quote=F, sep = "\t",
       col.names = F)












# #library(ggplot2)
# hg38 <- fread("recomb-hg38/genetic_map_GRCh38_merged.tab")
# chr1 <- hg38[chrom == "chr1"]
# 
# #chr1[, myrate := (pos_cm - shift(pos_cm))/(pos-shift(pos))]
# #chr1[, boxplot(myrate)]
# #which(chr1$recomb_rate > 100)
# #chr1[11977:11980]
# #ggplot(chr1[11000:12500], aes(x = pos, y = pos_cm)) + geom_point() + geom_line(aes(x = 13376798))
# 
# # clearly a major outlier - is this believable for recombination hotspots?
# # try dropping this SNP and recalculating?
# chr1 <- chr1[recomb_rate < 200]
# chr1[, myrate := (pos_cm/100 - shift(pos_cm)/100)]
# #ggplot(chr1[11000:12500], aes(x = pos, y = pos_cm)) + geom_point() + geom_line(aes(x = 13376798))
# 
# chr1[, length(unique(pos))] # need to know this number to input number of loci in slim
# out <- chr1[!is.na(myrate), .(pos, myrate)]
# out[, max(pos)-min(pos)]
# 
# fwrite(chr1[!is.na(myrate), .(myrate)], 
#        file = "hg38_chr1_slim_recmap.txt", 
#        quote=F, sep = "\t", 
#        col.names = F)
# 
# fwrite(chr1[!is.na(myrate), .(pos, pos_cm, myrate)], 
#        file = "hg38_chr1_slim_recmap_verbose.txt", 
#        quote=F, sep = "\t", 
#        col.names = F)
