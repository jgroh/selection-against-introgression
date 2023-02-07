library(data.table)

chromosome <- commandArgs(trailingOnly = T)[1]
dcode <- fread("aau1043_datas3")

dcode <- dcode[Chr != "chrX"]
dcode[, Chr := factor(Chr, levels = paste0('chr', 1:22))]
setkey(dcode, Chr)

# cumulative cM
#ggplot(dcode, aes(x = Begin, y = cM)) + geom_point() 

# rates
#ggplot(dcode, aes(x = Begin, y = cMperMb)) + geom_point()

# ===== Map1: realistic human rec map at evenly spaced physical positions =====

# want evenly space physical positions, so that deleterious loci are evenly spaced on the physical map
xout <- dcode[, seq(min(End), max(End), by = 1000)]

# cm positions of evenly spaced bp positions
dcode_interp <- dcode[, approx(x = End, y = cM, xout = seq(min(End), max(End), by = 1e3)), by = Chr]
setnames(dcode_interp, c("chr", "pos", "cM"))

# calculate Morgan distances between adjacent SNPs
# will have two distances: Morgan_dist used for calculation of r bar, 
# and Morgan_dist_slim which will describe recombination probabilities for Slim

dcode_interp[, Morgan_dist := (cM-shift(cM))/100, by = chr]
# not bothering with map function as map lengths are so small 

# replace NA value at start of chroms
dcode_interp[, Morgan_dist := c(Morgan_dist[2], Morgan_dist[2:nrow(.SD)]), by = chr]

# for slim:
# for subsequent chromosomes, set rec rate to 0.5 at first position 
dcode_interp[, Morgan_dist_slim := Morgan_dist]
dcode_interp[chr!="chr1", Morgan_dist_slim := c(0.5, Morgan_dist_slim[2:nrow(.SD)]), by = chr]

dcode_interp[, chr := factor(chr, levels = paste0('chr', 1:22))]
setkey(dcode_interp, chr)


# ----- calculate per locus contribution to r bar
# calculate one chromosome at a time to parallelize
dcode_sub <- dcode_interp[chr == paste0("chr", chromosome)]

L_all <- nrow(dcode_interp)
L_diff_chrom <- L_all - nrow(dcode_sub)


for(i in 1:nrow(dcode_sub)){

  Morgan_dist_same_chrom <- dcode_sub[i, cM/100] - dcode_sub[, cM/100]
  rijs_same_chrom <- (1/2) * ( 1 - exp(-2*abs(Morgan_dist_same_chrom)))
  ri_same_chrom_component <- (1/L_all^2)*sum(rijs_same_chrom)
  dcode_sub[i, r_bar_i_pt1 := ri_same_chrom_component][]
}
dcode_sub[, r_bar_i := r_bar_i_pt1 + (1/L_all^2)*(0.5)*L_diff_chrom]

plot(dcode_sub$r_bar_i)


# ----- map1
fwrite(dcode_sub[, .(chr, pos, cM, Morgan_dist_slim, r_bar_i)],
       file = paste0("hg38_wg_slim_recmap_verbose_chr", chromosome, ".txt.gz"),
       quote=F, sep = "\t",
       col.names = F)

# ===== Map2: constant recombination rate along chromosomes =====

fwrite(dcode[, .(len = max(cM)/100), by = Chr][, .(len)], 
       file = "hg38_chr_lengths.txt",
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
