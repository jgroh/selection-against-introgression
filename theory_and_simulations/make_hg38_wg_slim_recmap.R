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

# ===== 1kb map, for neutral sims =====
# want evenly space physical positions, so that deleterious loci are evenly spaced on the physical map
xout1kb <- dcode[, seq(min(End), max(End), by = 1e3)]

# cm positions of evenly spaced bp positions
dcode_interp1kb <- dcode[, approx(x = End, y = cM, xout = seq(min(End), max(End), by = 1e3)), by = Chr]
setnames(dcode_interp1kb, c("chr", "pos", "cM"))

# calculate Morgan distances between adjacent SNPs
# Morgan_dist_slim will give recombination probabilities for Slim 
# (not bothering with map function as map lengths are so small )

dcode_interp1kb[, Morgan_dist := (cM-shift(cM))/100, by = chr]

# replace NA value at start of chroms
dcode_interp1kb[, Morgan_dist := c(Morgan_dist[2], Morgan_dist[2:nrow(.SD)]), by = chr]

# for slim:
# for subsequent chromosomes, set rec rate to 0.5 at first position 
dcode_interp1kb[, Morgan_dist_slim := Morgan_dist]
dcode_interp1kb[chr!="chr1", Morgan_dist_slim := c(0.5, Morgan_dist_slim[2:nrow(.SD)]), by = chr]

setkey(dcode_interp1kb, chr, pos)
dcode_interp1kb[, chr := gsub("chr", "", chr)]


if(chromosome == 1){
  fwrite(dcode_interp1kb[, .(Morgan_dist_slim)], file = 'hg38_wg_slim_recmap_1kb.txt', 
         col.names = F, row.names = F, quote = F, sep = '\t')
  fwrite(dcode_interp1kb[, .(chr, pos, Morgan = cM/100, Morgan_dist_slim)], file = 'hg38_wg_slim_recmap_1kb_verbose.txt.gz', 
         col.names = T, row.names = F, quote = F, sep = '\t')
}


dcode_interp1kb[, .SD[1], by = chr]
hg38_wg_slim_recmap_1kb.txt 

# ===== 50kb map, used for selection sims and same resolution analysis of real data =====

# want evenly space physical positions, so that deleterious loci are evenly spaced on the physical map
xout <- dcode[, seq(min(End), max(End), by = 50e3)]

# cm positions of evenly spaced bp positions
dcode_interp <- dcode[, approx(x = End, y = cM, xout = seq(min(End), max(End), by = 50e3)), by = Chr]
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

#plot(dcode_sub$r_bar_i)


# ----- map1
fwrite(dcode_sub[, .(chr, pos, cM, Morgan_dist_slim, r_bar_i)],
       file = paste0("hg38_slim_recmap_chrs/hg38_wg_slim_recmap_verbose_chr", chromosome, ".txt.gz"),
       quote=F, sep = "\t",
       col.names = F)

# ===== Map2: constant recombination rate along chromosomes =====

#fwrite(dcode[, .(len = max(cM)/100), by = Chr][, .(len)], 
#       file = "hg38_chr_lengths.txt",
#       quote=F, sep = "\t", 
#       col.names = F)

