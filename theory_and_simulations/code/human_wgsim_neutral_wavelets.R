library(data.table)
library(gnomwav)

if(interactive()){
  frq <- fread('results/human_wgsim_neutral_equilibrium/replicate0_frqs.txt', col.names = c("rep","gen","frq"))
  map <- fread("hg38_wg_slim_recmap_1kb_verbose.txt.gz",  col.names = c("chr", 'pos_bp','Morgan', 'r'))
  
} else{
  args <- commandArgs(trailingOnly = T)
  frq <- fread(args[1], col.names = c("rep", "gen", "frq"))
  map <- fread(args[2], col.names = c("chr", 'pos_bp', 'Morgan', 'r'))
}

# drop first position from slim results, to merge with SNP positions
frq <- frq[, .SD[-1,], by = .(rep, gen)]

# subtract gen to correspond to thry
frq[, gen := gen - 2]

# ----- format recombination map
# change rec rates at starts of chromosomes (they were 0.5 for slim simulations)
map[, r := c(r[2], r[2:nrow(.SD)]), by = chr]

#ceiling(log2(map[, mean(Morgan_dist)]))
#ceiling(log2(map[, median(Morgan_dist)]))

# combine map and frqs
frq <- cbind(frq, map)


# for local testing
#frq <- frq[chr %in% c("chr21", "chr22")]


# ===== Interpolate ancestry and recombination for analyses in genetic units =====
#ceiling(log2(map[, mean(Morgan_dist)]))
#ceiling(log2(map[, median(Morgan_dist)]))

frq_interp <- frq[, approx(x = Morgan, y = frq, xout = seq(min(Morgan), max(Morgan), by = 2^-15)), by = .(rep, gen, chr)]
setnames(frq_interp, c("x","y"), c("Morgan","frq"))


rec_interp <- map[, approx(x = Morgan, y = pos_bp, xout = seq(min(Morgan), max(Morgan), by = 2^-15)), by = chr]
setnames(rec_interp, c("x","y"),c("Morgan","bp"))
rec_interp[, r := 2^-15/(bp-shift(bp))]
rec_interp[, r := c(r[2], r[2:nrow(rec_interp)])]
#ggplot(rec_interp[seq(1, .N, by = 1000)], aes(x = Morgan, y = r)) + geom_point() + facet_wrap(~chr)

# combine map on frqs on genetic map
frq_rec_interp <- merge(frq_interp, rec_interp)


# ===== Total Variances: ancestry and recombination =====

# physical units
totalvarsP <- frq[,  .(totalvar.frq = var(frq), totalvar.r = var(r)), by = .(rep, gen)]
totalvarsP[, units := "physical"]

# genetic units
totalvarsG <- frq_rec_interp[, .(totalvar.frq = var(frq), totalvar.r = var(r)), by = .(rep, gen)]
totalvarsG[, units := "genetic"]

totalvars <- rbind(totalvarsP, totalvarsG)

# ===== Wavelet Variances =====

# # ----- physical units -----
# wv_frq_rec_P_chrs <- frq[, gnom_var_decomp(data=.SD, chromosome = 'chr', rm.boundary = T, signals = c("frq", "r"), avg.over.chroms = F), by = .(rep, gen)]
# wv_frq_rec_P_chrs[, units := "physical"]
# 
# wv_frq_rec_P_wg <- frq[, gnom_var_decomp(data=.SD, chromosome = 'chr', rm.boundary = T, signals = c("frq", "r"), avg.over.chroms = T), by = .(rep, gen)]
# wv_frq_rec_P_wg[, units := "physical"]

# ----- genetic units -----
#wv_frq_rec_G_chrs <- frq_rec_interp[, gnom_var_decomp(data=.SD, chromosome = 'chr', rm.boundary = T, signals = c("frq", "r"), avg.over.chroms = F), by = .(rep, gen)]
#wv_frq_rec_G_chrs[, units := "genetic"]

wv_frq_rec_G_wg <- frq_rec_interp[, gnom_var_decomp(data=.SD, chromosome = 'chr', rm.boundary = T, signals = c("frq", "r"), avg.over.chroms = T), by = .(rep, gen)]
wv_frq_rec_G_wg[, units := "genetic"]

#wv_frq_rec_chrs_all <- rbind(wv_frq_rec_P_chrs, wv_frq_rec_G_chrs)
#wv_frq_rec_wg_all <- rbind(wv_frq_rec_P_wg, wv_frq_rec_G_wg)
wv_frq_rec_wg_all <- wv_frq_rec_G_wg

save(totalvars, wv_frq_rec_wg_all, file = gsub("_frqs.txt.gz", "_wavelet_results.RData", args[1]))



