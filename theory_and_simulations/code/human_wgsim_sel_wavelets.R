library(data.table)
#library(ggplot2)
library(cubature)
if(interactive()){
  source("/Users/jeff/workspace/gnomwav/R/correlation_decomp.R")
  source("/Users/jeff/workspace/gnomwav/R/multi_modwts.R")
  source("/Users/jeff/workspace/gnomwav/R/variance_decomp.R")
  frq <- fread('results/human_wgsim_sel1-1000/replicate0_frqs.txt', col.names = c("rep","gen","frq"))
  map <- fread("hg38_wg_slim_recmap_verbose.txt", col.names = c("chr", 'pos_bp', 'cM', 'Morgan_dist'))
  
} else{
  source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
  source("/Users/brogroh/gnomwav/R/multi_modwts.R")
  source("/Users/brogroh/gnomwav/R/variance_decomp.R")
  args <- commandArgs(trailingOnly = T)
  frq <- fread(args[1], col.names = c("rep", "gen", "frq"))
  map <- fread(args[2], col.names = c("chr", 'pos_bp', 'cM', 'Morgan_dist'))
  
}

# drop first position from slim results, to merge with SNP positions
frq <- frq[, .SD[-1,], by = .(rep, gen)]

# subtract gen to correspond to thry
frq[, gen := gen - 2]

# ----- format recombination map
map[, Morgan := cM/100][, cM := NULL]
setnames(map, "Morgan_dist", "rec")

# change rec rates at starts of chromosomes
map[chr != "chr1", rec := c(rec[2], rec[2:nrow(.SD)]), by = chr]

#ceiling(log2(map[, mean(Morgan_dist)]))
#ceiling(log2(map[, median(Morgan_dist)]))

# combine map and frqs
frq <- cbind(frq, map)


# for local testing
#frq <- frq[chr %in% c("chr21", "chr22")]

# ===== Interpolate ancestry and recombination for analyses in genetic units =====
#ceiling(log2(map[, mean(Morgan_dist)]))
#ceiling(log2(map[, median(Morgan_dist)]))

frq_interp <- frq[, approx(x = Morgan, y = frq, xout = seq(min(Morgan), max(Morgan), by = 2^-16)), by = .(rep, gen, chr)]
setnames(frq_interp, c("x","y"), c("Morgan","frq"))


rec_interp <- map[, approx(x = Morgan, y = pos_bp, xout = seq(min(Morgan), max(Morgan), by = 2^-16)), by = chr]
setnames(rec_interp, c("x","y"),c("Morgan","bp"))
rec_interp[, rec := 1/(bp-shift(bp))]
rec_interp[, rec := c(rec[2], rec[2:nrow(rec_interp)])]
#ggplot(rec_interp[seq(1, .N, by = 1000)], aes(x = Morgan, y = rec)) + geom_point() + facet_wrap(~chr)

# combine map on frqs on genetic map
frq_rec_interp <- merge(frq_interp,  rec_interp)


# ===== Total Variances: ancestry and recombination =====

# physical units
totalvarsP <- frq[,  .(totalvar.frq = var(frq), totalvar.rec = var(rec)), by = .(rep, gen)]
totalvarsP[, units := "physical"]
#ggplot(sel_totalvarP, aes(x = log10(gen), y = totalvar.frq)) + geom_point()

# genetic units
totalvarsG <- frq_rec_interp[, .(totalvar.frq = var(frq), totalvar.rec = var(rec)), by = .(rep, gen)]
totalvarsG[, units := "genetic"]

totalvars <- rbind(totalvarsP, totalvarsG)

# ===== Wavelet Variances =====

# ----- physical units -----
wv_frq_rec_P_chrs <- frq[, gnom_var_decomp(data=.SD, chromosome = 'chr', rm.boundary = F, signals = c("frq", "rec"), avg.over.chroms = F), by = .(rep, gen)]
wv_frq_rec_P_chrs[, units := "physical"]

wv_frq_rec_P_wg <- frq[, gnom_var_decomp(data=.SD, chromosome = 'chr', rm.boundary = F, signals = c("frq", "rec"), avg.over.chroms = T), by = .(rep, gen)]
wv_frq_rec_P_wg[, units := "physical"]

# ----- genetic units -----
wv_frq_rec_G_chrs <- frq_rec_interp[, gnom_var_decomp(data=.SD, chromosome = 'chr', rm.boundary = F, signals = c("frq", "rec"), avg.over.chroms = F), by = .(rep, gen)]
wv_frq_rec_G_chrs[, units := "genetic"]

wv_frq_rec_G_wg <- frq_rec_interp[, gnom_var_decomp(data=.SD, chromosome = 'chr', rm.boundary = F, signals = c("frq", "rec"), avg.over.chroms = T), by = .(rep, gen)]
wv_frq_rec_G_wg[, units := "genetic"]

wv_frq_rec_chrs_all <- rbind(wv_frq_rec_P_chrs, wv_frq_rec_G_chrs)
wv_frq_rec_wg_all <- rbind(wv_frq_rec_P_wg, wv_frq_rec_G_wg)

# ===== Wavelet Correlations ======

# look at total correlation
totalcorP <- frq[, .(cor = cor(frq, rec)), by = .(rep, gen)]
totalcorG <- frq_rec_interp[, .(cor = cor(frq, rec)), by = .(rep, gen)]

# ----- wavelet correlations  -----

# physical map
wavcorP <- frq[, cor_tbl(.SD, chromosome = 'chr', signals = c("frq", "rec"), rm.boundary = F), by = .(rep, gen)]
wavcorP[, units := "physical"]

# genetic map
wavcorG <- frq_rec_interp[, cor_tbl(.SD, chromosome = 'chr', signals = c("frq", "rec"), rm.boundary = F), by = .(rep, gen)]
wavcorG[, units := "genetic"]

wavcor_all <- rbind(wavcorP, wavcorG)

save(totalvars, wv_frq_rec_chrs_all, wv_frq_rec_wg_all, wavcor_all, file = gsub("_frqs.txt", "_wavelet_results.RData", args[1]))



