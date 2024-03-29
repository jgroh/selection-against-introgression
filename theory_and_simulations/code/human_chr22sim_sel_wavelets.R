library(data.table)
#library(ggplot2)
library(cubature)
library(wCorr)

if(interactive()){
  setwd("/Users/Jeff/workspace/selection-against-introgression/theory_and_simulations/")
  source("/Users/jeff/workspace/gnomwav/R/correlation_decomp.R")
  source("/Users/jeff/workspace/gnomwav/R/multi_modwts.R")
  source("/Users/jeff/workspace/gnomwav/R/variance_decomp.R")
  frq <- fread('results/human_chr22sim_sel900-1000_S1/replicate0_frqs.txt.gz', col.names = c("rep","gen","frq"))
  map <- fread("hg38_chr22_slim_recmap_verbose.txt.gz", col.names = c("chr", 'pos_bp', 'Morgan', 'r', "rbar_i"))
  
} else{
  source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
  source("/Users/brogroh/gnomwav/R/multi_modwts.R")
  source("/Users/brogroh/gnomwav/R/variance_decomp.R")
  args <- commandArgs(trailingOnly = T)
  frq <- fread(args[1], col.names = c("rep", "gen", "frq"))
  map <- fread(args[2], col.names = c("chr", 'pos_bp', 'Morgan', 'r', "rbar_i"))
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
plot(frq[gen==1, rbar_i])

# ===== Interpolate ancestry and recombination for analyses in genetic units =====
#ceiling(log2(map[, mean(Morgan_dist)]))
#ceiling(log2(map[, median(Morgan_dist)]))

frq_interp <- frq[, approx(x = Morgan, y = frq, xout = seq(min(Morgan), max(Morgan), by = 2^-16)), by = .(rep, gen, chr)]
setnames(frq_interp, c("x","y"), c("Morgan","frq"))


rec_interp <- map[, approx(x = Morgan, y = pos_bp, xout = seq(min(Morgan), max(Morgan), by = 2^-16)), by = chr]
setnames(rec_interp, c("x","y"),c("Morgan","bp"))
rec_interp[, r := 1/(bp-shift(bp))] # exact units don't matter for correlation purposes but it's Morgan/bp times a constant
rec_interp[, r := c(r[2], r[2:nrow(rec_interp)])]
#ggplot(rec_interp[seq(1, .N, by = 1000)], aes(x = Morgan, y = r)) + geom_point() + facet_wrap(~chr)

rbari_interp <- map[, approx(x = Morgan, y = rbar_i, xout = seq(min(Morgan), max(Morgan), by = 2^-16)), by = chr]
setnames(rbari_interp, c("x","y"),c("Morgan","rbar_i"))

# combine map on frqs on genetic map
frq_rec_interp <- merge(merge(frq_interp,  rec_interp), rbari_interp)

#ggplot(frq_rec_interp[gen == 1,], aes(x = Morgan, y = rbar_i)) + geom_point() + facet_wrap(~chr)

# ===== Mean ancestry proportion =====

meanfrq <- frq[, .(meanfrq = mean(frq)), by = .(rep, gen)]

# ===== Total Variances: ancestry and recombination =====

# physical units
totalvarsP <- frq[,  .(totalvar.frq = var(frq), 
                       totalvar.r = var(r),
                       totalvar.rbar_i = var(rbar_i)), by = .(rep, gen)]
totalvarsP[, units := "physical"]
#ggplot(sel_totalvarP, aes(x = log10(gen), y = totalvar.frq)) + geom_point()

# genetic units
totalvarsG <- frq_rec_interp[, .(totalvar.frq = var(frq), 
                                 totalvar.r = var(r),
                                 totalvar.rbar_i = var(rbar_i)), by = .(rep, gen)]
totalvarsG[, units := "genetic"]

totalvars <- rbind(totalvarsP, totalvarsG)


# ====== Total correlations =====

totalcorP <- frq[, .(cor_frq_r = cor(r, frq),
                     cor_frq_rbar_i = cor(rbar_i, frq)), 
                 by = .(rep, gen)]
totalcorP[, units := 'physical']

totalcorG <- frq_rec_interp[, .(cor_frq_r = cor(r, frq),
                                cor_frq_rbar_i = cor(rbar_i, frq)), 
                            by = .(rep, gen)]
totalcorG[, units := 'genetic']

totalcors <- rbind(totalcorG, totalcorP)


# ===== Wavelet Variances =====

# ----- physical units -----

wv_frq_rec_P_wg <- frq[, gnom_var_decomp(data=.SD, chromosome = NA, rm.boundary = F, signals = c("frq", "r", "rbar_i"), avg.over.chroms = T), by = .(rep, gen)]
wv_frq_rec_P_wg[, units := "physical"]

# ----- genetic units -----

wv_frq_rec_G_wg <- frq_rec_interp[, gnom_var_decomp(data=.SD, chromosome = NA, rm.boundary = F, signals = c("frq", "r", "rbar_i"), avg.over.chroms = T), by = .(rep, gen)]
wv_frq_rec_G_wg[, units := "genetic"]

wv_frq_rec_wg_all <- rbind(wv_frq_rec_P_wg, wv_frq_rec_G_wg)


# ===== Wavelet Correlations ======


# ----- by chromosome and averaging across chromosomes  -----

# ----- physical map

wavcorP_wg_r <- frq[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("frq", "r"), rm.boundary = F), by = .(rep, gen)]
wavcorP_wg_r[, vars := "frq_r"]

wavcorP_wg_rbar_i <- frq[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("frq", "rbar_i"), rm.boundary = F), by = .(rep, gen)]
wavcorP_wg_rbar_i[, vars := "frq_rbar_i"]

wavcorP_wg <- rbind(wavcorP_wg_r, wavcorP_wg_rbar_i)
wavcorP_wg[, units := "physical"]


# ----- genetic map

wavcorG_wg_r <- frq_rec_interp[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("frq", "r"), rm.boundary = F), by = .(rep, gen)]
wavcorG_wg_r[, vars := "frq_r"]

wavcorG_wg_rbar_i <- frq_rec_interp[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("frq", "rbar_i"), rm.boundary = F), by = .(rep, gen)]
wavcorG_wg_rbar_i[, vars := "frq_rbar_i"]

wavcorG_wg <- rbind(wavcorG_wg_r, wavcorG_wg_rbar_i)
wavcorG_wg[, units := "genetic"]

wavcor_wg <- rbind(wavcorP_wg, wavcorG_wg)

save(meanfrq, totalvars, totalcors, 
     wv_frq_rec_wg_all, 
     wavcor_wg, 
     file = gsub("_frqs.txt.gz", "_wavelet_results.RData", args[1]))

