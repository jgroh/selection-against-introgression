library(data.table)
library(cubature)
library(ggplot2)



# ===== Read and format simulation data =====

if(Sys.getenv("RSTUDIO") == "1"){
  
  source("~/workspace/gnomwav/R/multi_modwts.R")
  source("~/workspace/gnomwav/R/variance_decomp.R")
  source("~/workspace/gnomwav/R/theory.R")
  source("~/workspace/gnomwav/R/correlation_decomp.R")
  
  setwd("/Users/Jeff/workspace/selection-against-introgression/theory_and_simulations/results/neutral_sims/equilibrium/")
  n.sample <- 20
  haps <- fread("replicate0_hap_SNP_stat.txt", col.names = c("gen", "rep", "pos",  paste0("p0.", 1:n.sample), paste0("p1.", 1:n.sample), paste0("p2.", 1:n.sample)))
  frqs <- fread("replicate0_mean_SNP_stat.txt", col.names = c("gen", "rep", "pos", "p0", "p1", "p2"))
  ancestry <- fread("replicate0_hap_ancestry.txt", col.names = c("gen", "rep", "id", "left", "right", "source"))
  
} else {
  
  source("/Users/brogroh/gnomwav/R/multi_modwts.R")
  source("/Users/brogroh/gnomwav/R/variance_decomp.R")
  source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
  source("/Users/brogroh/gnomwav/R/theory.R")
  
  args <- commandArgs(trailingOnly = TRUE)
  n.sample <- as.numeric(args[1])
  ancestry <- fread(args[2],col.names = c("gen", "rep", "id", "left", "right", "source"))
  haps <- fread(args[3], col.names = c("gen", "rep", "pos",  paste0("p0.", 1:n.sample), paste0("p1.", 1:n.sample), paste0("p2.", 1:n.sample)))
  frqs <- fread(args[4], col.names = c("gen", "rep", "pos", "p0", "p1", "p2"))
}

haps <- melt(haps, id.vars = c('rep', 'gen', 'pos'), value.name = "allele")
haps[, c("pop", "id") := tstrsplit(variable, ".", fixed=T)][, variable := NULL]

haps[, Morgan := pos/1e8]
frqs[, Morgan := pos/1e8]

# subtract off generation to correspond with theory
haps[gen != 0, gen := gen-1]
frqs[gen != 0, gen := gen-1]
ancestry[gen != 0, gen := gen-1]

# check ancestry blocks make sense
#ggplot(ancestry[gen==10], aes(x = left, y = source)) + geom_line() + facet_wrap(~id)



# ===== Wavelet Variance of True Ancestry =====

xout = seq(2^-14, 1, by = 2^-14)

# ----- Method 1: directly measure ancestry at fine scales, most accurate comparison to theory (but less comparable to snp stat results due to differences in interpolation)

ancestry[, leftMorgan := left/1e8]
ancestry_grid1 <- ancestry[, approx(x = leftMorgan, y = source, xout = xout, method = 'constant', rule = 2), by = .(id,gen,rep)] 
setnames(ancestry_grid1, c('x','y'), c('Morgan', 'source'))


# wavelet variance 
ancestry_wv1 <- ancestry_grid1[, gnom_var_decomp(.SD,chromosome = NA, signals = "source"), by = .(rep,gen,id)]
ancestry_wv1 <- ancestry_wv1[grepl("d",level,fixed=T)]


# ----- Method 2. measure ancestry at SNP locations, then interpolate to the grid (same interpolation approach as for SNP stat)
# comparing this to method 1 above will show any distorting affect of the interpolation

sites <- frqs[abs(p0 - p1) >= 0.05 & gen == 0, pos, by = rep] 
# what is avg Morgan dist between snps
#1/length(sites$pos)

# get ancestry for haplotypes at SNP locations (constant interpolation here bc we're getting exact ancestry)
ancestry_snps <- ancestry[, approx(x = left, y = source, xout = sites$pos, method = 'constant'), by = .(id,gen,rep)]
setnames(ancestry_snps, c('x', 'y'), c('pos', 'source'))

# now interpolate to grid (linear interpolation should work better)
ancestry_snps[, Morgan := pos/1e8]
ancestry_grid2 <- ancestry_snps[, approx(x = Morgan, y= source, xout=xout, rule =2, method = 'linear'), by = .(id,gen,rep)]
setnames(ancestry_grid2, c('x','y'), c('Morgan','source'))

# wavelet variance 
ancestry_wv2 <- ancestry_grid2[, gnom_var_decomp(.SD,chromosome = NA,signals = "source"), by = .(rep,gen,id)]
ancestry_wv2 <- ancestry_wv2[grepl("d",level,fixed=T)]

ancestry_wv1[, ancestry_measure := "direct"]
ancestry_wv2[, ancestry_measure := "interpolated"]

hap_true_ancestry_wv <- rbind(ancestry_wv1, ancestry_wv2)

# average over haplotypes
hap_true_ancestry_wv <- hap_true_ancestry_wv[, .(single_hap = mean(variance.source)), by = .(rep, gen, level, ancestry_measure)]

# check match of both approaches with theory
#wvthry1.14 <- wavelet_variance_equilbrium(n.pop = 20000, n.sample=1, unit.scale = 2^-14, gen = c(1,10,100,1000), alpha= 0.5, level = 1:14)

# #can see how too coarse of SNP data results in fine scale variance being missed after interpolating ancestry to a fine grid
# ggplot(hap_true_ancestry_wv[gen!=0], aes(x = level, y = single_hap)) +
#   facet_wrap(~gen) +
#   geom_point() +
#   geom_line(aes(group = ancestry_measure, color = ancestry_measure)) +
#   geom_line(data = wvthry1.14[, .(level, gen, variance, id=11)], aes(x = level, y = variance), color = 'red', size=1.1) +
#   scale_x_discrete(breaks = c(paste0("d",1:14)), labels = as.character(-14:-1)) +
#   labs(x = expression(Scale: log[2](Morgan))) +
#   theme(aspect.ratio=1)



# ----- Wavelet Variance of mean true ancestry (mean of sample N=10) ----

# reformat so that haplotypes are in separate columns to take mean

ancestry_grid1[, pos := seq_len(.N), by = .(rep, gen, id)]
keycols <- c('gen', 'id', 'Morgan')
setkeyv(ancestry_grid1, cols = keycols)
ancestry_grid1[, id2 := rep(rep(paste0('id',1:n.sample), each = length(xout)), 5)]
ancestry_grid1_wide <- dcast(ancestry_grid1, rep + gen + Morgan ~ id2, value.var = 'source')


smpl_mean_true_ancestry <- ancestry_grid1_wide[, .(smpl_mean = rowSums(.SD)/n.sample), .SDcols = paste0('id',1:n.sample), by = .(rep, gen, Morgan)]

#ggplot(smpl_mean_true_ancestry, aes(x = Morgan, y = mean_ancestry)) + geom_line() + facet_wrap(~gen)

smpl_mean_true_ancestry_wv <- smpl_mean_true_ancestry[, gnom_var_decomp(.SD, signals = 'smpl_mean', chromosome = NA), by = .(rep, gen)]
setnames(smpl_mean_true_ancestry_wv, 'variance.smpl_mean',  'smpl_mean')
smpl_mean_true_ancestry_wv <- smpl_mean_true_ancestry_wv[level !='s14']
smpl_mean_true_ancestry_wv[, ancestry_measure := 'direct']

# combine results for output
true_ancestry_allWV <- rbind(melt(hap_true_ancestry_wv, variable.name = 'signal', 
     id.vars = c('rep','gen','level','ancestry_measure'), value.name = 'variance'),
     melt(smpl_mean_true_ancestry_wv, variable.name = 'signal',
     id.vars = c('rep','gen','level','ancestry_measure'), value.name = 'variance'))


#ggplot(true_ancestry_allWV[signal == 'smpl_mean'], aes(x = level, y = variance)) + geom_point() + facet_wrap(~gen)
# # check that this equals the weighted sum of parts from the single haplotype variance plus the wavelet covariance among haps
# # calculate wavelet covariance:
# true_ancestry_wc <- ancestry_grid1_wide[, cov_tbl(.SD, chromosome=NA, signals = c('id1','id2')), by = .(rep, gen)]
# setnames(true_ancestry_wc, "cov", "pair1_cov")
# 
# cmb34 <- combn(n.sample,2)
# # loop over remaining pairs
# for(pair in 2:ncol(cmb34)){
#   hap1 <- paste0('id',cmb34[1, pair])
#   hap2 <- paste0('id',cmb34[2, pair])
#   
#   # for pop 0
#   cov_tmp <- ancestry_grid1_wide[, cov_tbl(data=.SD, chromosome = NA, signals = c(hap1, hap2)), by = .(rep,gen)]
#   setnames(cov_tmp, "cov", paste0('pair', pair, '_cov'))
#   true_ancestry_wc <- merge(cov_tmp, true_ancestry_wc, by = c('rep', 'gen', 'level'))
# }
# 
# true_ancestry_wc <- melt(true_ancestry_wc, id.vars = c("rep", 'gen', 'level'), value.name = 'cov')
# true_ancestry_mean_wc <- true_ancestry_wc[, .(mean_cov = mean(cov)), by = .(rep, gen, level)]
# 
# # mean ancestry
# part1 <- ancestry_wv1[, .(part1 = mean(variance.source)), by = .(rep, gen, level)]
# true_ancestry_mean_data <- merge(true_ancestry_mean_wc, part1)
# true_ancestry_mean_data[, sum_of_parts := .9*mean_cov + (1/10)*part1]
# # 
# # ggplot(mean_true_ancestry_wv[grepl('d', level, fixed=T) & gen != 0], aes(x = level, y = variance.mean_ancestry)) +
# #   facet_wrap(~gen) + geom_point() +
# #   geom_line(group = 1) +
# #   geom_line(data=wvthry1.14_mean, aes(x = level, y = variance), color = 'red') +
# #   geom_line(data=true_ancestry_mean_wc[grepl('d',level, fixed=T)&gen!=0], aes(x = level, y = 0.9*mean_cov, group=1), color = 'orange') +
# #   geom_line(data = part1[gen != 0], aes(x = level, y = 0.1*part1, group = 1), color = 'purple') +
# #   geom_line(data=true_ancestry_mean_data[gen!=0], aes(x = level, y = sum_of_parts, group = 1), color = 'green')
# 
# # # black and green lines overlap perfectly, as expected.


# ===== Wavelet variance of snp stat computed on single haplotypes =====

# ----- Method 1. use allele frequencies from time of admixture to compute statistic 

# first ascertain sites with afd cutoff > 0.05 in gen 0
sites0 <- frqs[abs(p0 - p1) >= 0.05 & gen == 0, pos, by = rep] 
frqs_informative_gen0 <- frqs[gen == 0, .SD[pos %in% sites0[rep == .BY$rep, pos]], by = rep]
frqs_informative_gen0[, gen := NULL][, p2 := NULL]  # by setting gen to NULL here, when we merge with the haplotype table below, 
# this ensures that haplotypes are alwayspaired in the same row with the allele frqs from gen 0 

# merge haps and frqs to in order to compute ancestry statistic
haps_informative0 <- haps[allele != 2, .SD[pos %in% sites0[rep == .BY$rep, pos]], by = .(rep, gen)]
allsitedata0 <- merge(haps_informative0, 
                     frqs_informative_gen0, 
                     by = c("rep", "pos", "Morgan"))

allsitedata0[, h_hap := (allele - p0)/(p1-p0)]

# interpolate ancestry stat
haps_interp0 <- allsitedata0[, approx(x = Morgan, y = h_hap, xout=xout, rule=2), by = .(rep, gen, pop, id)]
setnames(haps_interp0, c("x", "y"), c("Morgan", "h_hap"))

# run wavelet variance decomp.
wv_haps0 <- haps_interp0[, gnom_var_decomp(.SD,chromosome = NA,signals = "h_hap"), by = .(rep,gen,pop,id)]
wv_haps0 <- wv_haps0[grepl("d", level, fixed=T)]
wv_haps0[, ascertainment := 'gen0']

# --- for plotting just within this script: reformat in order to compute variance correction
# prnt_wv_haps0 <- wv_haps0[pop!="p2", .(variance.prnt_haps = mean(variance.h_hap)), by = .(level, gen)]
# hyb_wv_haps0 <- wv_haps0[pop == "p2", .(variance.hyb_haps = mean(variance.h_hap)), by = .(gen, level)]
# wv_haps0.1 <- merge(prnt_wv_haps0, hyb_wv_haps0, by = c("level","gen"))
# wv_haps0.1[, adj_var := variance.hyb_haps - variance.prnt_haps]
# wv_haps0.1[, adj_propvar := adj_var/sum(adj_var), by = gen]
# wv_haps0.1[, ascertainment := "gen0"]
# 

# ----- Method 2. ascertain SNPs in contemporary generations
sites1 <- frqs[abs(p0 - p1) >= 0.05, pos, by = .(rep,gen)] 
frqs_informative1 <- frqs[, .SD[pos %in% sites1[rep == .BY$rep & gen == .BY$gen, pos]], by = .(rep, gen)]

# merge haps and frqs to in order to compute ancestry statistic
# note that haps contains only segregating sites in the selected haps, whereas frqs contains all of the sites
haps_informative1 <- haps[allele != 2, .SD[pos %in% sites1[rep == .BY$rep & gen == .BY$gen, pos]], by = .(rep, gen)]

allsitedata1 <- merge(haps_informative1, 
                      frqs_informative1, 
                      by = c("rep", "gen","pos", "Morgan"))

allsitedata1[, h_hap := (allele - p0)/(p1-p0)]

# interpolate ancestry stat
haps_interp1 <- allsitedata1[, approx(x = Morgan, y = h_hap, xout=xout, rule=2), by = .(rep, gen, pop, id)]
setnames(haps_interp1, c("x", "y"), c("Morgan", "h_hap"))

# run wavelet variance decomp
wv_haps1 <- haps_interp1[, gnom_var_decomp(.SD,chromosome = NA,signals = "h_hap"), by = .(rep,gen,pop,id)]
wv_haps1 <- wv_haps1[grepl("d", level, fixed=T)]
wv_haps1[, ascertainment:='contemporary']

# ----- output from script:

all_wv_haps <- rbind(wv_haps0, wv_haps1)
all_wv_haps <- all_wv_haps[, .(mean_wv.h_hap = mean(variance.h_hap)), by = .(rep, gen, pop, level, ascertainment)]

# --- for plotting just within this script: reformat in order to compute variance correction
# prnt_wv_haps1 <- wv_haps1[pop!="p2", .(variance.prnt_haps = mean(variance.h_hap)), by = .(level, gen)]
# hyb_wv_haps1 <- wv_haps1[pop == "p2", .(variance.hyb_haps = mean(variance.h_hap)), by = .(gen, level)]
# wv_haps1.1 <- merge(prnt_wv_haps1, hyb_wv_haps1, by = c("level","gen"))
# wv_haps1.1[, adj_var := variance.hyb_haps - variance.prnt_haps]
# wv_haps1.1[, adj_propvar := adj_var/sum(adj_var), by = gen]
# wv_haps1.1[, ascertainment := "contemporary"]
# 
# all_wv_haps.1 <- rbind(wv_haps0.1, wv_haps1.1)

# # check the accuracy of the correction

# wvthry1.14[, propvar := variance/sum(variance), by = gen]

# # raw variance is much off

# ggplot(all_wv_haps.1[gen %in% c(1,10,100,1000)], aes(x = level)) + facet_wrap(~gen) + theme(aspect.ratio = 1) +
#   geom_point(aes(y = variance.hyb_haps, shape = ascertainment), color = 'black') + geom_point(aes(y = variance.prnt_haps, shape = ascertainment), color = 'blue') +
#   scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
#   labs(x = expression(Scale: log[2] (Morgan)), y = 'Variance') +
#   geom_point(aes(y = variance.hyb_haps - variance.prnt_haps, shape=ascertainment), color = 'orange') +
#   geom_point(data = wvthry1.14, aes(y = variance), color = 'red')
# # 
# # but proportion-wise looks great
# ggplot(all_wv_haps.1[gen %in% c(1,10,100,1000)], aes(x = level)) + facet_wrap(~gen) + theme(aspect.ratio = 1) +
#   scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
#   labs(x = expression(Scale: log[2] (Morgan)), y = 'Proportion of variance') +
#   geom_point(aes(y = adj_propvar, shape = ascertainment), color = 'orange') +
#   geom_point(data = wvthry1.14, aes(y = propvar), color = 'red')
  

# ===== Wavelet Variance of sample mean snp statistic (N=10) ======

# reformat in order to compute mean
haps_interp0[, id2 := paste0('id', id), by = .(rep, gen, pop, id)]
haps_interp0_wide <- dcast(haps_interp0, rep + gen + pop + Morgan ~ id2, value.var = 'h_hap')
smpl_mean_h0 <- haps_interp0_wide[, .(smpl_mean_h = rowSums(.SD)/n.sample), .SDcols = paste0('id',1:n.sample), by = .(rep,gen,pop,Morgan)]
smpl_mean_h0[, ascertainment := 'gen0']

# reformat in order to compute mean
haps_interp1[, id2 := paste0('id', id), by = .(rep, gen, pop, id)]
haps_interp1_wide <- dcast(haps_interp0, rep + gen + pop + Morgan ~ id2, value.var = 'h_hap')
smpl_mean_h1 <- haps_interp1_wide[, .(smpl_mean_h = rowSums(.SD)/n.sample), .SDcols = paste0('id',1:n.sample), by = .(rep,gen,pop,Morgan)]
smpl_mean_h1[, ascertainment := 'contemporary']

# ---- output from script: wavelet variance of sample mean
smpl_mean_wv0 <- smpl_mean_h0[, gnom_var_decomp(.SD, chromosome = NA, signals = 'smpl_mean_h'), by = .(rep, gen, pop, ascertainment)]
smpl_mean_wv1 <- smpl_mean_h1[, gnom_var_decomp(.SD, chromosome = NA, signals = 'smpl_mean_h'), by = .(rep, gen, pop, ascertainment)]

smpl_mean_wv <- rbind(smpl_mean_wv0, smpl_mean_wv1)

# # look at wavelet variance of mean, can see the small bump for p2
# ggplot(smpl_mean_wv, aes(x = level, y = variance.smpl_mean_h, color = pop)) + facet_wrap(~gen) + geom_point()


# ===== Wavelet Covariances for pairs of haps =====
# Again, just use SNPs ascertained in gen0 for this

# we'll loop over all unordered pairs of individuals from each population. Initiate the data structure by calculating covariances for the first pair of individuals
h_wc <- haps_interp0_wide[, cov_tbl(data=.SD, chromosome = NA, signals = c('id1', 'id2')), by = .(rep,gen,pop)]
setnames(h_wc, 'cov', paste0('pair', 1, '_cov'))

cmb34 <- combn(n.sample,2)
# loop over remaining pairs
for(pair in 2:ncol(cmb34)){
  hap1 <- paste0('id',cmb34[1, pair])
  hap2 <- paste0('id',cmb34[2, pair])
  
  # for pop 0
  cov_tmp <- haps_interp0_wide[, cov_tbl(data=.SD, chromosome = NA, signals = c(hap1, hap2)), by = .(rep,gen,pop)]
  setnames(cov_tmp, "cov", paste0('pair', pair, '_cov'))
  h_wc <- merge(cov_tmp, h_wc, by = c('rep', 'gen','pop','level'))
}

# ----- output from script: average wavelet covariance among pairs of haplotypes 
h_wc <- melt(h_wc, id.vars = c("rep", 'gen', 'pop', 'level'), value.name = 'cov')
h_wc <- h_wc[, .(mean_cov = mean(cov)), by = .(rep, gen, level, pop)]
h_wc <- h_wc[grepl('d', level, fixed=T)]

# # check that wavelet variance of sample mean equals the sum of two parts
# part1 <- wv_haps0[, .(part1 = 0.1*mean(variance.h_hap)), by = .(rep, gen, pop, level)]
# part2 <- h_wc[, .(part2 = 0.9*mean_cov), by = .(rep, gen, level, pop)]
# both_parts <- merge(part1, part2)
# both_parts[, total := sum(part1, part2), by = .(rep, gen, pop, level)]
# 
# ggplot(mean_h_wv[gen != 0 & pop == 'p2' & level!='s14'], aes(x = level, y = variance.mean_h)) + facet_wrap(~gen) + geom_point() + 
#   geom_line(data=part1[gen !=0 & level!='s14' & pop=='p2'], aes(x = level, y = part1, group=1), color = 'orange') + 
#   geom_line(data=part2[gen !=0 & level !='s14'&pop=='p2'], aes(x = level, y = part2, group = 1), color = 'purple') + 
#   geom_line(data = both_parts[gen !=0 & pop == 'p2'], aes(x = level, y = total, group = 1), color = 'green')
# # looks good!


# ===== Wavelet Variance of Mean SNP Stat =====


# ----- Method 1. use allele frequencies from time of admixture for parental pops to compute statistic 
# this table has contemporary allele frequencies for p2 but allele frequencies from gen 0 from the parental pops
frqs_informative_prnt0_p2contemp <- merge(frqs_informative1[, .(rep, gen, pos, p2, Morgan)], frqs_informative_gen0[, .(rep, pos, p0, p1, Morgan)])

# interpolate 
frqs_informative_prnt0_p2contemp[, h_frq := (p2-p0)/(p1-p0)]
frqs_interp0 <- frqs_informative_prnt0_p2contemp[, approx(x = Morgan, y = h_frq, xout = xout, rule = 2), by = .(rep, gen)]
setnames(frqs_interp0, c('x','y'), c('Morgan', 'h_frq'))

# wavelet variance
frqs_wv0 <- frqs_interp0[, gnom_var_decomp(.SD, signals = 'h_frq', chromosome = NA), by = .(rep, gen)]
frqs_wv0 <- frqs_wv0[grepl('d', level, fixed=T)]
frqs_wv0[, ascertainment := 'gen0']


# ----- Method 2. use allele contemporary allele frequencies of parental pops to compute statistic 
frqs_informative1[, h_frq := (p2-p0)/(p1-p0) ]
frqs_interp1 <- frqs_informative1[, approx(x=Morgan, y = h_frq, xout=xout, rule = 2), by = .(rep, gen)]
setnames(frqs_interp1, c('x','y'), c('Morgan', 'h_frq'))

frqs_wv1 <- frqs_interp1[, gnom_var_decomp(.SD, signals = 'h_frq', chromosome = NA), by = .(rep, gen)]
frqs_wv1 <- frqs_wv1[grepl('d', level, fixed=T)]
frqs_wv1[, ascertainment := 'contemporary']

all_wv_frqs <- rbind(frqs_wv0, frqs_wv1)

# # examine effect of time of ascertainment
# ggplot(all_wv_frqs[gen %in% c(0,10,100,1000)], aes(x = level, y = variance.h_frq)) + facet_wrap(~gen, scales = 'free_y') + geom_point(aes(color = ascertainment)) +
#   theme(aspect.ratio = 1) +
#   scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
#   labs(x = expression(Scale: log[2] (Morgan)),
#         y = "Variance")


# ----- output from script: single hap, sample mean, and population mean wavelet variances
# reformatting
setnames(smpl_mean_wv, 'variance.smpl_mean_h', 'smpl_mean')
smpl_mean_wv <- melt(smpl_mean_wv, measure.vars = 'smpl_mean', value.name = 'variance', variable.name = 'signal')

all_wv_frqs[, pop := 'p2']
setnames(all_wv_frqs, 'variance.h_frq', 'pop_mean')
all_wv_frqs <- melt(all_wv_frqs, measure.vars = 'pop_mean', value.name = 'variance', variable.name = 'signal')

setnames(all_wv_haps, 'mean_wv.h_hap', 'single_hap')
all_wv_haps <- melt(all_wv_haps, measure.vars = 'single_hap', value.name = 'variance', variable.name = 'signal')

allSNPWV <- rbind(rbind(all_wv_haps, all_wv_frqs), smpl_mean_wv)

save(true_ancestry_allWV, allSNPWV, h_wc, file = gsub('_hap_ancestry.txt','_wavelet_results.RData',args[2]))
