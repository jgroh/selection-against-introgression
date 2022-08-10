library(data.table)
library(cubature)
library(ggplot2)
library(magrittr)


# ----- read and format simulation data -----

if(Sys.getenv("RSTUDIO") == "1"){
  
  source("~/workspace/gnomwav/R/multi_modwts.R")
  source("~/workspace/gnomwav/R/variance_decomp.R")
  source("~/workspace/gnomwav/R/theory.R")
  source("~/workspace/gnomwav/R/correlation_decomp.R")
  
  setwd("/Users/Jeff/workspace/selection-against-introgression/theory_and_simulations/results/snp_stat_const_rec_n10000/")
  n.sample <- 10
  haps <- fread("replicate0_haps.txt", col.names = c("gen", "rep", "pos",  paste0("p0.", 1:n.sample), paste0("p1.", 1:n.sample), paste0("p2.", 1:n.sample)))
  frqs <- fread("replicate0_frqs.txt", col.names = c("gen", "rep", "pos", "p0", "p1", "p2"))
  ancestry <- fread("replicate0_ancestry.txt", col.names = c("gen", "rep", "id", "left", "right", "source"))
  
} else {
  
  source("/Users/brogroh/gnomwav/R/multi_modwts.R")
  source("/Users/brogroh/gnomwav/R/variance_decomp.R")
  source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
  source("/Users/brogroh/gnomwav/R/theory.R")
  
  args <- commandArgs(trailingOnly = TRUE)
  n.sample <- as.numeric(args[1])
  haps <- fread(args[2], col.names = c("gen", "rep", "pos",  paste0("p0.", 1:n.sample), paste0("p1.", 1:n.sample), paste0("p2.", 1:n.sample)))
  frqs <- fread(args[3], col.names = c("gen", "rep", "pos", "p0", "p1", "p2"))
}

haps <- melt(haps, id.vars = c('rep', 'gen', 'pos'), value.name = "allele")
haps[, c("pop", "id") := tstrsplit(variable, ".", fixed=T)][, variable := NULL]

haps[, Morgan := pos/1e8]
frqs[, Morgan := pos/1e8]



# ===== Compare wavelet variance for true ancestry to expectation =====


sites <- frqs[abs(p0 - p1) >= 0.05 & gen == 0, pos, by = rep] 
# what is avg Morgan dist between snps
#1/length(sites$pos)

wvthry1.14 <- wavelet_variance_equilbrium(n.pop = 20000, n.sample=1, unit.scale = 2^-14, gen = c(10,100,1000), alpha= 0.5, level = 1:14)
#wvthry1.7 <- wavelet_variance_equilbrium(n.pop = 20000, n.sample=1, unit.scale = 2^-7, gen = c(10,100,1000), alpha= 0.5, level = 1:7)
#wvthry1.10 <- wavelet_variance_equilbrium(n.pop = 20000, n.sample=1, unit.scale = 2^-10, gen = c(10,100,1000), alpha= 0.5, level = 1:10)
#wvthry1.11 <- wavelet_variance_equilbrium(n.pop = 20000, n.sample=1, unit.scale = 2^-11, gen = c(10,100,1000), alpha= 0.5, level = 1:11)

xout = seq(2^-14, 1, by = 2^-14)

# get true ancestry for haplotypes at these sites (constant interpolation here bc we're getting exact ancestry)
ancestry_snps <- ancestry[, approx(x = left, y = source, xout = sites$pos, method = 'constant'), by = .(id,gen,rep)]
setnames(ancestry_snps, c('x', 'y'), c('pos', 'source'))

# now interpolate to grid (linear interpolation works better)
ancestry_snps[, Morgan := pos/1e8]
ancestry_interp <- ancestry_snps[, approx(x = Morgan, y= source, xout=xout, rule =2, method = 'linear'), by = .(id,gen,rep)]
setnames(ancestry_interp, c('x','y'), c('Morgan','ancestry'))

#ggplot(ancestry_interp[id==339632], aes(x = Morgan, y = ancestry)) + geom_point() + theme(aspect.ratio = 1)

# wavelet variance 
ancestry_wv <- ancestry_interp[, gnom_var_decomp(.SD,chromosome = NA,signals = "ancestry"), by = .(rep,gen,id)]
ancestry_wv <- ancestry_wv[grepl("d",level,fixed=T)]


# ggplot(ancestry_wv[gen == 1000], aes(x = level, y = variance.ancestry, group = id)) + 
#   geom_point() + 
#   #facet_wrap(~gen) + 
#   geom_line() +
#   #geom_line(data = wvthry1.7[gen == 1000, .(level, gen, variance, id=11)], aes(x = level, y = variance), color = 'red', size=1.1) + 
#   #geom_line(data = wvthry1.10[gen == 1000, .(level, gen, variance, id=11)], aes(x = level, y = variance), color = 'red', size=1.1) + 
#   #geom_line(data = wvthry1.11[gen == 1000, .(level, gen, variance, id=11)], aes(x = level, y = variance), color = 'red', size=1.1) + 
#   geom_line(data = wvthry1.14[gen == 1000, .(level, gen, variance, id=11)], aes(x = level, y = variance), color = 'red', size=1.1) + 
#   #scale_x_discrete(breaks = c(paste0("d",1:10)), labels = as.character(-10:-1)) + 
#   scale_x_discrete(breaks = c(paste0("d",1:14)), labels = as.character(-14:-1)) + 
#   #scale_x_discrete(breaks = c(paste0("d",1:7)), labels = as.character(-7:-1)) + 
#   labs(x = expression(Scale: log[2](Morgan))) + 
#   theme(aspect.ratio=1)



# ===== Wavelet variance of haplotypes =====

# first ascertain sites with afd cutoff > 0.5 in gen 0
sites <- frqs[abs(p0 - p1) >= 0.05 & gen == 0, pos, by = rep] 
frqs_informative_gen0 <- frqs[gen == 0, .SD[pos %in% sites[rep == .BY$rep, pos]], by = .(rep, gen)] # note, here option to ascertain sites variable in contemporary gens
frqs_informative_gen0[, gen := NULL][, p2 := NULL]  # if ascertaining in contemporary gens, remove the g:=NULL 

# merge haps and frqs to in order to compute ancestry statistic
haps_informative <- haps[allele != 2, .SD[pos %in% sites[rep == .BY$rep, pos]], by = .(rep, gen)]
allsitedata <- merge(haps_informative, 
                     frqs_informative_gen0, 
                     by = c("rep", "pos", "Morgan"))

allsitedata[, h_hap := (allele - p0)/(p1-p0)]

# interpolate ancestry stat
haps_interp <- allsitedata[, approx(x = Morgan, y = h_hap, xout=xout, rule=2), by = .(rep, gen, pop, id)]
setnames(haps_interp, c("x", "y"), c("Morgan", "h_hap"))

# run wavelet variance decomp
wv_haps <- haps_interp[, gnom_var_decomp(.SD,chromosome = NA,signals = "h_hap"), by = .(rep,gen,pop,id)]
wv_haps <- wv_haps[grepl("d", level, fixed=T)]

# reformat so that single chrom variance for parental chromosomes and admixed chromosomes are in separate columns
prnt_wv_haps <- wv_haps[pop!="p2", .(variance.prnt_haps = mean(variance.h_hap)), by = .(level, gen)]
hyb_wv_haps <- wv_haps[pop == "p2", .(variance.hyb_haps = mean(variance.h_hap)), by = .(gen, level)]

wv_haps2 <- merge(prnt_wv_haps, hyb_wv_haps, by = c("level","gen"))
wv_haps2[, adj_var := variance.hyb_haps - variance.prnt_haps]

#wv_haps2[, propvar.hyb_haps := variance.hyb_haps/sum(variance.hyb_haps), by = .(gen)]
#wv_haps2[, propvar.prnt_haps := variance.prnt_haps/sum(variance.prnt_haps), by = .(gen)]
wv_haps2[, adj_propvar := adj_var/sum(adj_var), by = gen]
wvthry1.14[, propvar := variance/sum(variance), by = gen]


ggplot(wv_haps2[gen %in% c(10,100,1000)], aes(x = level)) + facet_wrap(~gen) + theme(aspect.ratio = 1) + 
  geom_point(aes(y = variance.hyb_haps), color = 'black') + geom_point(aes(y = variance.prnt_haps), color = 'blue') + 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Variance') + 
  geom_point(aes(y = variance.hyb_haps - variance.prnt_haps), color = 'orange') + 
  geom_point(data = wvthry1.14, aes(y = variance), color = 'red')
  
ggplot(wv_haps2[gen %in% c(10,100,1000)], aes(x = level)) + facet_wrap(~gen) + theme(aspect.ratio = 1) + 
  scale_x_discrete(breaks = paste0("d",1:14), labels = as.character(-14:-1)) +
  labs(x = expression(Scale: log[2] (Morgan)), y = 'Proportion of variance') + 
  geom_point(aes(y = adj_propvar), color = 'orange') + 
  geom_point(data = wvthry1.14, aes(y = propvar), color = 'red')
  









# ----- calculate ancestry statistic at all loci for mean -----
frqs[, h_frq := (p2-p0)/(p1-p0)] # approximates the indicator for p1 ancestry







# interpolate ancestry frq stat 
frqs_informative_gen0 <- frqs[, .SD[pos %in% sites[rep == .BY$rep, pos]], by = .(rep, gen)]
#frqs_informative_gen0[, gen := NULL][, p2 := NULL][, h_frq := NULL]

# ****** NEED TO REVISIT are we calculating h using allele freq of contemporary pops
frqs_interp <- frqs_informative[, approx(x = Morgan, y = h_frq, xout=xout, rule=2), by = .(rep, gen)]
setnames(frqs_interp, c("x", "y"), c("Morgan", "h_frq"))

# ----- Wavelet Variance decomposition of frequency stat -----
wv_frq <- frqs_interp[, gnom_var_decomp(.SD, chromosome = NA, signals = "h_frq"), by = .(rep, gen)]
wv_frq <- wv_frq[grepl("d", level)]


# thry_frq <- wavelet_variance_equilbrium(n.pop = 20000, n.sample = 20000, unit.scale = 2^-10, gen = c(10,100,500,1000),level = 1:10, alpha = 0.5)
# # unclear currently how theory for frq of known admixture relates to mean of snp stat
# 
# ggplot(wv_frq[gen !=0], aes(x = level, y = variance.h_frq)) + geom_point() +
#  geom_line(aes(group = rep)) + facet_wrap(~gen, scales = "free_y") +
#  scale_x_discrete(breaks = paste0("d",1:10), labels = as.character(-10:-1)) +
#  labs(x = expression(Scale: log[2] (Morgan)),
#       y = "Variance") +
#  geom_line(data=thry_frq, aes(x = level, y = variance), color = "red") +
#  theme(aspect.ratio=1)


# ===== Wavelet Variance Decomp on individual haps =====

#haps[allele>1] # what to do about these? maybe just ignore non-biallelic sites for now
haps_informative <- haps[allele != 2, .SD[pos %in% sites[rep == .BY$rep, pos]], by = .(rep, gen)]

# combine haps and frqs to get informative sites (sites ascertained at time of admixture)
allsitedata <- merge(haps_informative, 
                frqs_informative_gen0, 
                by = c("rep", "gen", "pos", "Morgan"))

# used to double check whether parental frqs remain constant across gens
#allsitedata[pos == 206380 & id == 1]

# compute ancestry stat at informative loci
allsitedata[, h_hap := (allele - p0)/(p1-p0)]

# interpolate ancestry stat
haps_interp <- allsitedata[, approx(x = Morgan, y = h_hap, xout=xout, rule=2), by = .(rep, gen, pop, id)]
setnames(haps_interp, c("x", "y"), c("Morgan", "h_hap"))

#ggplot(haps_interp[pop != "p2" & gen == 0], aes(x = Morgan, y = h_hap, color = pop)) + geom_point() + 
#  facet_wrap(~id) + geom_line()


# ----- compute Wavelet Variance Decomp for individual haps -----
wv_haps <- haps_interp[, gnom_var_decomp(.SD,chromosome = NA,signals = "h_hap"), by = .(rep,gen,pop,id)]
wv_haps <- wv_haps[grepl("d", level, fixed=T)]
wv_haps[, propvar := variance.h_hap/sum(variance.h_hap), by = .(rep,gen,pop,id)]

# ggplot(wv_haps[gen != 0 & grepl("d", level, fixed=T), .(v = mean(variance.h_hap)), by = .(gen, pop, level)],
#        aes(x = level, y = v)) +
#   geom_point(aes(group = pop)) +
#   geom_line(aes(color = pop, group = pop)) +
#   facet_wrap(~gen)+
#   labs(x = expression(Scale: log[2] (Morgan)),
#        y = "Variance") +
#   theme(aspect.ratio=1)

# using allele frequencies at time of admixture, the wavelet variance of contemporary parental haps remains same
# ggplot(wv_haps[pop != "p2" & grepl("d", level, fixed=T),
#               .(v = mean(variance.h_hap)), by = .(pop,gen, level)],
#       aes(x = level, y = v)) + geom_point() +
#  geom_line(aes(group = interaction(gen,pop), color = as.factor(gen)))


# reformat so that single chrom variance for parental chromosomes and admixed chromosomes are in separate columns
prnt_wv_haps <- wv_haps[pop!="p2", .(variance.prnt_haps = mean(variance.h_hap)), by = .(level, gen)]

hyb_wv_haps <- wv_haps[pop == "p2", .(variance.hyb_haps = mean(variance.h_hap)), by = .(gen, level)]

wv_haps2 <- merge(prnt_wv_haps, hyb_wv_haps, by = c("level","gen"))



# true ancestry

# ancestry at the same sites as the SNPs
ancestry_snps <- ancestry[, approx(x = left, y = source, xout = sites$pos, method = 'constant'), by = .(id,gen,rep)]
setnames(ancestry_snps, c('x', 'y'), c('pos', 'source'))

# now interpolate to genetic distance grid
ancestry_snps[, Morgan := pos/1e8]
ancestry_interp <- ancestry_snps[, approx(x = Morgan, y= source, xout=xout, rule =2), by = .(id,gen,rep)]
setnames(ancestry_interp, c('x','y'), c('Morgan','ancestry'))


ancestry_wv <- ancestry_interp[, gnom_var_decomp(.SD,chromosome = NA,signals = "ancestry"), by = .(rep,gen,id)]
ancestry_wv <- ancestry_wv[grepl("d",level,fixed=T)]

wvthry <- wavelet_variance_equilbrium(n.pop = 20000, n.sample=1, unit.scale = 2^-7, gen = c(0,10,100,1000), alpha= 0.5, level = 1:7)

ggplot(ancestry_wv[gen %in% c(0,10,100,1000)], aes(x = level, y = variance.ancestry, group = id)) + 
  geom_point() + facet_wrap(~gen) + geom_line() +
  geom_line(data = wvthry[gen %in% c(0,10,100,1000), .(level, gen, variance, id=11)], aes(x = level, y = variance), color = 'red', size=1.1) + 
  #geom_line(data = wv_haps[gen %in% c(0,10,100,1000) & pop== 'p2'], aes(x = level, y = variance.h_hap, group = id), color = 'green') + 
  scale_x_discrete(breaks = c(paste0("d",1:7)), labels = as.character(-7:-1)) + 
  labs(x = expression(Scale: log[2](Morgan))) + 
  theme(aspect.ratio=1)

mean_ancestry_wv <- ancestry_wv[, .(ancestry_var = mean(variance.ancestry)), by = .(gen,level)]
mean_snp_stat_wv <- wv_haps[pop=='p2', .(mean_snp_stat_var = mean(variance.h_hap)), by = .(gen, level)]
mean_prnt_var <- wv_haps[gen == 0 & pop != 'p2', .(prnt_hap_var = mean(variance.h_hap)), by = level]

all_mean_wvs <- merge(mean_prnt_var, mean_wvs, by = 'level')

ggplot(all_mean_wvs) + facet_wrap(~gen) + theme(aspect.ratio = 1) + 
  geom_line(aes(x = level, y = ancestry_var, group = 1), color = 'black') +
  geom_line(aes(x = level, y = mean_snp_stat_var, group = 1), color = 'red') + 
  geom_line(aes(x = level, y = mean_snp_stat_var - prnt_hap_var, group = 1), color = 'green') + 
  scale_x_discrete(breaks = c(paste0("d",1:7)), labels = as.character(-7:-1)) + 
  labs(x = expression(Scale: log[2](Morgan))) 

str(all_mean_wvs)
all_mean_wvs[, mean_snp_stat_wv - prnt_hap_var]
# combine frq and hap data for output
allWV <- merge(wv_haps2, wv_frq)


#fwrite(allWV, file=paste0("results/admix_snp_stat/replicate",allWV[1,rep], "_wv_results.txt"),quote=F, sep="\t")


#wvtheory <- wavelet_variance_equilbrium(n.pop=20000, n.sample = 1, unit.scale = 2^-7, level = 1:7, alpha = 0.5, gen = c(10,100,500,1000))

#compare wv theory with actual wv on snp stat, shows there is much more noise in the snp-based statistic compared to the expectation for true ancestry state
# ggplot(hyb_wv_haps[gen != 0], aes(x = level, y = variance.hyb_haps)) + geom_point() +
#   geom_line(group = 1) + facet_wrap(~gen) +
#   geom_line(data = wvtheory,
#             aes(x = level, y = variance),
#             color  = 'red') +
#   theme(aspect.ratio=1) +
#   scale_x_discrete(breaks = paste0("d",1:7), labels = as.character(-7:-1)) +
#   labs(x = expression(Scale: log[2] (Morgan)),
#        y = "Variance")




# ===== Calculate Wavelet Covariances =====
haps_interp[, pos := seq_len(.N), by = .(rep, gen, pop, id)]

# reformat so that haplotypes are in separate columns
haps_interp2 <- dcast(haps_interp[, .(rep,gen,pop,id,Morgan,h_hap, pos=as.numeric(gsub("d","",pos)))], rep + gen + pop + pos ~ id, value.var = 'h_hap')

# average covariance over all pairs of individuals from same pop
# This gives a table of the pairs we'll average over.
cmb34 <- combn(n.sample, 2)

# we'll loop over all unordered pairs of individuals from each population. Initiate the data structure by calculating covariances for the first pair of individuals
covs_00 <- haps_interp2[pop == 'p0', cov_tbl(data=.SD, chromosome = NA, signals = c('1', '2')), by = .(rep,gen,pop)]
setnames(covs_00, 'cov', paste0('pair', 1, '_cov'))

covs_11 <- haps_interp2[pop == 'p1', cov_tbl(data=.SD, chromosome = NA, signals = c('1', '2')), by = .(rep,gen,pop)]
setnames(covs_11, 'cov', paste0('pair', 1, '_cov'))

covs_22 <- haps_interp2[pop == 'p2', cov_tbl(data=.SD, chromosome = NA, signals = c('1', '2')), by = .(rep,gen,pop)]
setnames(covs_22, 'cov', paste0('pair', 1, '_cov'))

for(pair in 2:ncol(cmb34)){
  id1 <- as.character(cmb34[1, pair])
  id2 <- as.character(cmb34[2, pair])
  
  # for pop 0
  covs0 <- haps_interp2[pop == 'p0', cov_tbl(data=.SD, chromosome = NA, signals = c(id1, id2)), by = .(rep,gen,pop)]
  setnames(covs0, "cov", paste0('pair', pair, '_cov'))
  covs_00 <- merge(covs0, covs_00, by = c('rep', 'gen', 'pop', 'level'))
  
  # for pop 1
  covs1 <- haps_interp2[pop == 'p1', cov_tbl(data=.SD, chromosome = NA, signals = c(id1, id2)), by = .(rep,gen,pop)]
  setnames(covs1, "cov", paste0('pair', pair, '_cov'))
  covs_11 <- merge(covs1, covs_11, by = c('rep', 'gen', 'pop', 'level'))
  
  # for pop 2
  covs2 <- haps_interp2[pop == 'p2', cov_tbl(data=.SD, chromosome = NA, signals = c(id1, id2)), by = .(rep,gen,pop)]
  setnames(covs2, "cov", paste0('pair', pair, '_cov'))
  covs_22 <- merge(covs2, covs_22, by = c('rep', 'gen', 'pop', 'level'))
  
}

cov_00_mean <- melt(covs_00, id.vars = c("rep", 'gen', 'pop', 'level'), value.name = 'cov')
cov_11_mean <- melt(covs_11, id.vars = c("rep", 'gen', 'pop', 'level'), value.name = 'cov')
cov_22_mean <- melt(covs_22, id.vars = c("rep", 'gen', 'pop', 'level'), value.name = 'cov')

cov_same_pop_mean <- rbind(rbind(cov_00_mean, cov_11_mean))


# ggplot(cov_same_pop_mean[level !="s7"], aes(x = level, y = cov,col = pop, group = interaction(pop, variable))) + 
#   geom_point() + facet_wrap(~gen, scales = 'free_y') + 
#   geom_line()

wv_cov_terms <- cov_same_pop_mean[, .(prnt_wav_cov = mean(cov)), by = .(level, rep, gen)]

allWV <- merge(wv_cov_terms, allWV, by = c("rep", "gen", "level"))

#ggplot(allWV, aes(x = level, y = prnt_wav_cov)) + geom_point() + facet_wrap(~ gen)


# ========= Conditioning on haplotype state at time of admixture =====
# here we'll try to calculate the product hl*hl' conditioning on source populations for locus 1 and locus 2
# these will be used to try to evaluate the various terms in the integral for the wavelet variance (see theory in supp mat)

# reformat so locus positions are in separate columns to make calculation easier
# haps_interp[, pos := paste0("d", pos)]
# haps_interp3 <- dcast(haps_interp, rep + gen + pop + id ~ pos, value.var = 'h_hap')
# 
# cmb34 <- combn(n.sample, 2)
# cmb5 <- expand.grid(1:n.sample, 1:n.sample)
# 
# L <- 128
# all_hh <- data.table()
# 
# for (g in 0){ # look at the product for parental chromosomes in gen zero, should not make a huge difference across gens
#   
#   # these vectors will hold the mean of hl * hl' for the same chromosome with entries corresponding to distance between l and l'
#   hh_00_same <- vector()
#   hh_11_same <- vector()
#   
#   # these vectors will hold the mean of hl * hl' for different chromosomes with entries corresponding to distance between l and l'
#   hh_00_diff <- vector()
#   hh_11_diff <- vector()
#   hh_01_diff <- vector()
#   
#   for(dist in 1:(L-1)){ # loop over possible distances between l and l'
#     
#     # this vector will hold values of the product and we'll take it's mean at the end
#     vals_00_dist <- vector()
#     vals_11_dist <- vector()
#     vals_01_dist <- vector()
#     
#     # we could average over all starting positions but this seemed to take forever. plus we'll be averaging over replicate sims
#     l1 <- 1
#     l2 <- l1 + dist
#     cols <- paste0('d', c(l1,l2))
#     
#     hh_00_same[dist] <- haps_interp3[pop=="p0" & gen == g, mean(get(cols[1])*get(cols[2]))]
#     hh_11_same[dist] <- haps_interp3[pop=="p1" & gen == g, mean(get(cols[1])*get(cols[2]))]
#     
#     np <- ncol(cmb34)
#     for(h in 1:np){
#       # ids of haplotypes
#       i <- cmb34[1, h]
#       j <- cmb34[2, h]
#       
#       vals_00_dist  <- c(vals_00_dist, haps_interp3[id == i & gen == g & pop == "p0", get(cols[1])]*haps_interp3[id == j & gen == g & pop == "p0", get(cols[2])])
#       vals_11_dist <- c(vals_11_dist, haps_interp3[id == i & gen == g & pop == "p1", get(cols[1])]*haps_interp3[id == j & gen == g & pop == "p1", get(cols[2])])
#     }
#     
#     hh_00_diff[dist] <- mean(vals_00_dist)
#     hh_11_diff[dist] <- mean(vals_11_dist)
#     
#     # for term hh_01_diff, loop over combinations one drawn from each pop
#     np2 <- nrow(cmb5)
#     for(h in 1:np2){
#       i <- cmb5[h, 1]
#       j <- cmb5[h, 2]
#       
#       vals_01_dist <- c(vals_01_dist, haps_interp3[id == i & gen == g & pop == "p0", get(cols[1])]*haps_interp3[id == j & gen == g & pop == "p1", get(cols[2])])
#     }
#     
#     hh_01_diff[dist] <- mean(vals_01_dist)
#     
#   }
#   
#   hh <- data.table(l1l2_dist = 1:(L-1), gen=g, rep= sites[1, rep], hh_00_same = hh_00_same, hh_11_same=hh_11_same, hh_00_diff=hh_00_diff, hh_11_diff=hh_11_diff, hh_01_diff = hh_01_diff) 
#   all_hh <- rbind(all_hh, hh)
# }
# 
# #ggplot(melt(all_hh, id.vars = c("l1l2_dist", "gen", "rep")), aes(x = l1l2_dist, y = value, color = variable)) + geom_point()
# save(allWV, all_hh, file = paste0("results/admix_snp_stat/replicate",allWV[1,rep], "_wv_results.RData"))
# 
