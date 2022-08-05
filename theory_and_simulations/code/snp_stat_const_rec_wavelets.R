library(data.table)
library(cubature)
library(ggplot2)
library(magrittr)


# ----- read and format simulation data

if(Sys.getenv("RSTUDIO") == "1"){
  
  source("~/workspace/gnomwav/R/multi_modwts.R")
  source("~/workspace/gnomwav/R/variance_decomp.R")
  source("~/workspace/gnomwav/R/theory.R")
  source("~/workspace/gnomwav/R/correlation_decomp.R")
  
  setwd("/Users/Jeff/workspace/selection-against-introgression/theory_and_simulations/results/admix_snp_stat/")
  n.sample <- 20
  haps <- fread("replicate0_haps.txt", col.names = c("rep", "gen", "pos",  paste0("p0.", 1:20), paste0("p1.", 1:20), paste0("p2.", 1:20)))
  frqs <- fread("replicate0_frqs.txt", col.names = c("rep", "gen", "pos", "p0", "p1", "p2"))
  
} else {
  
  source("/Users/brogroh/gnomwav/R/multi_modwts.R")
  source("/Users/brogroh/gnomwav/R/variance_decomp.R")
  source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
  source("/Users/brogroh/gnomwav/R/theory.R")
  
  args <- commandArgs(trailingOnly = TRUE)
  n.sample <- as.numeric(args[1])
  haps <- fread(args[2], col.names = c("rep", "gen", "pos",  paste0("p0.", 1:n.sample), paste0("p1.", 1:n.sample), paste0("p2.", 1:n.sample)))
  frqs <- fread(args[3], col.names = c("rep", "gen", "pos", "p0", "p1", "p2"))
}

# ----- format data ------
haps <- melt(haps, id.vars = c('rep', 'gen', 'pos'), value.name = "allele")
haps[, c("pop", "id") := tstrsplit(variable, ".", fixed=T)][, variable := NULL]

haps[, Morgan := pos/1e8]
frqs[, Morgan := pos/1e8]

# confirm hybrid allele frequencies are intermediate
#frqs[abs(p0 - p1) > 0.8]  %>% ggplot(aes(x = pos)) + geom_point(aes(y = p0, color = 'black')) + 
 # geom_point(aes(y=p1, color = 'red')) + geom_point(aes(y=p2, color = 'green')) + 
 # facet_wrap(~gen)


# ----- calculate ancestry statistic -----
frqs[, h_frq := (p2-p0)/(p1-p0)] # approximates the indicator for p1 ancestry

# what is the density of snps
#frqs[abs(p0 - p1) >= 0.5 & gen == 0, mean(Morgan-shift(Morgan), na.rm=T), by = .(rep,gen)] 
#frqs[abs(p0 - p1) >= 0.75 & gen == 0, mean(Morgan-shift(Morgan), na.rm=T), by= .(rep,gen)] 
#frqs[abs(p0 - p1) >= 0.9 & gen == 0, mean(Morgan-shift(Morgan), na.rm=T), by= .(rep,gen)] 
#frqs[abs(p0 - p1) == 1 & gen == 0, mean(Morgan-shift(Morgan), na.rm=T), by= .(rep,gen)] 


# are informative SNPs informative through time
# sites0 <- frqs[abs(p0 - p1) >= .9 & gen == 0, pos]
# sites1000 <- frqs[abs(p0 - p1) >= .9 & gen == 1000, pos]
# 
# frqs[pos %in% sites1000, abs(p0-p1), by = .(gen,pos,rep)] %>%
#   ggplot(aes(x = log10(gen), y = V1)) + 
#   geom_point() + geom_line(aes(group = pos))


## *** Option to either restrict to sites that are informative in gen 0 or sites informative at time of analysis
# ascertain informative SNPs in gen 0
sites <- frqs[abs(p0 - p1) >= 0.75 & gen == 0, pos, by = rep] 

# interpolate ancestry frq stat 
xout = seq(2^-7, 1, by = 2^-7)
frqs_informative <- frqs[, .SD[pos %in% sites[rep == .BY$rep, pos]], by = .(rep, gen)]

frqs_interp <- frqs_informative[, approx(x = Morgan, y = h_frq, xout=xout, rule=2), by = .(rep, gen)]
setnames(frqs_interp, c("x", "y"), c("Morgan", "hyb_frq"))

# ----- Wavelet Variance decomposition of frequency stat -----
wv_frq <- frqs_interp[, gnom_var_decomp(.SD, chromosome = NA, signals = "hyb_frq"), by = .(rep, gen)]
wv_frq <- wv_frq[grepl("d", level)]


#thry_frq <- wavelet_variance_equilbrium(n.pop = 20000, n.sample = 20000, unit.scale = 2^-7, gen = c(10,100,500,1000),level = 1:7, alpha = 0.5)
# unclear currently how theory for frq of known admixture relates to mean of snp stat

# ggplot(wv_frq[gen !=0], aes(x = level, y = variance.h_frq)) + geom_point() +
#  geom_line(aes(group = rep)) + facet_wrap(~gen, scales = "free_y") +
#  scale_x_discrete(breaks = paste0("d",1:7), labels = as.character(-7:-1)) +
#  labs(x = expression(Scale: log[2] (Morgan)),
#       y = "Variance") +
#  geom_line(data=thry_frq, aes(x = level, y = variance), color = "red") +
#  theme(aspect.ratio=1)


# ===== Wavelet Variance Decomp on individual haps =====

#haps[allele>1] # what to do about these? maybe just ignore non-biallelic sites for now
haps_informative <- haps[allele != 2, .SD[pos %in% sites[rep == .BY$rep, pos]], by = .(rep, gen)]

# combine haps and frqs to get informative sites (sites ascertained at time of admixture)
allsitedata <- merge(haps_informative, 
                frqs_informative, 
                by = c("rep", "gen", "pos", "Morgan"))

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

# all 
# ggplot(wv_haps[gen != 0 & grepl("d", level, fixed=T), .(v = mean(variance.h_hap)), by = .(gen, pop, level)],
#        aes(x = level, y = v)) +
#   geom_point(aes(group = pop)) +
#   geom_line(aes(color = pop, group = pop)) +
#   facet_wrap(~gen)+
#   scale_x_discrete(breaks = paste0("d",1:7), labels = as.character(-7:-1)) +
#   labs(x = expression(Scale: log[2] (Morgan)),
#        y = "Variance") +
#   theme(aspect.ratio=1)

# ggplot(wv_haps[pop != "p2" & grepl("d", level, fixed=T),
#               .(v = mean(variance.h_hap)), by = .(gen, level)],
#       aes(x = level, y = v)) + geom_point() +
#  geom_line(aes(group = gen, color = as.factor(gen)))

# *** option to analyze WV on parental chromosomes only from gen zero or from generation of analysis
# only gen 0
#prnt_wv_haps <- wv_haps[pop!="p2" & gen ==0, .(prnt_var = mean(variance.h_hap)), by = level]


# reformat so that single chrom variance for parental chromosomes and admixed chromosomes are in separate columns
prnt_wv_haps <- wv_haps[pop!="p2", .(variance.prnt_haps = mean(variance.h_hap)), by = .(level, gen)]

hyb_wv_haps <- wv_haps[pop == "p2", .(variance.hyb_haps = mean(variance.h_hap)), by = .(gen, level)]

wv_haps2 <- merge(prnt_wv_haps, hyb_wv_haps, by = c("level","gen"))


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
haps_interp[, pos := paste0("d", pos)]
haps_interp3 <- dcast(haps_interp, rep + gen + pop + id ~ pos, value.var = 'h_hap')

cmb34 <- combn(n.sample, 2)
cmb5 <- expand.grid(1:n.sample, 1:n.sample)

L <- 128
all_hh <- data.table()

for (g in 0){ # look at the product for parental chromosomes in gen zero, should not make a huge difference across gens
  
  # these vectors will hold the mean of hl * hl' for the same chromosome with entries corresponding to distance between l and l'
  hh_00_same <- vector()
  hh_11_same <- vector()
  
  # these vectors will hold the mean of hl * hl' for different chromosomes with entries corresponding to distance between l and l'
  hh_00_diff <- vector()
  hh_11_diff <- vector()
  hh_01_diff <- vector()
  
  for(dist in 1:(L-1)){ # loop over possible distances between l and l'
    
    # this vector will hold values of the product and we'll take it's mean at the end
    vals_00_dist <- vector()
    vals_11_dist <- vector()
    vals_01_dist <- vector()
    
    # we could average over all starting positions but this seemed to take forever. plus we'll be averaging over replicate sims
    l1 <- 1
    l2 <- l1 + dist
    cols <- paste0('d', c(l1,l2))
    
    hh_00_same[dist] <- haps_interp3[pop=="p0" & gen == g, mean(get(cols[1])*get(cols[2]))]
    hh_11_same[dist] <- haps_interp3[pop=="p1" & gen == g, mean(get(cols[1])*get(cols[2]))]
    
    np <- ncol(cmb34)
    for(h in 1:np){
      # ids of haplotypes
      i <- cmb34[1, h]
      j <- cmb34[2, h]
      
      vals_00_dist  <- c(vals_00_dist, haps_interp3[id == i & gen == g & pop == "p0", get(cols[1])]*haps_interp3[id == j & gen == g & pop == "p0", get(cols[2])])
      vals_11_dist <- c(vals_11_dist, haps_interp3[id == i & gen == g & pop == "p1", get(cols[1])]*haps_interp3[id == j & gen == g & pop == "p1", get(cols[2])])
    }
    
    hh_00_diff[dist] <- mean(vals_00_dist)
    hh_11_diff[dist] <- mean(vals_11_dist)
    
    # for term hh_01_diff, loop over combinations one drawn from each pop
    np2 <- nrow(cmb5)
    for(h in 1:np2){
      i <- cmb5[h, 1]
      j <- cmb5[h, 2]
      
      vals_01_dist <- c(vals_01_dist, haps_interp3[id == i & gen == g & pop == "p0", get(cols[1])]*haps_interp3[id == j & gen == g & pop == "p1", get(cols[2])])
    }
    
    hh_01_diff[dist] <- mean(vals_01_dist)
    
  }
  
  hh <- data.table(l1l2_dist = 1:(L-1), gen=g, rep= sites[1, rep], hh_00_same = hh_00_same, hh_11_same=hh_11_same, hh_00_diff=hh_00_diff, hh_11_diff=hh_11_diff, hh_01_diff = hh_01_diff) 
  all_hh <- rbind(all_hh, hh)
}

#ggplot(melt(all_hh, id.vars = c("l1l2_dist", "gen", "rep")), aes(x = l1l2_dist, y = value, color = variable)) + geom_point()
save(allWV, all_hh, file = paste0("results/admix_snp_stat/replicate",allWV[1,rep], "_wv_results.RData"))

