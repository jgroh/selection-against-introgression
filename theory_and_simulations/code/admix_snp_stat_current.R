library(data.table)
library(tidyverse)
library(cubature)

# ----- read and format simulation data

if(Sys.getenv("RSTUDIO") == "1"){
  
  source("~/workspace/gnomwav/R/multi_modwts.R")
  source("~/workspace/gnomwav/R/variance_decomp.R")
  source("~/workspace/gnomwav/R/theory.R")
  
  setwd("/Users/Jeff/workspace/selection-against-introgression/theory_and_simulations/results/admix_snp_stat/")
  haps <- fread("replicate1_haps.txt", col.names = c("rep", "gen", "pos",  paste0("p0.", 1:8), paste0("p1.", 1:8), paste0("p2.", 1:8)))
  frqs <- fread("replicate1_frqs.txt", col.names = c("rep", "gen", "pos", "p0", "p1", "p2"))
  
} else {
  
  source("/Users/brogroh/gnomwav/R/multi_modwts.R")
  source("/Users/brogroh/gnomwav/R/variance_decomp.R")
  source("/Users/brogroh/gnomwav/R/theory.R")
  
  args <- commandArgs(trailingOnly = TRUE)
  haps <- fread(args[1], col.names = c("rep", "gen", "pos",  paste0("p0.", 1:8), paste0("p1.", 1:8), paste0("p2.", 1:8)))
  frqs <- fread(args[2], col.names = c("rep", "gen", "pos", "p0", "p1", "p2"))
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
frqs[, h_frq := (p2-p0)/(p1-p0)]

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
setnames(frqs_interp, c("x", "y"), c("Morgan", "h_frq"))

# ----- Wavelet Variance decomposition of frequency stat -----
wv_frq <- frqs_interp[, gnom_var_decomp(.SD, chromosome = NA, signals = "h_frq"), by = .(rep, gen)]
wv_frq <- wv[grepl("d", level)]


#thry_frq <- wavelet_variance_equilbrium(n.pop = 20000, n.sample = 20000, unit.scale = 2^-7, gen = c(10,100,500,1000),level = 1:7, alpha = 0.5)
# unclear currently how theory for frq of known admixture relates to mean of snp stat

#ggplot(wv[gen !=0], aes(x = level, y = variance.h_frq)) + geom_point() + 
#  geom_line(aes(group = rep)) + facet_wrap(~gen) + 
#  scale_x_discrete(breaks = paste0("d",1:7), labels = as.character(-7:-1)) + 
#  labs(x = expression(Scale: log[2] (Morgan)), 
#       y = "Variance") + 
#  geom_line(data=thry_frq, aes(x = level, y = variance), color = "red") + 
#  theme(aspect.ratio=1)



# ===== Wavelet Variance Decomp on individual haps =====

#haps[allele>1] # what to do about these? maybe just ignore non-biallelic sites for now
haps_informative <- haps[allele != 2, .SD[pos %in% sites[rep == .BY$rep, pos]], by = .(rep, gen)]

# combine haps and frqs to get informative sites
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

#ggplot(wv_haps[pop != "p2" & grepl("d", level, fixed=T), 
#               .(v = mean(variance.h_hap)), by = .(gen, level)], 
#       aes(x = level, y = v)) + geom_point() +
#  geom_line(aes(group = gen, color = as.factor(gen)))

# *** option to analyze WV on parental chromosomes only from gen zero or from generation of analysis
# only gen 0
#prnt_wv_haps <- wv_haps[pop!="p2" & gen ==0, .(prnt_var = mean(variance.h_hap)), by = level]

# contemporary gens
prnt_wv_haps <- wv_haps[pop!="p2", .(prnt_var = mean(variance.h_hap)), by = .(level, gen)]

hyb_wv_haps <- wv_haps[pop == "p2", .(hyb_var = mean(variance.h_hap)), by = .(gen, level)]

wv_haps2 <- merge(prnt_wv_haps, hyb_wv_haps, by = c("level","gen"))

# combine frq and hap data for output
allWV <- merge(wv_haps2, wv_frq)





#wvtheory <- wavelet_variance_equilbrium(n.pop=20000, n.sample = 1, unit.scale = 2^-7, level = 1:7, alpha = 0.5, gen = c(10,100,500,1000))

# compare wv theory with actual wv on snp stat
# ggplot(hyb_wv_haps[gen != 0], aes(x = level, y = var)) + geom_point() +
#   geom_line(group = 1) + facet_wrap(~gen) +
#   geom_line(data = wvtheory,
#             aes(x = level, y = variance), 
#             color  = 'red') + 
#   theme(aspect.ratio=1) +
#   scale_x_discrete(breaks = paste0("d",1:7), labels = as.character(-7:-1)) + 
#   labs(x = expression(Scale: log[2] (Morgan)), 
#        y = "Variance")


#ggplot(wv_haps2[gen != 0], aes(x = level, y = adj_var)) + geom_point() + 
#  geom_line(group = 1) +
#  facet_wrap(~gen) + 
#  geom_line(data = wvtheory, 
#            aes(x = level, y = variance), color = 'red') +
#  theme(aspect.ratio=1) + 
#  scale_x_discrete(breaks = paste0("d",1:7), labels = as.character(-7:-1)) + 
#  labs(x = expression(Scale: log[2] (Morgan)), 
#       y = "Variance", color = NULL)




