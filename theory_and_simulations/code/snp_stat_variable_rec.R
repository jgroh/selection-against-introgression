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
  
  setwd("/Users/Jeff/workspace/selection-against-introgression/theory_and_simulations/results/snp_stat_variable_rec_sel/")
  n.sample <- 10
  #haps <- fread("replicate0_haps.txt", col.names = c("rep", "gen", "pos",  paste0("p0.", 1:20), paste0("p1.", 1:20), paste0("p2.", 1:20)))
  frqs <- fread("replicate0_frqs.txt", col.names = c("rep", "gen", "pos", "p0", "p1", "p2"))
  
  recLines <- readLines("../../variable_rec_map_slim.txt")
  recmap <- data.table(pos = as.numeric(unlist(strsplit(recLines[1], split = ','))), rec = as.numeric(unlist(strsplit(recLines[2], split = ','))))
  
} else {
  
  source("/Users/brogroh/gnomwav/R/multi_modwts.R")
  source("/Users/brogroh/gnomwav/R/variance_decomp.R")
  source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
  source("/Users/brogroh/gnomwav/R/theory.R")
  
  args <- commandArgs(trailingOnly = TRUE)
  frqs <- fread(args[1], col.names = c("rep", "gen", "pos", "p0", "p1", "p2"))
  recLines <- readLines(args[2])
  recmap <- data.table(pos = as.numeric(unlist(strsplit(recLines[1], split = ','))), rec = as.numeric(unlist(strsplit(recLines[2], split = ','))))
 
  sim <- args[3]
  #haps <- fread(args[3], col.names = c("rep", "gen", "pos",  paste0("p0.", 1:n.sample), paste0("p1.", 1:n.sample), paste0("p2.", 1:n.sample)))
}

# ----- format data ------
#haps <- melt(haps, id.vars = c('rep', 'gen', 'pos'), value.name = "allele")
#haps[, c("pop", "id") := tstrsplit(variable, ".", fixed=T)][, variable := NULL]

#haps[, Morgan := pos/1e8]
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


# ===== interpolate on physical scale =====
# interpolate ancestry frq stat 
xoutMorgan = seq(2^-7, 1, by = 2^-7)
xoutBP <- seq(1, 1e8, by = 500000)

frqs_informative <- frqs[, .SD[pos %in% sites[rep == .BY$rep, pos]], by = .(rep, gen)]
#ggplot(frqs_informative, aes(x = Morgan, y = h_frq)) + geom_point() + facet_wrap(~gen)


frqs_interp_P <- frqs_informative[, approx(x = pos, y = h_frq, xout=xoutBP, rule = 2), by = .(rep, gen)]
frqs_interp_G <- frqs_informative[, approx(x = Morgan, y = h_frq, xout=xoutMorgan, rule=2), by = .(rep, gen)]

setnames(frqs_interp_P, c("x", "y"), c("pos", "hyb_frq"))
setnames(frqs_interp_G, c("x", "y"), c("Morgan", "hyb_frq"))
#ggplot(frqs_interp_G, aes(x = Morgan, y = hyb_frq)) + geom_point() + facet_wrap(~gen)
#ggplot(frqs_interp_P[gen %in% c(3, 10, 100, 1000)], aes(x = pos, y = hyb_frq)) + geom_point() + facet_wrap(~gen)

# interpolate recombination at same points
recmap[, Morgan := cumsum(rec*1e5)]
rec_interp_G <- recmap[, approx(x=Morgan, y = rec, xout=xoutMorgan, rule=2)]
setnames(rec_interp_G, c("x", "y"), c("Morgan", "rec"))

rec_interp_P <- recmap[, approx(x=pos, y = rec, xout=xoutBP, rule=2)]
setnames(rec_interp_P, c("x", "y"), c("pos", "rec"))

all_frq_rec_P <- merge(rec_interp_P, frqs_interp_P)
all_frq_rec_G <- merge(rec_interp_G, frqs_interp_G)

#ggplot(all_frq_rec, aes(x = Morgan, y = hyb_frq)) + geom_point() + facet_wrap(~gen) 
#ggplot(all_frq_rec, aes(x = Morgan, y = rec)) + geom_point() + facet_wrap(~gen)

wavcor_P <- all_frq_rec_P[, cor_tbl(.SD, chromosome=NA, signals = c('hyb_frq','rec'), rm.boundary = T), by = .(rep, gen)]
wavcor_P[, units := 'physical']
wavcor_G <-  all_frq_rec_G[, cor_tbl(.SD, chromosome=NA, signals = c('hyb_frq','rec'), rm.boundary = T), by = .(rep, gen)]
wavcor_G[, units := 'genetic']

#ggplot(wavcor, aes(x = level, y = cor)) + geom_point() + facet_wrap(~gen)


# ----- Wavelet Variance decomposition of frequency stat and recomb -----
wv_frq_rec_G <- all_frq_rec_G[, gnom_var_decomp(.SD, chromosome = NA, signals = c("hyb_frq", "rec")), by = .(rep, gen)]
wv_frq_rec_G[, units := 'genetic']
wv_frq_rec_P <- all_frq_rec_P[, gnom_var_decomp(.SD, chromosome = NA, signals = c("hyb_frq", "rec")), by = .(rep, gen)]
wv_frq_rec_P[, units := 'physical']

output <- merge(rbind(wv_frq_rec_G, wv_frq_rec_P), rbind(wavcor_G, wavcor_P))
output[, sim := sim]

fwrite(output, file = "", quote = F, sep = "\t", row.names = F)

