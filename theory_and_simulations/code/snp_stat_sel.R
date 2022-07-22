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
  
  setwd("/Users/Jeff/workspace/selection-against-introgression/theory_and_simulations/results/snp_stat_sel/")
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

# interpolate ancestry frq stat 
xout = seq(2^-7, 1, by = 2^-7)
frqs_informative <- frqs[, .SD[pos %in% sites[rep == .BY$rep, pos]], by = .(rep, gen)]

frqs_interp <- frqs_informative[, approx(x = Morgan, y = h_frq, xout=xout, rule=2), by = .(rep, gen)]
setnames(frqs_interp, c("x", "y"), c("Morgan", "hyb_frq"))


# interpolate recombination at same points
recmap[, Morgan := cumsum(rec*1e5)]
rec_interp <- recmap[, approx(x=Morgan, y = rec, xout=xout, rule=2)]
setnames(rec_interp, c("x", "y"), c("Morgan", "rec"))


all_frq_rec <- merge(frqs_interp, rec_interp)

#ggplot(all_frq_rec, aes(x = Morgan, y = hyb_frq)) + geom_point() + facet_wrap(~gen) 
#ggplot(all_frq_rec, aes(x = Morgan, y = rec)) + geom_point() + facet_wrap(~gen)

wavcor <- all_frq_rec[, cor_tbl(.SD, chromosome=NA, signals = c('hyb_frq','rec'), rm.boundary = T), by = .(rep, gen)]

#ggplot(wavcor, aes(x = level, y = cor)) + geom_point() + facet_wrap(~gen)


# ----- Wavelet Variance decomposition of frequency stat and recomb -----
wv_frq_rec <- all_frq_rec[, gnom_var_decomp(.SD, chromosome = NA, signals = c("hyb_frq", "rec")), by = .(rep, gen)]


output <- merge(wv_frq_rec, wavcor)

fwrite(output, file = "", quote = F, sep = "/t", row.names = F)

