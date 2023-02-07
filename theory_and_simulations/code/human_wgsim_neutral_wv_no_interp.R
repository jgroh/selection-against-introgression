library(data.table)
#library(ggplot2)
library(cubature)

if(interactive()){
  source("/Users/jeff/workspace/gnomwav/R/correlation_decomp.R")
  source("/Users/jeff/workspace/gnomwav/R/multi_modwts.R")
  source("/Users/jeff/workspace/gnomwav/R/variance_decomp.R")
  frq <- fread('results/human_wgsim_const_rec/replicate9_frqs.txt', col.names = c("rep", "gen", "frq"))
  chrlengths <- fread("hg38_chr_lengths.txt", col.names = 'Morgans')
} else{
  source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
  source("/Users/brogroh/gnomwav/R/multi_modwts.R")
  source("/Users/brogroh/gnomwav/R/variance_decomp.R")
  
  args <- commandArgs(trailingOnly = T)
  frq <- fread(args[1], col.names = c("rep", "gen", "frq"))
  chrlengths <- fread(args[2], col.names = 'Morgans')
}

# subtract gen to correspond to thry
frq[, gen := gen - 2]

# get chromosome IDS of loci
nloci <- chrlengths[, round(Morgans/2^-16)]

chrID = rep(1, nloci[1]+1)
for(i in 2:22){
  chrID = c(chrID,  rep(i, nloci[i]+1))
}

frq[, chr := chrID, by = gen]

# run wavelet variance
wv <- frq[, gnom_var_decomp(.SD, chromosome = 'chr', signals = 'frq', rm.boundary = F, avg.over.chroms = T), by = .(rep, gen)]

fwrite(wv, file = gsub("_frqs.txt", "_wavelet_results.txt", args[1]), col.names = T, quote = F, row.names = F, sep = "\t")
