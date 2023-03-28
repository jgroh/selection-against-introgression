library(data.table)
library(gnomwav)

if(interactive()){
  frq <- fread('results/human_wgsim_neutral_equilibrium_no_interp/replicate9_frqs.txt', col.names = c("rep", "gen", "frq"))
  chrlengths <- fread("hg38_chr_lengths.txt", col.names = 'Morgans')

} else{
  args <- commandArgs(trailingOnly = T)
  frq <- fread(args[1], col.names = c("rep", "gen", "frq"))
  chrlengths <- fread(args[2], col.names = 'Morgans')
}

# subtract gen to correspond to thry
frq[, gen := gen - 2]

# get chromosome IDS of loci
nloci <- chrlengths[, round(Morgans/2^-15)]

chrID = rep(1, nloci[1]+1)
for(i in 2:22){
  chrID = c(chrID,  rep(i, nloci[i]+1))
}

frq[, chr := chrID, by = gen]

# run wavelet variance
wv <- frq[, gnom_var_decomp(.SD, chromosome = 'chr', signals = 'frq', rm.boundary = F, avg.over.chroms = T), by = .(rep, gen)]

fwrite(wv, file = gsub("_frqs.txt.gz", "_wavelet_results.txt", args[1]), col.names = T, quote = F, row.names = F, sep = "\t")
