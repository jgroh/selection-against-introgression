library(data.table)
library(waveslim)

args <- commandArgs(trailingOnly = T)
file <- args[1]
rep_gen <- args[2]
outfile <- args[3]

#file <- "/Users/jeff/workspace/selection-against-introgression/theory_and_simulations/results/equilibrium_test/test_haps.txt"

a <- fread(file, header=F)

# calculate exact variance among individuals

v <- var(rowMeans(a))

# calculate wavelet variances for population mean
d <- dwt(colMeans(a), "haar", n.levels = 10)
WV <- data.table(population = sapply(d[-1], function(x)sum(x^2)/1024), level = 1:10)

# calculate wavelet variances for individuals and append to table
for(i in 1:nrow(a)){
  d <- dwt(unlist(a[i,]), "haar", n.levels = 10)
  wv.ind <- data.table(x = sapply(d[-1], function(x)sum(x^2)/1024), level = 1:10)
  setnames(wv.ind, "x", paste0("ind_",i))
  WV <- merge(WV, wv.ind, all =T)
}

# average over individuals
WV[, meanInd  := rowMeans(WV[, !c("level", "population")])]

output <- WV[, c("level", "population", "meanInd")]

output[, var := v]
output[, rep_gen := rep_gen]

write.table(output, file = outfile, quote =F, row.names =F, col.names =T, sep = "\t")
