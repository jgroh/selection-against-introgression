library(data.table)
library(waveslim)

# read selection coefficient, generation, replicate
args <- commandArgs(trailingOnly = TRUE)
sim <- args[1]

params <- unlist(strsplit(sim, "_"))
s <- gsub("s","", params[1])
l.s <- gsub("replicate","", params[2])
gen <- gsub(".txt","",gsub("replicate","", params[3]))

# read genotypes
# snps are rows, individuals columns
f <- fread(file = paste0("sims/",sim))

# selected allele is coded as '2' so change to 1
f <- f[, lapply(.SD, function(x){x[x>1] <- 1; x})]

# reformat to long
f[,"position" := 1:1024]
f <- melt(f, id.vars = "position", variable.name = "id", 
     value.name = "allele")

# wavelet variance by individual
mwt <- f[,wave.variance(brick.wall(modwt(allele, "haar", n.levels = 10),"haar")), by = id]
mwt[, scale:= 1:11, by= id]
mwt <- mwt[scale != 11]

avg.wv <- mwt[, mean(wavevar), by = scale]
setnames(avg.wv, "V1", "wavevar")
avg.wv[, s := s]
avg.wv[, l.s := l.s]
avg.wv[, gen := gen]

save("avg.wv", file = paste0("sims/",sim,".RData"))
