library(data.table)
library(waveslim)
library(magrittr)

# read selection coefficient, generation, replicate
args <- commandArgs(trailingOnly = TRUE)
sim <- args[1]

params <- unlist(strsplit(sim, "_"))
s <- gsub("s","", params[1])
l.s <- gsub("replicate","", params[2])
gen <- gsub("replicate","", params[3])

# read genotypes
# snps are rows, individuals columns
f <- fread(file = paste0("sims/",sim,".txt"))

# selected allele is coded as '2' so change to 1
f <- f[, lapply(.SD, function(x){x[x>1] <- 1; x})]

# reformat to long
f[,"position" := 1:1024]
f <- melt(f, id.vars = "position", variable.name = "id", 
     value.name = "allele")

# look at allele frequency
# f[, mean(allele), by = position] %>% 
#   ggplot(aes(x = position, y = V1)) + geom_point()
 
# look at haplotypes
#f[id == unique(id)[100]] %>% 
#  ggplot(aes(x = position, y= allele)) + geom_point() +
#  geom_line()

# wavelet coefficients for individuals
mwtInd <- f[,brick.wall(wf="haar",x=modwt(allele, "haar", n.levels = 10)), by = id]
mwtInd[, position := 1:1024, by = id]
mwtAvgInd <- mwtInd[, mean(d9), by = position] 
setnames(mwtAvgInd, "V1", "d9")

# wavelet coefficidnts for population mean
frq <- f[, mean(allele), by=position]
setnames(frq, "V1", "frequency")
mwtMean <- frq[,brick.wall(wf="haar",x=modwt(frequency, "haar", n.levels = 10))]
mwtMean <- data.table(mwtMean$d9)
mwtMean[, position := 1:1024]
setnames(mwtMean, "V1", "d9")

# plot scale 9 wavelet coefficients
#mwtInd[id %in% unique(id)[1:20]] %>% ggplot(aes(x = position, y = d9)) +
#  geom_point() +
#  geom_line(aes(group = id, color = id)) +
#  theme(legend.position = "none") +
#  geom_line(data=mwtAvgInd, size = 1.5, alpha = 1) +
#  labs(y = "Scale 9 coefficients") + 
#  theme_classic()


# wvCustom <- function(x){
#   return(sum(x^2, na.rm=T)/length(x[!is.na(x)]))
# }
# wvIndTmp <- mwtInd[, lapply(.SD, wvCustom), by = id]
# plot(as.numeric(wvIndTmp[, lapply(.SD,mean), .SDcols = paste0("d", 1:10)]))

# wavelet variance by individual haplotype
wvInd <- f[,wave.variance(brick.wall(modwt(allele, "haar", n.levels = 10),"haar")), by = id]
wvInd[, scale:= 1:11, by= id]
wvInd <- wvInd[scale != 11] # remove large scale smooth 
wvInd[, signal := "mean_ind"]

# average over haplotypes
wvAvgInd <- wvInd[, mean(wavevar), by = scale]
setnames(wvAvgInd, "V1", "wavevar")
wvAvgInd[,signal := "mean_ind"][, s := s][, l.s := l.s][, gen := gen][]

# wavelet variance for population mean
wvPopMean <- as.data.table(frq[,wave.variance(brick.wall(modwt(frequency, "haar", n.levels = 10),"haar"))])
wvPopMean[, scale:= 1:11]
wvPopMean <- wvPopMean[scale != 11] # remove large scale smooth 
wvPopMean[, lower:=NULL]
wvPopMean[, upper:=NULL]
wvPopMean[,signal := "pop_mean"]
wvPopMean[, s := s][, l.s := l.s][, gen := gen][]

# output 
simWV <- merge(wvPopMean, wvAvgInd, all = T)

save(simWV, file = paste0("sims/",sim,".RData"))

# plot wavelet spectrum
#simWV[signal == "pop_mean"] %>% 
#  ggplot(aes(x = scale, y = wavevar, group = signal)) + 
#  geom_point() + geom_line() + theme_classic()
