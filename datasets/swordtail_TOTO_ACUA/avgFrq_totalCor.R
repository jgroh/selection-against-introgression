# script combines scaffold data from each year and runs wavelets
library(data.table)
library(tools)
library(waveslim)

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]

scaffFiles <- dir(paste0("ACUA_",year),full.names=T)
scaffs <- basename(file_path_sans_ext(scaffFiles))

# each file has the same object name for the scaffold 
# so we load into separate environments
loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 

# combine data into one file for the year
gnom <- rbindlist(lapply(as.list(scaffFiles), 
                          function(x){loadFrom(x, "chrom1")}))
gnom[, pop_mean_minor_parent := 1 - pop_mean]
gnom[, ind_minor_parent := 1 - ind_frq]

# correlation on individuals
# ind_cor <- gnom[, cor(r, ind_minor_parent), by = .(ID)]
# setnames(ind_cor, "V1", "correlation")

# compute avg. minor parent frq and total correlation
avgFrq <- gnom[ID==ID[1], mean(pop_mean_minor_parent)]
totalCor <- gnom[ID==ID[1], cor(pop_mean_minor_parent, r)]

save(avgFrq, totalCor, file = paste0("ACUA_",year,"/avgFrq_totalCor.RData"))
