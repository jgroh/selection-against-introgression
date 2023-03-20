# script combines scaffold data from each year and runs wavelets
library(data.table)
library(tools)
library(waveslim)

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]

scaffFiles <- dir(paste0("ACUA_",year), pattern = "ScyDAA6*",full.names=T)
scaffs <- basename(file_path_sans_ext(scaffFiles))

# each file has the same object name for the scaffold 
# so we load into separate environments
loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 

# combine data into one file for the year
gnom <- rbindlist(lapply(as.list(scaffFiles), 
                          function(x){loadFrom(x, "chromAnc50kb")}))
gnom[, minPrntAnc := 1 - meanFreq]

gnom[, mean(minPrntAnc)]

# compute avg. minor parent frq and total correlation
avgFrq <- gnom[ID==ID[1], mean(minPrntAnc)]

save(avgFrq, file = paste0("ACUA_",year,"/avgFrq.txt"))
