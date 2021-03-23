library(data.table)
loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 

dlist <- as.list(dir(path = "sims/", pattern="*.RData", full.names = T))

allSimWV <- rbindlist(lapply(dlist, function(x){loadFrom(x, "simWV")}))

write.table(allSimWV, file = "all-sim-wv.txt", quote = F, row.names=F)