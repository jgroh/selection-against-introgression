library(data.table)
options(scipen=999)

rmap <- fread("Halldorsson_etal_2019_data/aau1043_datas3")
rmap[, Morgan := cM/100]

bed1kb <- rmap[, .(seq(min(End), max(End), by = 1e3) - 500,
                                  seq(min(End), max(End), by = 1e3) + 500), by = Chr]
fwrite(bed1kb, file = "", quote=F, sep="\t", col.names=F)
