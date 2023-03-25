library(data.table)
options(scipen=999)

rmap <- fread("Halldorsson_etal_2019_data/aau1043_datas3")
rmap[, Morgan := cM/100]

bed50kb <- rmap[, .(seq(min(End), max(End), by = 50e3) - 25e3,
                                  seq(min(End), max(End), by = 50e3) + 25e3), by = Chr]
fwrite(bed50kb, file = "", quote=F, sep="\t", col.names=F)
