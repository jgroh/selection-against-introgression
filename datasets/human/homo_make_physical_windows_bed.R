library(data.table)

rmap <- fread("aau1043_datas3")
#rmap <- fread("recomb-hg38/genetic_map_GRCh38_merged.tab")
rmap[, Morgan := cM/100]

bed1kb <- rmap[Chr != "chrX", .(seq(min(End), max(End), by = 1e3) - 500,
                                  seq(min(End), max(End), by = 1e3) + 500), by = Chr]
fwrite(bed1kb, file = "", quote=F, sep="\t", col.names=F)
