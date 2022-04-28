library(data.table)

rmap <- fread("recomb-hg38/genetic_map_GRCh38_merged.tab")
rmap[, Morgan := pos_cm/100]

bed1kb <- rmap[chrom != "chrX", .(seq(min(pos), max(pos), by = 1e3) - 500,
                                  seq(min(pos), max(pos), by = 1e3) + 500), by = chrom]
fwrite(bed1kb, file = "", quote=F, sep="\t", col.names=F)