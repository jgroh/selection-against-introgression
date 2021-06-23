library(data.table)
# make bed file from recombination map

args <- commandArgs(trailingOnly = TRUE)
mapFile <- args[1]
outFile <- args[2]

map <- fread(mapFile)

map[, `:=` (cM_interval=pos_cM-shift(pos_cM),
            start=shift(pos_bp)-1,
            end=pos_bp-1), by = chr]

map[, r:= (cM_interval/100)/(end-start)]

map <- map[, .SD[-1], by=chr][,.(chr,start,end,marker,r)]

#boxplot(log10(map[,r]), ylab = "log10 per bp recomb rate")

write.table(map, file=outFile, quote=F, sep='\t', col.names=F, row.names=F)



