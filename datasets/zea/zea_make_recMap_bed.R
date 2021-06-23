library(data.table)
# make bed file from recombination map

args <- commandArgs(trailingOnly = TRUE)
mapFile <- args[1]
outFile <- args[2]

map <- fread(mapFile)

map[, `:=` (cM_interval=pos_cM-shift(pos_cM),
            start=shift(pos_bp)-1,
            end=(pos_bp-1)), by = chr]

# add 1 to end of last interval for zero-based indexing of bed file
idx <- map[, .(idx = .I[end==max(end)]), by=chr]$idx
map[idx, end := end+1]

# get average per bp recombination rate for interval
map[, r:= (cM_interval/100)/(end-start)]

# remove first value and select column order
map <- map[, .SD[-1], by=chr][, .(chr,start,end,marker,r)]

#boxplot(log10(map[,r]), ylab = "log10 per bp recomb rate")

write.table(map, file=outFile, quote=F, sep='\t', col.names=F, row.names=F)



