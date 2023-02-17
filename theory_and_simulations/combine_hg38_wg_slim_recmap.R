library(data.table)
#library(ggplot2)

map <- rbindlist(lapply(list.files(path = "hg38_slim_recmap_chrs/", full.names=T), fread))
setnames(map, c("chr", "pos", "cM", "Morgan_dist", "rbar_i"))

#ggplot(map[seq(1, nrow(map), by = 100)], aes(x = pos, y = rbar_i)) + geom_point() + facet_wrap(~chr)

map[,chr := as.numeric(gsub("chr", "", chr))]
setkey(map, chr, pos)

fwrite(map[chr==22, .(Morgan_dist)], 
       file = "hg38_chr22_slim_recmap.txt", quote = F, sep = "\t", col.names = F)

fwrite(map[chr==22, .(chr, pos, Morgan = cM/100, Morgan_dist, rbar_i)], 
       file = "hg38_chr22_slim_recmap_verbose.txt.gz", quote = F, sep = "\t", col.names = T)

fwrite(map[, .(Morgan_dist)], 
       file = "hg38_wg_slim_recmap.txt", quote = F, sep = "\t", col.names = F)

fwrite(map[, .(chr, pos, Morgan = cM/100, Morgan_dist, rbar_i)], 
       file = "hg38_wg_slim_recmap_verbose.txt.gz", quote = F, sep = "\t", col.names = T)

fwrite(map[, .(len = max(cM)/100), by = chr][, .(len)], 
       file = "hg38_chr_lengths.txt", quote = F, sep = "\t", col.names = F)
                                                        
