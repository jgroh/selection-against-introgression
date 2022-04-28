library(data.table)

rmap <- fread("recomb-hg38/genetic_map_GRCh38_merged.tab")
rmap[, Morgan := pos_cm/100]

bedM <- rmap[chrom != "chrX", .(Morgan = seq(min(Morgan), max(Morgan), by = 2^-16)), by = chrom]

bedM[, wstart := Morgan - 0.5*2^-16]
bedM[, wend := wstart + 2^-16]

bedM[, bp_start := round(approx(x = rmap[chrom == .BY, Morgan], 
                                y = rmap[chrom == .BY, pos], 
                                xout = bedM[chrom == .BY, wstart])$y), by = chrom]
bedM[, bp_end := round(approx(x = rmap[chrom == .BY, Morgan], 
                              y = rmap[chrom == .BY, pos], 
                              xout = bedM[chrom == .BY, wend])$y), by = chrom]
bedM <- bedM[!is.na(bp_start) & !is.na(bp_end)]

fwrite(bedM[, .(chrom, bp_start, bp_end, Morgan)], file = "", quote=F, sep="\t", col.names=F)
