library(data.table)
options(scipen=999)

rmap <- fread("Halldorsson_etal_2019_data/aau1043_datas3")
rmap[, Morgan := cM/100]

bedM <- rmap[, .(Morgan = seq(min(Morgan), max(Morgan), by = 2^-16)), by = Chr]

bedM[, wstart := Morgan - 0.5*2^-16]
bedM[, wend := wstart + 2^-16]

bedM[, bp_start := round(approx(x = rmap[Chr == .BY, Morgan], 
                                y = rmap[Chr == .BY, End], 
                                xout = bedM[Chr == .BY, wstart])$y), by = Chr]

bedM[, bp_end := round(approx(x = rmap[Chr == .BY, Morgan], 
                              y = rmap[Chr == .BY, End], 
                              xout = bedM[Chr == .BY, wend])$y), by = Chr]

bedM <- bedM[!is.na(bp_start) & !is.na(bp_end)]

fwrite(bedM[, .(Chr, bp_start, bp_end, Morgan)], file = "", quote=F, sep="\t", col.names=F)
