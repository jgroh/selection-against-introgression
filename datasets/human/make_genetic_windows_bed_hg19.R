library(data.table)
options(scipen=999)

bedM <- NULL

for(i in c(1:22, "X")){
  
  chr <- i
  s <- suppressWarnings(fread(paste0("Sankararaman_etal_2014_data/summaries.release/CEU.hapmap/summaries/chr-", chr, ".thresh-90.length-0.00.gz")))
  
  xoutM <- seq(min(s$V3)/100, max(s$V3/100), by = 2^-16)
  
  bed <- data.table(chr = paste0("chr", chr),
                    start = round(approx(x = s$V3/100, y = s$V4, xout = xoutM - 2^-17)$y),
                    end = round(approx(x = s$V3/100, y = s$V4, xout = xoutM + 2^-17)$y), 
                    Morgan = xoutM)
  bedM <- rbind(bedM, bed)
}

bedM <- bedM[!is.na(start) & !is.na(end)]
fwrite(bedM[, .(chr, start, end, Morgan)], file = "", quote=F, sep="\t", col.names=F)
