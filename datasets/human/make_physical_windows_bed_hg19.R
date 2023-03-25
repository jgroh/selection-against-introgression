library(data.table)
options(scipen=999)

bed50kb <- NULL

for(i in c(1:22, "X")){
  chr <- i
  file1 <- paste0("Sankararaman_etal_2014_data/summaries.release/CEU.hapmap/summaries/chr-", chr, ".thresh-90.length-0.00.gz")
  sank_calls <- suppressWarnings(fread(file1))
  
  mn <-   min(sank_calls$V4)
  mx <-  max(sank_calls$V4)
  
  bed <- data.table(chr = paste0('chr', chr), start = seq(mn, mx, by = 50e3) - 25e3, end = seq(mn, mx, by = 50e3) + 25e3)
  bed50kb <- rbind(bed50kb, bed)
  bed50kb <- bed50kb[start > 0]
}

fwrite(bed50kb, file = "", quote=F, sep="\t", col.names=F)
