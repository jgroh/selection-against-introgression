library(data.table)
options(scipen=999)

bed1kb <- NULL

for(i in c(1:22, "X")){
  chr <- i
  file1 <- paste0("Sankararaman_etal_2014_data/chr-", chr, ".thresh-90.length-0.00.gz")
  file2 <- paste0("Steinrucken_etal_2018/March2018/CEU_lax_chr", chr, "/chr", chr, "_frqs.txt")
  
  sank_calls <- suppressWarnings(fread(file1))
  stein_calls <- fread(file2, col.names = c("chr", "pos", "freq"))
  
  mn <-  max(min(stein_calls$pos), min(sank_calls$V4))
  mx <- min(max(stein_calls$pos), max(sank_calls$V4))
  
  bed <- data.table(chr = paste0('chr', chr), start = seq(mn, mx, by = 1e3) - 500, end = seq(mn, mx, by = 1e3) + 500)
  bed1kb <- rbind(bed1kb, bed)
}

fwrite(bed1kb, file = "", quote=F, sep="\t", col.names=F)
