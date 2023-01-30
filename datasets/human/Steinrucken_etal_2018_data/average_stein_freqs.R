library(data.table)

chromosome <- commandArgs(trailingOnly=T)[1]

d <- do.call(cbind, lapply(list.files(paste0("March2018/CEU_lax_chr", chromosome), pattern = ".filtered", full.names = T), fread))

fwrite(data.table(rowMeans(d)), file = paste0("March2018/CEU_lax_chr", chromosome, "/chr", chromosome, "_frq.txt"), quote=F, row.names=F, col.names=F)
