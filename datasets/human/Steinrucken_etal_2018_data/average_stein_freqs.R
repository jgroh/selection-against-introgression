library(data.table)

chromosome <- commandArgs(trailingOnly=T)[1]

# positions 
pos <- fread(paste0("March2018/CEU_lax_chr", chromosome, "/chr", chromosome, ".pos"))


# marginal posterior values at sites
d1 <- do.call(cbind, lapply(list.files(paste0("March2018/CEU_lax_chr", chromosome), pattern = ".filtered", full.names = T), fread))

# threshold posterior values
d2 <- ifelse(d1 > 0.42, 1, 0)

output <- data.table(chr = chromosome, pos = pos$V1, frq = rowMeans(d1), frq_thresh = rowMeans(d2))
outfile <- paste0("March2018/CEU_lax_chr", chromosome, "/chr", chromosome, "_frqs.txt")

fwrite(output, file = outfile, quote=F, sep = "\t", row.names=F, col.names=F)
