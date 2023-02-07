library(data.table)

snps <- fread("local_ancestry/results/thinnedSNPs/HILO_MAIZE55/K2/whole_genome.var.sites")

snp_bed <- snps[, .(V1, V2-1, V2)]

fwrite(x= snp_bed, file = "", sep = "\t", col.names=F, row.names=F, quote=F)
