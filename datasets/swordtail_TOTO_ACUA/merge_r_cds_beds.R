library(data.table)

units <- commandArgs(trailingOnly = T)[1]

cds_file <- paste0("xbir_cds_", units, "_unit_windows.bed")
r_file <- paste0("xbir_r_", units, "_unit_windows.bed")

cds <- fread(cds_file, col.names = c("chr", "start", "end", "cds_density"))
cm <- fread(r_file, col.names = c("chr", "start", "end", "r"))

write.table(merge(cds,cm), 
			file = paste0("xbir_r_cds_", units, "_unit_windows.txt"),
		   	quote=F, sep = "\t", row.names=F,col.names=T)
