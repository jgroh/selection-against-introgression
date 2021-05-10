library(data.table)

cds <- fread("xbir_1kb_CDS.bed", col.names = c("chr", "start", "end", "coding_bp"))
cm <- fread("xbir_1kb_rec.bed", col.names = c("chr", "start", "end", "cM"))

write.table(merge(cds,cm), file = "xbir_1kb_CDS_and_cM.txt", quote=F, sep = "\t", row.names=F,col.names=T)

