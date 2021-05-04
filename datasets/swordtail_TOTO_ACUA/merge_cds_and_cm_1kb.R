library(data.table)

cds <- fread("cds_per_kb.bed", col.names = c("chr", "start", "end", "coding_bp"))
cm <- fread("allChr_1kb_rec.bed", col.names = c("chr", "start", "end", "cM"))

write.table(merge(cds,cm), file = "allChr_cds_and_cm_per_1kb.txt", quote=F, sep = "\t", row.names=F,col.names=T)

