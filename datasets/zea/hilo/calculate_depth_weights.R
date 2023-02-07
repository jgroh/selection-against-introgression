library(data.table)

ID <- commandArgs(trailingOnly=T)[1]

depth <- fread('cat /dev/stdin', col.names = c("chr", "pos", "depth"))

w <- depth[,mean((2*depth)/(depth+1))]
fwrite(data.table(ID, w), file = "", sep = "\t", quote=F, col.names=F, row.names=F)

