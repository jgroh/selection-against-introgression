library(data.table)
args <- commandArgs(trailingOnly = T)

outfile <- args[1]
data <- fread(args[2])

for(i in 2:length(args)){
  data <- rbind(data, fread(args[i]))
}

write.table(data, file = outfile, sep = "\t", row.names = F, col.names = T)