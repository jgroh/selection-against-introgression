library(data.table)

frqList <- list()
frqListM <- list()

for(i in 1:22){
  #filename <- "chr-9.thresh-90.length-0.00" # for local testing
  
  filename <- paste0("summaries.release/EUR-ASN.hapmap/summaries/chr-", i, ".thresh-90.length-0.00")
  
  # col 2 is chromosome, 4 is physical position, col 11 is average marginal posterior probability
  dt <- fread(filename)
  
  # interpolate to physical grid. first try 1kb for comparison to skov analysis
  frq <- dt[, approx(xout = seq(min(V4), max(V4), by = 1e3), x = V4, y = V11), by = V2]
  setnames(frq, c("chr", "pos", "frq"))
  
  rec <- dt[, approx(xout = seq(min(V4), max(V4), by = 1e3), x = V4, y = Morgan), by = V2]
  setnames(rec, c("chr", "pos", "Morgan"))
  rec[, rec := (Morgan - shift(Morgan))/(pos-shift(pos))]
  rec[, rec := c(.SD[2, rec], .SD[2:nrow(.SD), rec])]
  frq <- merge(frq, rec)
  frqList[[i]] <- frq
  
  # interpolate to genetic grid
  dt[, Morgan := V3 /100]
  frqM <- dt[, approx(xout = seq(min(Morgan), max(Morgan), by = 2^-16), x = Morgan, y = V11), by = V2]
  setnames(frqM, c("chr", "Morgan", "frq"))
  frqM
  
  recM <- dt[, approx(xout = seq(min(Morgan), max(Morgan), by = 2^-16), x = Morgan, y= V4), by = V2]
  setnames(recM, c("chr", "Morgan", "pos"))
  recM[, rec := (Morgan - shift(Morgan))/(pos-shift(pos))]
  recM[, rec := c(.SD[2, rec], .SD[2:nrow(.SD), rec])]
  frqM <- merge(frqM, recM[, .(chr, Morgan, rec)])
  frqListM[[i]] <- frqM
}
  
frqs <- rbindlist(frqList)
frqsM <- rbindlist(frqListM)

fwrite(frqs, file = "interpolated_frq_physical_map.txt", quote=F, sep = "\t", col.names=T, row.names = F)
fwrite(frqsM, file = "interpolated_frq_genetic_map.txt", quote=F, sep = "\t", col.names=T, row.names = F)

