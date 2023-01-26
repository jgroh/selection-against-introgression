library(data.table)


if(Sys.getenv("RSTUDIO") == "1"){
  chromosome <- 4
  #skov_fragments <- fread("Skov_etal_2020_data/41586_2020_2225_MOESM3_ESM.txt")
  sank_calls <- fread("Sankararaman_etal_2014_data/chr-4.thresh-90.length-0.00.gz")
  stein_calls <- fread("Steinrucken_etal_2018/March2018/CEU_lax_chr4/chr4_frqs.txt", col.names = c("chr", "pos", "freq"))
  #skov_fragments <- fread("Skov_etal_2020_data/41586_2020_2225_MOESM3_ESM.txt")
} else {
  args <- commandArgs(trailingOnly = TRUE)
  chromosome <- args[1]
  sank_calls <- fread(args[2])
  stein_calls <- fread(args[3], col.names = c("chr", "pos", "freq"))
  #skov_fragments <- fread(args[4])
}

# column names of sank_calls
# col2: chromosome, col3: genetic position, col4: physical position, col11: average posterior probability of N
setnames(sank_calls, paste0("V", 1:17))
setnames(sank_calls, c("V2", "V3","V4", "V11"), c("chr", "cM", "pos", "freq"))

# ----- interpolate to physical windows
xout <- seq(max(min(stein_calls$pos), min(sank_calls$pos)), min(max(stein_calls$pos), max(sank_calls$pos)), by=1e3)
frq1kb <- sank_calls[, approx(x = pos, y = freq, xout = seq(min(pos), max(pos), by = 1e3)), by = chr]
setnames(frq1kb, c("chr", "pos", "sank_frq"))

frq1kb[, stein_frq := approx(x = stein_calls$pos, y = stein_calls$freq, xout = xout)$y, by = chr]

frq1kb[, Morgan := approx(x = sank_calls$pos, y = sank_calls$cM, xout = xout)$y/100, by = chr]

frq1kb[, rec := (Morgan-shift(Morgan))/1e3]
frq1kb[, rec := approx(x = pos, y = rec, xout = pos, rule = 2)$y, by = chr]

# ----- interpolate to genetic windows
xoutM <- seq(min(sank_calls$cM)/100, max(sank_calls$cM/100), by = 2^-16)

frqM <- sank_calls[, approx(x = cM/100, y = freq, xout = xoutM), by = chr]
setnames(frqM, c("chr", "Morgan", "sank_frq"))


sank_calls[, stein_freq := approx(x = stein_calls$pos, y = stein_calls$freq, xout = sank_calls$pos)$y, by = chr]
frqM[, stein_frq := approx(x = sank_calls$cM/100, y = sank_calls$stein_freq, xout = Morgan)$y, by = chr]

#stein_calls[, Morgan := approx(x = sank_calls$pos, y = sank_calls$cM/100, xout = pos)$y, by = chr]
#frqM[, stein_frq := approx(x = stein_calls$Morgan, y = stein_calls$freq, xout = Morgan)$y, by = chr]

# === output ===== 
fwrite(frq1kb, file = paste0("archaic_freqs_hg19/chr", chromosome, "_frq_physical_windows.txt"), quote = F, sep = "\t")
fwrite(frqM, file = paste0("archaic_freqs_hg19/chr", chromosome, "_frq_genetic_windows.txt"), quote = F, sep = "\t")


