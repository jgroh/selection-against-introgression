library(data.table)
library(waveslim)

args <- commandArgs(trailingOnly = TRUE)
windows <- args[1]
analysis <- args[2]

if(Sys.getenv("RSTUDIO") == "1"){
  source("~/workspace/gnomwav/R/multi_modwts.R")
  source("~/workspace/gnomwav/R/variance_decomp.R")
  source("~/workspace/gnomwav/R/correlation_decomp.R")
} else{
  source("/Users/brogroh/gnomwav/R/multi_modwts.R")
  source("/Users/brogroh/gnomwav/R/variance_decomp.R")
  source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
}


# ===== read and format data =====


if(windows == "physical"){
  # read gene density files
  gd <- fread("gene_density_physical_windows.txt", col.names = c("chr", "start", "end", "gd"))
  gd[, pos := start + 500]
  
  Bvals <- fread("B_vals_physical_windows.txt", col.names = c("chr", "start", "end", "B"))
  Bvals[, pos := start + 500]
  
  # read frequency files
  chr_files <- dir(path = "archaic_freqs/", pattern = "chr.*_frq_physical_windows.txt", full.names=T)
  gnom <- rbindlist(lapply(chr_files, fread))
  
  # combine gd and freq files
  gnom[, chr := paste0("chr", chr)]
  gnom <- merge(gnom, merge(gd, Bvals, by = c("chr", "pos", "start", "end")), by = c("chr", "pos"))
  gnom[, gdr := gd/rec]
  gnom[gdr == Inf, gdr := gnom[gdr != Inf, max(gdr)]]
  gnom[gd==0 & rec==0, gdr := gnom[, mean(gdr, na.rm=T)]][]
  setkey(gnom, chr, pos)
  
} else if (windows == "genetic"){
  # read gene density files
  gd <- fread("gene_density_genetic_windows.txt", col.names = c("chr", "start", "end", "Morgan", "gd"))
  
  # read frequency files
  chr_files <- dir(path = "archaic_freqs/", pattern = "chr.*_frq_genetic_windows.txt", full.names=T)
  gnom <- rbindlist(lapply(chr_files, fread))
  
  Bvals <- fread("B_vals_genetic_windows.txt")
  
  # combine gd and freq files
  gnom[, chr := paste0("chr", chr)]
  gnom <- merge(gd, gnom, by = c("chr", "Morgan"))
  gnom[, gdr := gd/rec]
  gnom[gdr == Inf, gdr := gnom[gdr != Inf, max(gdr)]]
  gnom[gd == 0 & rec == 0, gdr := gnom[, mean(gdr, na.rm=T)]][]
  setkey(gnom, chr, pos)
  
}

if(Sys.getenv("RSTUDIO") == "1"){
  gnom <- gnom[chr %in% paste0("chr", 10:15)][seq(1, .N, by=100)]
} 

outfile <- paste0("wavelet_results/", analysis, "_", windows, ".txt")

# ==== Wavelet Variance ====
if(analysis == "wv"){
  wv <- gnom[, gnom_var_decomp(.SD, chromosome = "chr", signals = c("freq", "rec", "gd"), rm.boundary = T, avg.over.chroms = T)]
  
  fwrite(wv, file=outfile, sep = "\t", quote=F)
} 

# ===== Wavelet Correlations =====
if(analysis == "wc_freq_B"){
  gnom[, B := gnom[, approx(x = pos, y = B, xout=pos, rule = 2), by = chr]$y]

  wc <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "B"), rm.boundary = T)]
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_freq_rec"){
  wc <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T)]
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_freq_gd"){
  wc <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "gd"), rm.boundary = T)]
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_freq_gdr"){
  wc <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "gdr"), rm.boundary = T)]
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_rec_gd"){
  wc <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("rec", "gd"), rm.boundary = T)]
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

# ===== Linear Model Analysis =====
if(analysis == "lm"){
  z <- wvlt_lm(data = gnom, chromosome = "chr", yvar = "freq", xvars = c("B"))
  fwrite(z, file=outfile, sep = "\t", quote=F)
}

