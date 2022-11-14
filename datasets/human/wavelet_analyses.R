library(data.table)
library(waveslim)

if(interactive()){
  source("~/workspace/gnomwav/R/multi_modwts.R")
  source("~/workspace/gnomwav/R/variance_decomp.R")
  source("~/workspace/gnomwav/R/correlation_decomp.R")
  windows <- 'physical'
  analysis <- 'wv'

} else{
  source("/Users/brogroh/gnomwav/R/multi_modwts.R")
  source("/Users/brogroh/gnomwav/R/variance_decomp.R")
  source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
  args <- commandArgs(trailingOnly = TRUE)
  windows <- args[1]
  analysis <- args[2]
}


# ===== read and format data =====


if(windows == "physical"){
  # read gene density files
  gd <- fread("gene_density_physical_windows.txt", col.names = c("chr", "start", "end", "gd"))
  gd[, pos := start + 500]
  
  # read frequency files
  chr_files <- dir(path = "archaic_freqs/", pattern = "chr.*_frq_physical_windows.txt", full.names=T)
  gnom <- rbindlist(lapply(chr_files, fread))
  
  # combine gd and freq files
  gnom[, chr := paste0("chr", chr)]
  gnom <- merge(gnom, gd,  by = c("chr", "pos"))

  # combine with B values (need to interpolate B values )
  Bvals <- fread("B_vals_physical_windows.txt", col.names = c("chr", "start", "end", "B"))
  Bvals[, pos := start + 500]
  
  B_interp <- Bvals[, approx(x = pos, y = B, xout = pos, rule = 2), by = chr]
  setnames(B_interp, c("x", "y"), c("pos", "B"))
  gnom <- merge(gnom, B_interp)
  
  gnom[, log10rec := log10(rec)]
  gnom[is.infinite(log10rec), log10rec := min(gnom[!is.infinite(log10rec), log10rec])]
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
  
  # combine gd and freq files
  gnom[, chr := paste0("chr", chr)]
  gnom <- merge(gd, gnom, by = c("chr", "Morgan"))
  
  # Bvals
  Bvals <- fread("B_vals_genetic_windows.txt", col.names = c("chr", "start", "end", "Morgan", "B"))
  B_interp <- Bvals[, approx(x = Morgan, y = B, xout = Morgan, rule = 2), by = chr]
  setnames(B_interp, c("x", "y"), c("Morgan", "B"))
  gnom <- merge(gnom, B_interp)
  
  gnom[, log10rec := log10(rec)]
  gnom[is.infinite(log10rec), log10rec := min(gnom[!is.infinite(log10rec), log10rec])]
  
  gnom[, gdr := gd/rec]
  gnom[gdr == Inf, gdr := gnom[gdr != Inf, max(gdr)]]
  gnom[gd == 0 & rec == 0, gdr := gnom[, mean(gdr, na.rm=T)]]
  setkey(gnom, chr, pos)
  
}

if(interactive()){
  gnom <- gnom[chr %in% paste0("chr", 10:15)][seq(1, .N, by=100)]
} 

outfile <- paste0("wavelet_results/", analysis, "_", windows, ".txt")

# ==== Wavelet Variance ====
if(analysis == "wv"){
  
  wv <- gnom[, gnom_var_decomp(.SD, chromosome = "chr", signals = c("freq", "rec", "gd", "B"), rm.boundary = F, avg.over.chroms = T)]
  
  fwrite(wv, file=outfile, sep = "\t", quote=F)
} 

# ===== Wavelet Correlations =====
if(analysis == "wc_freq_B"){
  wc <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "B"), rm.boundary = F)]
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_freq_rec"){
  wc <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = F)]
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_freq_log10rec"){
  wc <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "log10rec"), rm.boundary = F)]
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_freq_gd"){
  wc <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "gd"), rm.boundary = F)]
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_freq_gdr"){
  wc <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "gdr"), rm.boundary = F)]
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_rec_gd"){
  wc <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("rec", "gd"), rm.boundary = F)]
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

# ===== Linear Model Analysis =====
if(analysis == "lm"){
  z1 <- modwt_lm_rsqrd(data = gnom, chromosome = "chr", yvar = "freq", xvars = c("rec", "gd"))
  z1[, model := "rec_cds_density"]
  z2 <- modwt_lm_rsqrd(data = gnom, chromosome = "chr", yvar = "freq", xvars = c("log10rec", "gd"))
  z2[, model := "log10rec_cds_density"]
  
  z <- rbind(z1, z2)
  fwrite(z, file=outfile, sep = "\t", quote=F)
}

