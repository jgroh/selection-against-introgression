library(data.table)
#library(gnomwav)
library(waveslim)
library(wCorr)


if(interactive()){
  windows <- 'physical'
  analysis <- 'wv'
  assembly <- 'hg38'
  thresh <- 'thresh'

} else{

  	source("~/gnomwav/R/multi_modwts.R") 
	source("~/gnomwav/R/variance_decomp.R")
	source("~/gnomwav/R/correlation_decomp.R")
  	args <- commandArgs(trailingOnly = TRUE)
  	windows <- args[1]
  	analysis <- args[2]
  	assembly <- args[3]
  	thresh <- args[4]
}


# ===== read and format data =====

cds_file <- paste0("gene_density_", windows, "_windows_", assembly, ".txt")
Bvals_file <- paste0("B_vals_", windows, "_windows_", assembly, ".txt")
print(Bvals_file)

if(windows == "physical"){
  # read cds files
  cds <- fread(cds_file, col.names = c("chr", "start", "end", "cds"))
  cds[, pos := start + 25e3]
  
  # read frequency files
  chr_files <- dir(path = paste0("archaic_freqs_", assembly, "/"), pattern = "chr.*_frq_physical_windows.txt", full.names=T)
  gnom <- rbindlist(lapply(chr_files, fread))

  # combine gd and freq files
  gnom[, chr := paste0("chr", chrom)]
  gnom <- merge(gnom, cds,  by = c("chr", "pos"))

  # combine with B values (need to interpolate B values )
  Bvals <- fread(Bvals_file, col.names = c("chr", "start", "end", "B"))
  Bvals[, pos := start + 25e3]
  
  B_interp <- Bvals[chr != 'chrX', approx(x = pos, y = B, xout = pos, rule = 2), by = chr]
  setnames(B_interp, c("x", "y"), c("pos", "B"))

  gnom <- merge(gnom, B_interp, by = c("chr", "pos"))

  # log transform of recomb
  gnom[, log10rec := log10(rec)]
 
  #ggplot(gnom[chr== 'chr20'], aes(x = pos, y = log10rec)) + geom_point()
  
  # interpolate to remove inf values
  gnom_complete_obs <- gnom[!is.na(log10rec) & !is.infinite(log10rec) & !is.nan(log10rec)]
  gnom[, log10rec := approx(x = gnom_complete_obs[chr == .BY, pos], 
                            y = gnom_complete_obs[chr == .BY, log10rec],
                            xout = gnom[chr == .BY, pos], rule = 2)$y, by = chr]
                            
  # measure that takes into account both cds and recomb: coding bp per Morgan
  gnom[, cdsM := cds/(Morgan_end-Morgan_start)]
  gnom_complete_obs2 <- gnom[!is.na(cdsM) & !is.infinite(cdsM) & !is.nan(cdsM)]
  gnom[, cdsM := approx(x = gnom_complete_obs2[chr == .BY, pos], 
                            y = gnom_complete_obs2[chr == .BY, cdsM],
                            xout = gnom[chr == .BY, pos], rule = 2)$y, by = chr]
  
  #ggplot(gnom[chr== 'chr20'], aes(x = pos, y = cdsM)) + geom_point()
  setkey(gnom, chr, pos)
  
} else if (windows == "genetic"){
  # read gene density files
  cds <- fread(cds_file, col.names = c("chr", "start", "end", "Morgan", "cds"))
  
  # read frequency files
  chr_files <- dir(path = paste0("archaic_freqs_", assembly, "/"), pattern = "chr.*_frq_genetic_windows.txt", full.names=T)
  gnom <- rbindlist(lapply(chr_files, fread))
  
  # combine gd and freq files
  gnom[, chr := paste0("chr", chrom)]
  gnom <- merge(cds, gnom, by = c("chr", "start", "end", "Morgan"))
  
  # Bvals
  Bvals <- fread(Bvals_file, col.names = c("chr", "start", "end", "Morgan", "B"))
  Bvals <- Bvals[chr != 'chrX'] # if X chrom is present will throw error since there are no values on the X 
  B_interp <- Bvals[, approx(x = Morgan, y = B, xout = Morgan, rule = 2), by = chr]
  setnames(B_interp, c("x", "y"), c("Morgan", "B"))
  gnom <- merge(gnom, B_interp, by = c("chr", "Morgan"))
  
  # fix NA freq values
  # these are places where the Morgan window has zero bp length
  # gnom[is.na(sank_freq), end-start]
  gnom_complete_obs0 <- gnom[!is.na(skov_freq) & !is.infinite(skov_freq) & !is.nan(skov_freq)]
  gnom[, skov_freq := approx(x = gnom_complete_obs0[chr == .BY, Morgan], 
                       y = gnom_complete_obs0[chr == .BY, skov_freq],
                       xout = gnom[chr == .BY, Morgan], rule = 2)$y, by = chr]
  
  
  # fix infinite rec values with interpolation
  # this occurs bc bp interpolation gives same bp for adjacent Morgan intervals
  gnom_complete_obs1 <- gnom[!is.na(rec) & !is.infinite(rec) & !is.nan(rec)]
  gnom[, rec := approx(x = gnom_complete_obs1[chr == .BY, Morgan], 
                            y = gnom_complete_obs1[chr == .BY, rec],
                            xout = gnom[chr == .BY, Morgan], rule = 2)$y, by = chr]
  # log transform of recomb
  gnom[, log10rec := log10(rec)]
  
  #ggplot(gnom[chr== 'chr20'], aes(x = pos, y = log10rec)) + geom_point()

  # interpolate to remove inf values
  gnom_complete_obs2 <- gnom[!is.na(log10rec) & !is.infinite(log10rec) & !is.nan(log10rec)]
  gnom[, log10rec := approx(x = gnom_complete_obs2[chr == .BY, Morgan], 
                            y = gnom_complete_obs2[chr == .BY, log10rec],
                            xout = gnom[chr == .BY, Morgan], rule = 2)$y, by = chr]
  
  # now, just cds base pairs takes into account cds density and recomb rate,
  # make new variable just to run parallel with physical window analysis
  gnom[, cdsM := cds]
  
  #ggplot(gnom[chr== 'chr20'], aes(x = pos, y = cds)) + geom_point()
  setkey(gnom, chr, pos)
  
}

# threshold posterior?

if(thresh == 'thresh'){
  gnom[, sank_freq := sank_freq_thresh]
  gnom[, stein_freq := stein_freq_thresh]
} else{
  gnom[, sank_freq := sank_freq_nothresh]
  gnom[, stein_freq := stein_freq_nothresh]
}

outfile <- paste0("wavelet_results/", analysis, "_", windows, "_", assembly, "_", thresh, ".txt")

# ==== Wavelet Variance ====
if(analysis == "wv"){
  
  wv <- gnom[, gnom_var_decomp(.SD, chromosome = "chr", signals = c("skov_freq", "sank_freq", "stein_freq", "rec", "log10rec", "cds", "cdsM", "B"), rm.boundary = T, avg.over.chroms = T)]
  
  fwrite(wv, file=outfile, sep = "\t", quote=F)
} 

# ===== Wavelet Correlations =====
if(analysis == "wc_calls"){
  wc_skov_sank <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("skov_freq", "sank_freq"), rm.boundary = T)]
  wc_skov_sank[, study := "skov_vs_sank"]
  
  wc_skov_stein <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("skov_freq", "stein_freq"), rm.boundary = T)]
  wc_skov_stein[, study := "skov_vs_stein"]
  
  wc_sank_stein <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("sank_freq", "stein_freq"), rm.boundary = T)]
  wc_sank_stein[, study := "sank_vs_stein"]
  
  wc <- rbindlist(list(wc_skov_sank, wc_sank_stein, wc_skov_stein))
  fwrite(wc, file=outfile, sep = "\t", quote=F)
  
}

if(analysis == "wc_freq_B"){
  wc_skov <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("skov_freq", "B"), rm.boundary = T)]
  wc_skov[, study := "skov"]
 
  wc_sank <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("sank_freq", "B"), rm.boundary = T)]
  wc_sank[, study := "sank"]

  wc_stein <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("stein_freq", "B"), rm.boundary = T)]
  wc_stein[, study := "stein"]
  
  wc <- rbindlist(list(wc_skov, wc_sank, wc_stein))
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_freq_B_chrs"){
  wc_skov <- gnom[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("skov_freq", "B"), rm.boundary = T), by = chr]
  wc_skov[, study := "skov"]
 
  wc_sank <- gnom[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("sank_freq", "B"), rm.boundary = T), by = chr]
  wc_sank[, study := "sank"]

  wc_stein <- gnom[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("stein_freq", "B"), rm.boundary = T), by = chr]
  wc_stein[, study := "stein"]
  
  wc <- rbindlist(list(wc_skov, wc_sank, wc_stein))
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}
if(analysis == "wc_freq_rec"){
  wc_skov <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("skov_freq", "rec"), rm.boundary = T)]
  wc_skov[, study := "skov"]
  
  wc_sank <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("sank_freq", "rec"), rm.boundary = T)]
  wc_sank[, study := "sank"]
  
  wc_stein <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("stein_freq", "rec"), rm.boundary = T)]
  wc_stein[, study := "stein"]
  
  wc <- rbindlist(list(wc_skov, wc_sank, wc_stein))
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_freq_rec_chrs"){
  wc_skov <- gnom[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("skov_freq", "rec"), rm.boundary = T), by = chr]
  wc_skov[, study := "skov"]
  
  wc_sank <- gnom[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("sank_freq", "rec"), rm.boundary = T), by = chr]
  wc_sank[, study := "sank"]
  
  wc_stein <- gnom[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("stein_freq", "rec"), rm.boundary = T), by = chr]
  wc_stein[, study := "stein"]
  
  wc <- rbindlist(list(wc_skov, wc_sank, wc_stein))
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}
if(analysis == "wc_freq_log10rec"){
  wc_skov <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("skov_freq", "log10rec"), rm.boundary = T)]
  wc_skov[, study := "skov"]
  
  wc_sank <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("sank_freq", "log10rec"), rm.boundary = T)]
  wc_sank[, study := "sank"]
  
  wc_stein <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("stein_freq", "log10rec"), rm.boundary = T)]
  wc_stein[, study := "stein"]
  
  wc <- rbindlist(list(wc_skov, wc_sank, wc_stein))
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_freq_cds_chrs"){
  wc_skov <- gnom[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("skov_freq", "cds"), rm.boundary = T), by = chr]
  wc_skov[, study := "skov"]
  
  wc_sank <- gnom[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("sank_freq", "cds"), rm.boundary = T), by = chr]
  wc_sank[, study := "sank"]
  
  wc_stein <- gnom[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("stein_freq", "cds"), rm.boundary = T), by = chr]
  wc_stein[, study := "stein"]
  
  wc <- rbindlist(list(wc_skov, wc_sank, wc_stein))
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_freq_cdsM_chrs"){ # cds per Morgan
  wc_skov <- gnom[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("skov_freq", "cdsM"), rm.boundary = T), by = chr]
  wc_skov[, study := "skov"]
  
  wc_sank <- gnom[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("sank_freq", "cdsM"), rm.boundary = T), by = chr]
  wc_sank[, study := "sank"]
  
  wc_stein <- gnom[, gnom_cor_decomp(.SD, chromosome = NA, signals = c("stein_freq", "cdsM"), rm.boundary = T), by = chr]
  wc_stein[, study := "stein"]
  
  wc <- rbindlist(list(wc_skov, wc_sank, wc_stein))
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

if(analysis == "wc_rec_cds"){
  wc_skov <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("rec", "cds"), rm.boundary = T)]
  wc_skov[, study := "skov"]
  
  wc_sank <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("rec", "cds"), rm.boundary = T)]
  wc_sank[, study := "sank"]
  
  wc_stein <- gnom[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("rec", "cds"), rm.boundary = T)]
  wc_stein[, study := "stein"]
  
  wc <- rbindlist(list(wc_skov, wc_sank, wc_stein))
  fwrite(wc, file=outfile, sep = "\t", quote=F)
}

# ===== Linear Model Analysis =====
if(analysis == "lm"){
  
  z1 <- modwt_lm_rsqrd(data = gnom, chromosome = "chr", yvar = "skov_freq", xvars = c("log10rec"))
  z1[, study := "skov"][, model := "log10rec"][]
  
  z2 <- modwt_lm_rsqrd(data = gnom, chromosome = "chr", yvar = "sank_freq", xvars = c("log10rec"))
  z2[, study := "sank"][, model := "log10rec"][]
  
  z3 <- modwt_lm_rsqrd(data = gnom, chromosome = "chr", yvar = "stein_freq", xvars = c("log10rec"))
  z3[, study := "stein"][, model := "log10rec"][]
  
  z1.1 <- modwt_lm_rsqrd(data = gnom, chromosome = "chr", yvar = "skov_freq", xvars = c("log10rec", "cds"))
  z1.1[, study := "skov"][, model := "log10rec_cds"]
  
  z2.1 <- modwt_lm_rsqrd(data = gnom, chromosome = "chr", yvar = "sank_freq", xvars = c("log10rec", "cds"))
  z2.1[, study := "sank"][, model := "log10rec_cds"]
  
  z3.1 <- modwt_lm_rsqrd(data = gnom, chromosome = "chr", yvar = "stein_freq", xvars = c("log10rec", "cds"))
  z3.1[, study := "stein"][, model := "log10rec_cds"]
  
  z <- rbindlist(list(z1,z2,z3,z1.1,z2.1,z3.1))
  fwrite(z, file=outfile, sep = "\t", quote=F)
}

