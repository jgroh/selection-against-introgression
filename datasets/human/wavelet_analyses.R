library(data.table)
library(waveslim)

#source("~/workspace/gnomwav/R/multi_modwts.R")
#source("~/workspace/gnomwav/R/variance_decomp.R")
#source("~/workspace/gnomwav/R/correlation_decomp.R")

source("/Users/brogroh/gnomwav/R/multi_modwts.R")
source("/Users/brogroh/gnomwav/R/variance_decomp.R")
source("/Users/brogroh/gnomwav/R/correlation_decomp.R")

# read gene density files
gd1kb <- fread("gene_density_1kb_windows.txt", col.names = c("chr", "start", "end", "gd"))
gd1kb[, pos := start + 500]
gdM <- fread("gene_density_genetic_windows.txt", col.names = c("chr", "start", "end", "Morgan", "gd"))

# read frequency files
chr_files1kb <- dir(path = "archaic_freqs/", pattern = "chr.*_frq_1kb_windows.txt", full.names=T)
chr_filesM <- dir(path = "archaic_freqs/", pattern = "chr.*_frq_genetic_windows.txt", full.names=T)
gnom1kb <- rbindlist(lapply(chr_files1kb, fread))
gnomM <- rbindlist(lapply(chr_filesM, fread))


# combine gd and freq files
gnom1kb[, chr := paste0("chr", chr)]
gnom1kb[, gdr := gd/rec]
setkey(gnom1kb, chr, pos)
gnom1kb <- merge(gd1kb, gnom1kb, by = c("chr", "pos"))

gnomM[, chr := paste0("chr", chr)]
gnomM[, gdr := gd/rec]
setkey(gnomM, chr, pos)
gnomM <- merge(gdM, gnomM, by = c("chr", "Morgan"))

# ===== Wavelet Variances =====
wv1kb <- gnom1kb[, gnom_var_decomp(.SD, chromosome = "chr", signals = c("freq", "rec", "gd"), rm.boundary = T, avg.over.chroms = T)]

wvM <- gnomM[, gnom_var_decomp(.SD, chromosome = "chr", signals = c("freq", "rec", "gd"), rm.boundary = T, avg.over.chroms = T)]


# ===== Wavelet Correlations =====
wc1kb_freq_rec <- gnom1kb[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T)]
wcM_freq_rec <- gnomM[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T)]

wc1kb_freq_gdr <- gnom1kb[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "gdr"), rm.boundary = T)]
wcM_freq_gdr <- gnomM[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "gdr"), rm.boundary = T)]


# ===== Lm results =====

lm1kb <- wvlt_lm(data = gnom1kb, chromosome = "chr", yvar = "freq", xvars = c("rec", "gd"))
lmM <- wvlt_lm(data = gnomM, chromosome = "chr", yvar = "freq", xvars = c("rec", "gd"))


# ===== output =====

fwrite(wv1kb, file = "wavlet_results/wv1kb.txt", sep = "\t", quote = F)
fwrite(wvM, file = "wavlet_results/wvM.txt", sep = "\t", quote = F)

fwrite(wc1kb_freq_rec, file = "wavlet_results/wc1kb_freq_rec.txt", sep = "\t", quote = F)
fwrite(wcM_freq_rec, file = "wavlet_results/wcM_freq_rec.txt", sep = "\t", quote = F)

fwrite(wc1kb_freq_gdr, file = "wavlet_results/wc1kb_freq_gdr.txt", sep = "\t", quote = F)
fwrite(wcM_freq_gdr, file = "wavlet_results/wcM_freq_gdr.txt", sep = "\t", quote = F)

fwrite(lm1kb, file = "wavlet_results/lm1kb.txt", sep = "\t", quote = F)
fwrite(lmM, file = "wavlet_results/lmM.txt", sep = "\t", quote = F)
