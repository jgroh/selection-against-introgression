# 1. ========== Load Dependencies and Data ==========
library(data.table)
library(tools)
library(waveslim)
library(stringi)

source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
source("/Users/brogroh/gnomwav/R/multi_modwts.R")
source("/Users/brogroh/gnomwav/R/variance_decomp.R")

# for testing locally
#source("~/workspace/gnomwav/R/correlation_decomp.R")
#source("~/workspace/gnomwav/R/multi_modwts.R")
#source("~/workspace/gnomwav/R/variance_decomp.R")

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]
#year <- 2013

# ===== Read and format data =====

# read r and cds density files

g_features <-  fread("xbir_r_cds_genetic_unit_windows.txt")
p_features <- fread("xbir_r_cds_physical_unit_windows.txt")

# read chrom lengths file
chrLenP <- fread("xbir10x_chrlengths.txt", col.names = c("chr", "len"))

# interpolated ancestry files for each chromosome
scaffFiles <- paste0("ACUA_", year, "/", paste0(chrLenP$chr, ".RData"))
# scaffFiles <- dir(path=paste0("ACUA_",year),pattern=".RData",full.names=T) # older version, useful if not all chromosome files present
# scaffFiles <- scaffFiles[1:4] # for testing interactively

scaffs <- basename(file_path_sans_ext(scaffFiles))


# each file has the same object name for the chromosome ancestry table
# so we load into separate environments
loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 

# combine data into one file for the year
gnomG <- rbindlist(lapply(as.list(scaffFiles), 
                  function(x){loadFrom(x, "chromAncInterpMorgan")}))

gnomP <- rbindlist(lapply(as.list(scaffFiles), 
                             function(x){loadFrom(x, "chromAnc1kb")}))

# For ACUA, minor parent is malinche
gnomG[, indivFreq := 1 - indivFreq]
gnomG[, meanFreq := 1 - meanFreq]

gnomP[, indivFreq := 1 - indivFreq]
gnomP[, meanFreq := 1 - meanFreq]

# merge ancestry data and genomic features data
# due to how the ancestry tables were constructured, merge on different variables for physical and genetic units
# for physical scale, merge on position
p_features[, position := as.integer(start + 500)]
gnomP <- merge(gnomP, p_features, by = c("chr", "position"))

g_features[, end := c(end[1:(.N-1)], end[.N]-1), by = chr]
gnomG <- merge(gnomG, g_features, by = c("chr", "end"))

# log transform recombination rates
gnomP[, log10r := log10(r)]
gnomG[, log10r := log10(r)]


# ===== Variance decomp =====

# individual-level ancestry

# physical units
indWV_P <- gnomP[, gnom_var_decomp(.SD, chromosome = "chr", 
                        signals = c("indivFreq"),
                        rm.boundary = TRUE, 
                        avg.over.chroms = TRUE),
      by = ID]


# genetic units
indWV_G <- gnomG[, gnom_var_decomp(.SD, chromosome = "chr", 
                                   signals = c("indivFreq"),
                                   rm.boundary = TRUE, 
                                   avg.over.chroms = TRUE),
                 by = ID]

# average over individuals
indWV_P <- indWV_P[, lapply(.SD, mean), 
                   .SDcols = c("variance.indivFreq", "propvar.indivFreq"), by = level]


# average over individuals
indWV_G <- indWV_G[, lapply(.SD, mean), 
                   .SDcols = c("variance.indivFreq", "propvar.indivFreq"), by = level]

indWV_G[, units := "genetic"]
indWV_P[, units := "physical"]


# ------ mean ancestry, recombination, gene density -----
meanWV_P <- gnomP[ID==ID[1], gnom_var_decomp(.SD, chromosome = "chr", 
                                   signals = c("meanFreq", "cds_density", "r", "log10r", "dist_to_marker"),
                                   rm.boundary = TRUE, 
                                   avg.over.chroms = TRUE)]

# ------ mean ancestry, recombination, gene density -----
meanWV_G <- gnomG[ID==ID[1], gnom_var_decomp(.SD, chromosome = "chr", 
                                             signals = c("meanFreq", "cds_density", "r", "log10r", "dist_to_marker"),
                                             rm.boundary = TRUE, 
                                             avg.over.chroms = TRUE)]

meanWV_P[, units := "physical"]
meanWV_G[, units := "genetic"]

wavvar_G <- merge(meanWV_G, indWV_G)
wavvar_P <- merge(meanWV_P, indWV_P)

wavvar <- rbind(wavvar_G, wavvar_P)

# # ====== Correlation decomp =====
# ---- calculate correlations between ancestry and both recombination and gene density. since recombination and gene
# density are themselves correlated, also calculate partial correlations

# ----- physical units
oldnames <- c("cor_n", "cor_jack", "lower95ci", "upper95ci")

cortbl1P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "r"))]
setnames(cortbl1P, old = oldnames, new = paste0(oldnames, ".meanFreq_r"))

cortbl2P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "cds_density"))]
setnames(cortbl2P, old = oldnames, new = paste0(oldnames, ".meanFreq_cds"))

cortbl3P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("r", "cds_density"))]
setnames(cortbl3P, old = oldnames, new = paste0(oldnames, ".r_cds"))


# 
# # for partial correlation of freq and recombination, regress each against cds and take residuals
# gnomP[ID==ID[1], frq_cds_res := residuals(lm(meanFreq~cds_density))]
# gnomP[ID==ID[1], r_cds_res := residuals(lm(r~cds_density))]
# gnomP[ID==ID[1], log10r_cds_res := residuals(lm(log10r~cds_density))]
# 
# cortbl1P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "r"))]
# setnames(cortbl1P, old = oldnames, new = paste0(oldnames, ".meanFreq_r"))
# 
# cortbl1.1P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("frq_cds_res", "r_cds_res"))]
# setnames(cortbl1.1P, old = oldnames, new = paste0(oldnames, ".meanFreq_r_prtl"))
# 
# cortbl1.2P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("frq_cds_res", "log10r_cds_res"))]
# setnames(cortbl1.2P, old = oldnames, new = paste0(oldnames, ".meanFreq_log10r_prtl"))
# 
# # for partial correlation of freq and cds, regress each against recomb
# gnomP[ID==ID[1], frq_r_res := residuals(lm(meanFreq~r))]
# gnomP[ID==ID[1], frq_log10r_res := residuals(lm(meanFreq~log10r))]
# gnomP[ID==ID[1], cds_r_res := residuals(lm(cds_density~r))]
# gnomP[ID==ID[1], cds_log10r_res := residuals(lm(cds_density~log10r))]
# 
# cortbl2P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("frq_r_res", "cds_r_res"))]
# setnames(cortbl2P, old = oldnames, new = paste0(oldnames, ".meanFreq_cdsDensity.resid_r"))
# 
# cortbl2.1P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("frq_log10r_res", "cds_log10r_res"))]
# setnames(cortbl2.1P, old = oldnames, new = paste0(oldnames, ".meanFreq_cdsDensity.resid_log10r"))
# 
# # for correlation between recombination and gene density, use straight up correlation
# cortbl3P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("r", "cds_density"))]
# setnames(cortbl3P, old = oldnames, new = paste0(oldnames, ".r_cdsDensity"))
# 
# cortbl3.1P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("log10r", "cds_density"))]
# setnames(cortbl3.1P, old = oldnames, new = paste0(oldnames, ".log10r_cdsDensity"))
# 
# wavcorP <- merge(merge(merge(cortbl1P, cortbl1.1P), merge(cortbl2P, cortbl2.1P)), merge(cortbl3P, cortbl3.1P))

wavcorP <- merge( merge(cortbl1P, cortbl2P),  cortbl3P)
wavcorP[, units := "physical"]



# ----- genetic units

cortbl1G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "r"))]
setnames(cortbl1G, old = oldnames, new = paste0(oldnames, ".meanFreq_r"))

cortbl2G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "cds_density"))]
setnames(cortbl2G, old = oldnames, new = paste0(oldnames, ".meanFreq_cds"))

cortbl3G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("r", "cds_density"))]
setnames(cortbl3G, old = oldnames, new = paste0(oldnames, ".r_cds"))

# 
# # for partial correlation of freq and recombination, regress each against cds and take residuals
# gnomG[ID==ID[1], frq_cds_res := residuals(lm(meanFreq~cds_density))]
# gnomG[ID==ID[1], r_cds_res := residuals(lm(r~cds_density))]
# gnomG[ID==ID[1], log10r_cds_res := residuals(lm(log10r~cds_density))]
# 
# cortbl1.1G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("frq_cds_res", "r_cds_res"))]
# setnames(cortbl1.1G, old = oldnames, new = paste0(oldnames, ".meanFreq_r_prtl"))
# 
# cortbl1.2G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("frq_cds_res", "log10r_cds_res"))]
# setnames(cortbl1.2G, old = oldnames, new = paste0(oldnames, ".meanFreq_log10r_prtl"))
# 
# # for partial correlation of freq and cds, regress each against recomb
# gnomG[ID==ID[1], frq_r_res := residuals(lm(meanFreq~r))]
# gnomG[ID==ID[1], frq_log10r_res := residuals(lm(meanFreq~log10r))]
# gnomG[ID==ID[1], cds_r_res := residuals(lm(cds_density~r))]
# gnomG[ID==ID[1], cds_log10r_res := residuals(lm(cds_density~log10r))]
# 
# cortbl2G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("frq_r_res", "cds_r_res"))]
# setnames(cortbl2G, old = oldnames, new = paste0(oldnames, ".meanFreq_cdsDensity.resid_r"))
# 
# cortbl2.1G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("frq_log10r_res", "cds_log10r_res"))]
# setnames(cortbl2.1G, old = oldnames, new = paste0(oldnames, ".meanFreq_cdsDensity.resid_log10r"))
# 
# # for correlation between recombination and gene density, use straight up correlation
# cortbl3G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("r", "cds_density"))]
# setnames(cortbl3G, old = oldnames, new = paste0(oldnames, ".r_cdsDensity"))
# 
# cortbl3.1G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("log10r", "cds_density"))]
# setnames(cortbl3.1G, old = oldnames, new = paste0(oldnames, ".log10r_cdsDensity"))

#wavcorG <- merge(merge(merge(cortbl1G, cortbl1.1G), merge(cortbl2G, cortbl2.1G)), merge(cortbl3G, cortbl3.1G))

wavcorG <- merge( merge(cortbl1G, cortbl2G),  cortbl3G)

wavcorG[, units := "genetic"]


wavcor <- rbind(wavcorG, wavcorP)


# # ===== Linear Model analysis =====

# ----physical units
rsqrd1_p <- wvlt_lm(data = gnomP[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
                         xvars = c("r", "cds_density"))
rsqrd1_p[, model := "r_cdsDensity"]

rsqrd2_p <- wvlt_lm(data = gnomP[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
              xvars = c("log10r", "cds_density"))
rsqrd2_p[, model := "log10r_cdsDensity"]

rsqrd_p <- rbind(rsqrd1_p, rsqrd2_p)
rsqrd_p[, units := "physical"]

# ---- genetic units
rsqrd1_g <- wvlt_lm(data = gnomG[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
                          xvars = c("r", "cds_density"))
rsqrd1_g[, model := "r_cdsDensity"]

rsqrd2_g <- wvlt_lm(data = gnomG[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
                          xvars = c("log10r", "cds_density"))
rsqrd2_g[, model := "log10r_cdsDensity"]

rsqrd_g <- rbind(rsqrd1_g, rsqrd2_g)
rsqrd_g[, units := "genetic"]


rsqrd <- rbind(rsqrd_p, rsqrd_g)

# # ====== Save output =======
save(wavvar, wavcor, rsqrd, file = paste0("ACUA_", year, "/wavelet_results.RData"))
