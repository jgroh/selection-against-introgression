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
#year <- 2018

# ===== Read and format data =====

# read r and cds density files

g_features <-  fread("xbir_r_cds_genetic_unit_windows.txt")
p_features <- fread("xbir_r_cds_physical_unit_windows.txt")

# read chrom lengths file
chrLenP <- fread("xbir10x_chrlengths.txt", col.names = c("chr", "len"))

# interpolated ancestry files for each chromosome
scaffFiles <- paste0("ACUA_2018/",paste0(chrLenP$chr, ".RData"))
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
                                   signals = c("meanFreq", "cds_density", "log10r"),
                                   rm.boundary = TRUE, 
                                   avg.over.chroms = TRUE)]

# ------ mean ancestry, recombination, gene density -----
meanWV_G <- gnomG[ID==ID[1], gnom_var_decomp(.SD, chromosome = "chr", 
                                             signals = c("meanFreq", "cds_density", "log10r"),
                                             rm.boundary = TRUE, 
                                             avg.over.chroms = TRUE)]

meanWV_P[, units := "physical"]
meanWV_G[, units := "genetic"]

wavvar_G <- merge(meanWV_G, indWV_G)
wavvar_P <- merge(meanWV_P, indWV_P)

wavvar <- rbind(wavvar_G, wavvar_P)

# # ====== Correlation decomp =====

# physical units
cortbl1P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "log10r"))]
setnames(cortbl1P, "cor", "cor_meanFreq_log10r")
cortbl2P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "cds_density"))]
setnames(cortbl2P, "cor", "cor_meanFreq_cdsDensity")
cortbl3P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("log10r", "cds_density"))]
setnames(cortbl3P, "cor", "cor_log10r_cdsDensity")

wavcorP <- merge(cortbl3P, merge(cortbl1P,cortbl2P), by = "level")
wavcorP[, units := "physical"]

# genetic units
cortbl1G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "log10r"))]
setnames(cortbl1G, "cor", "cor_meanFreq_log10r")
cortbl2G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "cds_density"))]
setnames(cortbl2G, "cor", "cor_meanFreq_cdsDensity")
cortbl3G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("log10r", "cds_density"))]
setnames(cortbl3G, "cor", "cor_log10r_cdsDensity")

wavcorG <- merge(cortbl3G, merge(cortbl1G,cortbl2G), by = "level")
wavcorG[, units := "genetic"]

wavcor <- rbind(wavcorG, wavcorP)


#save(wavcor, file = paste0("ACUA_", year, "/wavcor_results.RData"))

# # ===== Linear Model analysis =====

rsqrd_p <- wvlt_lm_rsqrd(data = gnomP[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
              xvars = c("log10r", "cds_density"))

rsqrd_g <- wvlt_lm_rsqrd(data = gnomG[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
                         xvars = c("log10r", "cds_density"))

rsqrd <- rbind(rsqrd_p, rsqrd_g)

# # ====== Save output =======
save(wavvar, wavcor, rsqrd, file = paste0("ACUA_", year, "/wavelet_results.RData"))
