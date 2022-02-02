# 1. ========== Load Dependencies and Data ==========
library(data.table)
library(tools)
library(waveslim)
library(stringi)
library(magrittr)

#source("~/gnomwav/R/correlation_decomp.R")
#source("~/gnomwav/R/multi_modwts.R")
#source("~/gnomwav/R/variance_decomp.R")

# for testing locally
 source("~/workspace/gnomwav/R/correlation_decomp.R")
 source("~/workspace/gnomwav/R/multi_modwts.R")
 source("~/workspace/gnomwav/R/variance_decomp.R")

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]
#year <- 2018

# ===== Read and format data =====

# read cds density file
cdsCM <- fread("xbir_1kb_CDS_and_cM.txt")

# read chrom lengths file
chrLenP <- fread("xbir10x_chrlengths.txt", col.names = c("chr", "len"))

# interpolated ancestry files for each chromosome
scaffFiles <- paste0("ACUA_2018/",paste0(chrLenP$chr, ".RData"))
# scaffFiles <- dir(path=paste0("ACUA_",year),pattern=".RData",full.names=T) # older version, useful if not all chromosome files present

scaffFiles <- scaffFiles[1:5] # for testing locally
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
cdsCM[, position := as.integer(start + 500)] # this is the midpoint of the 1kb intervals where I interpolate ancestry
gnomP <- merge(gnomP, cdsCM, by = c("chr", "position"))

# log transform recombination rates
gnomP[, cmTr := log10(cM)]


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
                                   signals = c("meanFreq", "coding_bp", "cmTr"),
                                   rm.boundary = TRUE, 
                                   avg.over.chroms = TRUE)]
# ------ mean ancestry, recombination, gene density -----
meanWV_G <- gnomG[ID==ID[1], gnom_var_decomp(.SD, chromosome = "chr", 
                                             signals = c("meanFreq"),
                                             rm.boundary = TRUE, 
                                             avg.over.chroms = TRUE)]

meanWV_P[, units := "physical"]
meanWV_P[, signal := NULL]
meanWV_G[, units := "genetic"]
meanWV_G[, signal := NULL]

wavvar_G <- merge(meanWV_G, indWV_G)
wavvar_P <- merge(meanWV_P, indWV_P)

# # ====== Correlation decomp =====
# 
cortbl1 <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "cmTr"))]
cortbl2 <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "coding_bp"))]
cortbl3 <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("cmTr", "coding_bp"))]

wavcor <- merge(cortbl3, merge(cortbl1,cortbl2))


# # ===== Linear Model analysis =====

rsqrd <- wvlt_lm_rsqrd(data = gnomP[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
              xvars = c("cmTr", "coding_bp"))

# # ====== Save output =======
# 
# 
# save(wavvar, wavcor, rsqrd, file = paste0("ACUA_", year, "/wavelet_results.RData"))
