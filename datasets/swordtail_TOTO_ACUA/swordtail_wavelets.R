# 1. ========== Load Dependencies and Data ==========
library(data.table)
library(tools)
library(waveslim)
library(stringi)
library(wCorr)

source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
source("/Users/brogroh/gnomwav/R/multi_modwts.R")
source("/Users/brogroh/gnomwav/R/variance_decomp.R")

# for testing locally
# source("~/workspace/gnomwav/R/correlation_decomp.R")
# source("~/workspace/gnomwav/R/multi_modwts.R")
# source("~/workspace/gnomwav/R/variance_decomp.R")

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
#scaffFiles <- paste0("ACUA_", year, "/", paste0(chrLenP$chr, ".RData"))
scaffFiles <- paste0(year, "/", paste0(chrLenP$chr, ".RData"))
# scaffFiles <- scaffFiles[1:5] # for testing interactively

scaffs <- basename(file_path_sans_ext(scaffFiles))


# each file has the same object name for the chromosome ancestry table
# so we load into separate environments
loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]}

# combine data into one file for the year
gnomG <- rbindlist(lapply(as.list(scaffFiles),
                  function(x){loadFrom(x, "chromAncInterpMorgan")}))

gnomP <- rbindlist(lapply(as.list(scaffFiles),
                             function(x){loadFrom(x, "chromAnc50kb")}))

# For ACUA, minor parent is malinche
gnomG[, indivFreq := 1 - indivFreq]
gnomG[, meanFreq := 1 - meanFreq]

gnomP[, indivFreq := 1 - indivFreq]
gnomP[, meanFreq := 1 - meanFreq]

# merge ancestry data and genomic features data
# due to how the ancestry tables were constructured, merge on different variables for physical and genetic units
# for physical scale, merge on position
p_features[, position := as.integer(start + 25e3)]
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
                        rm.boundary = FALSE,
                        avg.over.chroms = TRUE),
      by = ID]


# genetic units
indWV_G <- gnomG[, gnom_var_decomp(.SD, chromosome = "chr",
                                   signals = c("indivFreq"),
                                   rm.boundary = FALSE,
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
                                   signals = c("meanFreq", "cds_density", "r", "log10r"),
                                   rm.boundary = FALSE,
                                   avg.over.chroms = TRUE)]

# ------ mean ancestry, recombination, gene density -----
meanWV_G <- gnomG[ID==ID[1], gnom_var_decomp(.SD, chromosome = "chr",
                                             signals = c("meanFreq", "cds_density", "r", "log10r"),
                                             rm.boundary = FALSE,
                                             avg.over.chroms = TRUE)]

meanWV_P[, units := "physical"]
meanWV_G[, units := "genetic"]

wavvar_G <- merge(meanWV_G, indWV_G)
wavvar_P <- merge(meanWV_P, indWV_P)

wavvar <- rbind(wavvar_G, wavvar_P)

# # ====== Correlation decomp =====
# ---- calculate correlations between ancestry and both recombination and gene density.

# ----- physical units

cortbl1P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "r"))]
cortbl1P[, vars := "meanFreq_r"]

cortbl1P.1 <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "log10r"))]
cortbl1P.1[, vars := "meanFreq_log10r"]

cortbl2P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "cds_density"))]
cortbl2P[, vars := "meanFreq_cds_density"]

cortbl3P <- gnomP[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("r", "cds_density"))]
cortbl3P[, vars := "r_cds_density"]


wavcorP <- rbindlist(list(cortbl1P, cortbl1P.1,cortbl2P,cortbl3P))
wavcorP[, units := "physical"]



# ----- genetic units

cortbl1G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "r"))]
cortbl1G[, vars := "meanFreq_r"]

cortbl1G.1 <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "log10r"))]
cortbl1G.1[, vars := "meanFreq_logr"]

cortbl2G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("meanFreq", "cds_density"))]
cortbl2G[, vars := "meanFreq_cds_density"]

cortbl3G <- gnomG[ID==ID[1], gnom_cor_decomp(data = .SD, chromosome = "chr", signals = c("r", "cds_density"))]
cortbl3G[, vars := "r_cds_density"]

wavcorG <- rbindlist(list(cortbl1G, cortbl1G.1,cortbl2G,cortbl3G))
wavcorG[, units := "genetic"]

wavcor <- rbind(wavcorG, wavcorP)


# # ===== Linear Model analysis =====

# ----physical units
rsqrd1_p <- modwt_lm_rsqrd(data = gnomP[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
                      xvars = c("r"))
rsqrd1_p[, model := "r"]

rsqrd1.1_p <- modwt_lm_rsqrd(data = gnomP[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
                           xvars = c("log10r"))
rsqrd1.1_p[, model := "log10r"]

rsqrd2_p <- modwt_lm_rsqrd(data = gnomP[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
                         xvars = c("r", "cds_density"))
rsqrd2_p[, model := "r_cdsDensity"]

rsqrd2.1_p <- modwt_lm_rsqrd(data = gnomP[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
              xvars = c("log10r", "cds_density"))
rsqrd2.1_p[, model := "log10r_cdsDensity"]

rsqrd_p <- rbindlist(list(rsqrd1_p, rsqrd1.1_p, rsqrd2_p, rsqrd2.1_p))
rsqrd_p[, units := "physical"]



# ---- genetic units
rsqrd1_g <- modwt_lm_rsqrd(data = gnomG[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
                           xvars = c("r"))
rsqrd1_g[, model := "r"]

rsqrd1.1_g <- modwt_lm_rsqrd(data = gnomG[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
                             xvars = c("r"))
rsqrd1.1_g[, model := "log10r"]

rsqrd2_g <- modwt_lm_rsqrd(data = gnomG[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
                           xvars = c("r", "cds_density"))
rsqrd2_g[, model := "r_cdsDensity"]

rsqrd2.1_g <- modwt_lm_rsqrd(data = gnomG[ID==ID[1]], chromosome="chr", yvar = "meanFreq",
                             xvars = c("log10r", "cds_density"))
rsqrd2.1_g[, model := "log10r_cdsDensity"]

rsqrd_g <- rbindlist(list(rsqrd1_g, rsqrd1.1_g, rsqrd2_g, rsqrd2.1_g))
rsqrd_g[, units := "genetic"]


rsqrd <- rbind(rsqrd_p, rsqrd_g)

# # ====== Save output =======
save(wavvar, wavcor, rsqrd, file = paste0("ACUA_", year, "/wavelet_results.RData"))
