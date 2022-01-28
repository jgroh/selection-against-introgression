# 1. ========== Load Dependencies and Data ==========
library(data.table)
library(tools)
library(waveslim)
library(stringi)
#library(ggplot2kb
library(magrittr)

source("../../wavelet_functions.R")

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]
#year <- 2018

# read cds density file
cdsCM <- fread("xbir_1kb_CDS_and_cM.txt")

# read chrom lengths file
chrLenP <- fread("xbir10x_chrlengths.txt", col.names = c("chr", "len"))

# interpolated ancestry files for each chromosome
scaffFiles <- paste0("ACUA_2018/",paste0(chrLenP$chr, ".RData"))
# scaffFiles <- dir(path=paste0("ACUA_",year),pattern=".RData",full.names=T) # older version, useful if not all chromosome files present
# scaffFiles <- scaffFiles[1:3] # for testing locally
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
gnomG[,indivFreq := 1 - indivFreq]
gnomG[,meanFreq := 1 - meanFreq]

gnomP[,indivFreq := 1 - indivFreq]
gnomP[,meanFreq := 1 - meanFreq]

# merge ancestry data and genomic features data
cdsCM[, position := as.integer(start + 500)] # this is the midpoint of the 1kb intervals where I interpolate ancestry
gnomP <- merge(gnomP, cdsCM, by = c("chr", "position"))


# 2. ===== Data Transformation ===== 

gnomP[, position := seq_len(.N), by = .(ID,chr)]
gnomG[, position := seq_len(.N), by = .(ID,chr)]

gnomP[,meanFreq := mean(indivFreq), by = .(position,chr)]
gnomG[, meanFreq := mean(indivFreq), by = .(position,chr)]


# log transform recombination rates
gnomP[, cmTr := log10(cM)]

#ggplot(gnomP[sample(1:8e6,1000)], aes(x = log(coding_bp+1), y = meanFreqTr, colour = chr)) + 
#  geom_point() + geom_smooth(method = "lm")

#ggplot(gnomP[chr == chr[1] & position < 1000000], aes(x= position, y = coding_bp)) + geom_point()


# 3. ===== Total Variances =====

# 3.1. ----- mean ancestry -----
totalVar_meanFreqP <- gnomP[ID==ID[1],var(meanFreq)]
totalVar_meanFreqG <- gnomG[ID==ID[1],var(meanFreq)]

# 3.2. ----- individual-level ancestry -----
totalVar_indivFreqP <- gnomP[, .(totalVar = var(indivFreq)), by=ID]
totalVar_indivFreqG <- gnomG[, .(totalVar = var(indivFreq)), by=ID]

# 3.1. ----- recombination (log10 cM) -----
totalVar_cmTrP <- gnomP[ID==ID[1],var(cmTr)]

# ----- coding bp density -----
totalVar_cdsP <- gnomP[ID==ID[1],var(coding_bp)]


# 4. ===== Chromosome-level variance =====
# total genetic length as weights for chromosomes 
chrLenG <- gnomG[, .(weight = max(Morgan)), by = chr]
chrLenG[, weight := weight/sum(weight)]

# total physical length as weights
chrLenP[, weight := len/sum(len)]

# 4.1 ---------- mean ancestry, genetic units ----------

# chromosome means of population mean minor parent ancestry
chrMeans_meanFreqG <- gnomG[, .(meanFreq = mean(meanFreq)), by = .(chr)]
chrMeans_meanFreqG <- merge(chrMeans_meanFreqG, chrLenG)

# total mean ancestry
totalMean_meanFreqG <- gnomG[ID==ID[1], mean(meanFreq)]

# chromosome-level weighted variance of pop mean ancestry. Because weights sum to one, we just take weighted sum of squares
chrVar_meanFreqG <- chrMeans_meanFreqG[, sum(weight*(meanFreq - totalMean_meanFreqG)^2)]

#proportion of total genomic variance accounted for by variation among chromosome means
chrPropVar_meanFreqG <- chrVar_meanFreqG/totalVar_meanFreqG

#save(totalVarP, totalVarG, file = paste0("ACUA_",year,"/totalVar.RData"))


# 4.2. ----- individual-level ancestry, genetic units -----

# chromosome means per individual
chrMeans_indivFreqG <- gnomG[, .(chrMeanFreq = mean(indivFreq)), by = .(ID,chr)]
chrMeans_indivFreqG <- merge(chrMeans_indivFreqG, chrLenG)

# add total mean ancestry per individual
chrMeans_indivFreqG[, totalMean := weighted.mean(chrMeanFreq, weight), by = ID]

# chromosome-level weighted variance of ancestry per individual. Because weights sum to one, we just take weighted sum of squares
chrVar_indivFreqG <- chrMeans_indivFreqG[, .(chrVar = sum(weight*(chrMeanFreq - totalMean)^2)), by = ID]

#proportion of total genomic variance accounted for by variation among chromosome means
chrPropVar_indivFreqG <- merge(totalVar_indivFreqG, chrVar_indivFreqG)[, .(chrPropVar = chrVar/totalVar), by = ID]


# 4.3. ----- mean ancestry, physical units -----
# chromosome means of population mean minor parent ancestry
chrMeans_meanFreqP <- gnomP[, .(meanFreq = mean(meanFreq)), by = .(chr)]
chrMeans_meanFreqP <- merge(chrMeans_meanFreqP, chrLenP)

# total mean ancestry
totalMean_meanFreqP <- gnomP[ID==ID[1], mean(meanFreq)]

# chromosome-level weighted variance of pop mean ancestry. Because weights sum to one, we just take weighted sum of squares
chrVar_meanFreqP <- chrMeans_meanFreqP[, sum(weight*(meanFreq - totalMean_meanFreqP)^2)]

#proportion of total genomic variance accounted for by variation among chromosome means
chrPropVar_meanFreqP <- chrVar_meanFreqP/totalVar_meanFreqP

# 4.4. ----- individual-level ancestry, physical units -----
# chromosome means per individual
chrMeans_indivFreqP <- gnomP[, .(chrMeanFreq = mean(indivFreq)), by = .(ID,chr)]
chrMeans_indivFreqP <- merge(chrMeans_indivFreqP, chrLenP)

# add total mean ancestry per individual
chrMeans_indivFreqP[, totalMean := weighted.mean(chrMeanFreq, weight), by = ID]

# chromosome-level weighted variance of ancestry per individual. Because weights sum to one, we just take weighted sum of squares
chrVar_indivFreqP <- chrMeans_indivFreqP[, .(chrVar = sum(weight*(chrMeanFreq - totalMean)^2)), by = ID]

#proportion of total genomic variance accounted for by variation among chromosome means
chrPropVar_indivFreqP <- merge(totalVar_indivFreqP, chrVar_indivFreqP)[, .(chrPropVar = chrVar/totalVar), by = ID]


# 4.5 ----- recombination rate, physical units -----

# chromosome means 
chrMeans_cmTrP <- gnomP[, .(meanCmTr = mean(cmTr)), by = .(chr)]
chrMeans_cmTrP <- merge(chrMeans_cmTrP, chrLenP)

# total genomic mean recombination
totalMean_cmTrP <- gnomP[ID==ID[1], mean(cmTr)]

# chromosome-level weighted variance of recombination. Because weights sum to one, we just take weighted sum of squares
chrVar_cmTrP <- chrMeans_cmTrP[, sum(weight*(meanCmTr - totalMean_cmTrP)^2)]

#proportion of total genomic variance accounted for by variation among chromosome means
chrPropVar_cmTrP <- chrVar_cmTrP/totalVar_cmTrP


# 4.5 ----- coding bp, physical units -----

# chromosome means 
chrMeans_cdsP <- gnomP[, .(meanCds = mean(coding_bp)), by = .(chr)]
chrMeans_cdsP <- merge(chrMeans_cdsP, chrLenP)

# total genomic mean recombination
totalMean_cdsP <- gnomP[ID==ID[1], mean(coding_bp)]

# chromosome-level weighted variance of recombination. Because weights sum to one, we just take weighted sum of squares
chrVar_cdsP <- chrMeans_cdsP[, sum(weight*(meanCds - totalMean_cdsP)^2)]

#proportion of total genomic variance accounted for by variation among chromosome means
chrPropVar_cdsP <- chrVar_cdsP/totalVar_cdsP



# ===== Wavelet Transforms =====

# ----- mean ancestry, genetic scale -----
maxlevsG <- gnomG[ID==ID[1], floor(log2(nrow(.SD))), by = chr]$V1
allcolsG <- unique(c(paste0("d", 1:max(maxlevsG)), paste0("s", maxlevsG)))

# modwt of mean ancestry on genetic scale
meanModwtG <- gnomG[ID==ID[1], brickWallModwt(meanFreq, allcolsG), by = chr] %>%
  melt(., id.vars = "chr", variable.name = "level", value.name = "coefficient")
meanModwtG[, position:= seq_len(.N), by = chr]

# ----- individual-level ancetry, genetic scale -----

# modwt per individual
indivModwtG <- gnomG[, brickWallModwt(meanFreq, allcolsG), by = .(ID,chr)] %>%
  melt(., id.vars = c("ID","chr"), variable.name = "level", value.name = "coefficient")
indivModwtG[, position:= seq_len(.N), by = .(ID,chr)]


#  ----- mean ancestry, physical units -----

maxlevsP <- gnomP[ID==ID[1], floor(log2(nrow(.SD))), by = chr]$V1
allcolsP <- unique(c(paste0("d", 1:max(maxlevsP)), paste0("s", maxlevsP)))

# modwt of mean ancestry
meanModwtP <- gnomP[ID==ID[1], brickWallModwt(meanFreq, allcolsP), by = chr] %>%
  melt(., id.vars = "chr", variable.name = "level", value.name = "coefficient")
meanModwtP[, position:= seq_len(.N), by = chr]


# 5.5. ----- individual-level ancestry, physical units -----

# modwt per individual
indivModwtP <- gnomP[, brickWallModwt(meanFreq, allcolsP), by = .(ID,chr)] %>%
  melt(., id.vars = c("ID","chr"), variable.name = "level", value.name = "coefficient")
indivModwtP[, position:= seq_len(.N), by = .(ID,chr)]


# ----- Recombination rate, physical units -------

cmTrModwtP <- gnomP[ID==ID[1], brickWallModwt(cmTr, allcolsP), by = chr] %>%
  melt(., id.vars = "chr", variable.name = "level", value.name = "coefficient")
cmTrModwtP[, position:= seq_len(.N), by = chr]

# ----- Recombination rate, physical units -------

cdsModwtP <- gnomP[ID==ID[1], brickWallModwt(coding_bp, allcolsP), by = chr] %>%
  melt(., id.vars = "chr", variable.name = "level", value.name = "coefficient")
cdsModwtP[, position:= seq_len(.N), by = chr]


# 5. ========== Wavelet variances ==========

# 5.1 ----- mean ancestry, genetic scale -----

# wavelet variance estimates from modwt
meanWavVarG_Chrs <- gnomG[ID==ID[1], waveletVarianceModwt(meanFreq, allcolsG, na.condition="na"), 
                        by =  chr] %>%
  melt(., id.vars = "chr", variable.name = "level", value.name = "variance")

# weighted average variance magnitudes per scale across chromosomes (exclude chromosomes without a certain scale present)
chrWaveletWeightsG <- meanModwtG[!is.na(coefficient), .(n.coeff=nrow(.SD)), by = .(chr, level)]
chrWaveletWeightsG[, wavelet.weight:= n.coeff/sum(n.coeff), by = level][, n.coeff:=NULL][]
meanWavVarG_Chrs <- merge(meanWavVarG_Chrs, chrWaveletWeightsG)

meanWavVarG_magnitude <- meanWavVarG_Chrs[, .(variance = weighted.mean(variance, wavelet.weight, na.rm=T)), by = level]

# proportion of genomic variance by scale per chromosome 
meanWavVarG_Chrs[is.na(variance), variance := 0] # reassign NA values to zero
meanWavVarG_Chrs[, propVar := variance/sum(variance), by = chr]

# weight chromosomes by total genetic length
meanWavVarG_Chrs <- merge(chrLenG, meanWavVarG_Chrs, by = "chr")

# average proportion over chromosomes, weighting by length
meanWavVarG <- meanWavVarG_Chrs[, .(propVar = weighted.mean(propVar, weight)), by = level]
#ggplot(meanWavVarG, aes(x = level, y = propVar)) + geom_point()

# adjust proportions to account for additional chromosome-level variance
meanWavVarG[, propVar := propVar*(1-chrPropVar_meanFreqG)]

# total genomic proportion of ancestry variance
propVarDecomp_meanFreqG <- rbind(meanWavVarG, data.table(level = "chr", propVar = chrPropVar_meanFreqG))

# should sum to 1
#propVarDecomp_meanFreqG[, sum(propVar)]
#ggplot(propVarDecomp_meanFreqG, aes(level, propVar))+geom_point()

# 5.1.1 ------ combine tables for both raw variance decomposition and proportion of genomic variance decomosition ----- 

rawVarDecomp_meanFreqG <- rbind(meanWavVarG_magnitude, data.table(level = "chr", variance = chrVar_meanFreqG))
rawVarDecomp_meanFreqG[, signal := "mean"]
#ggplot(rawVarDecomp_meanFreqG, aes(level, variance)) + geom_point()

# combine genomic contribution variance decomposition and variance magnitude decomposition
meanFreqAllVarDecompG <- merge(propVarDecomp_meanFreqG, rawVarDecomp_meanFreqG)


# 5.2. ----- individual-level ancetry, genetic scale -----

# wavelet variance estimates from modwt
indivWavVarG_Chrs <- gnomG[, waveletVarianceModwt(indivFreq, allcolsG, na.condition="na"), 
                          by =  .(ID,chr)] %>%
  melt(., id.vars = c("ID","chr"), variable.name = "level", value.name = "variance")

# average over chromosomes within individual using weights of how many wavelets per chromosome
# (exclude chromosomes without a certain scale present)
indivWavVarG_Chrs <- merge(indivWavVarG_Chrs, chrWaveletWeightsG)
indivWavVarG_magnitude <- indivWavVarG_Chrs[, .(variance= weighted.mean(variance, wavelet.weight, na.rm=T)), by= .(ID,level)]

# average over individuals
indivWavVarG_magnitude <- indivWavVarG_magnitude[, .(variance = mean(variance)), by = level]

# proportion of variance by scale per chromosome per individual
indivWavVarG_Chrs[is.na(variance), variance := 0] # reassign NA values to zero
indivWavVarG_Chrs[, propVar := variance/sum(variance), by = .(ID,chr)]

# weight chromosomes by total genetic length
indivWavVarG_Chrs <- merge(chrLenG, indivWavVarG_Chrs, by=  "chr")

# average proportion over chromosomes within individuals
indivWavVarG <- indivWavVarG_Chrs[, .(propVar = weighted.mean(propVar, weight)), by = .(ID,level)]
#ggplot(meanWavVarG, aes(x = level, y = propVar)) + geom_point()


# adjust proportions to account for additional chromosome-level variance
indivWavVarG <- merge(chrPropVar_indivFreqG, indivWavVarG)
indivWavVarG[, propVar := propVar*(1-chrPropVar), by = ID]
indivWavVarG[, chrPropVar := NULL]

chrPropVar_indivFreqG[, level := "chr"]
setnames(chrPropVar_indivFreqG, "chrPropVar", "propVar")

# total genomic proportion of ancestry variance
propVarDecomp_indivFreqG <- rbind(indivWavVarG, chrPropVar_indivFreqG)

# should sum to 1 within each individual
#propVarDecomp_indivFreqG[, sum(propVar), by = ID]

# average over individuals
propVarDecomp_indivFreqG <- propVarDecomp_indivFreqG[, .(propVar = mean(propVar)), by = level]
propVarDecomp_indivFreqG[, signal := "individual"]

# 5.2.1 ------ combine tables for both raw variance decomposition and proportion of genomic variance decomosition ----- 

rawVarDecomp_indivFreqG <- rbind(indivWavVarG_magnitude, data.table(level = "chr", variance = chrVar_indivFreqG[, mean(chrVar)]))
rawVarDecomp_indivFreqG[, signal := "individual"]

#ggplot(rbind(rawVarDecomp_indivFreqG, rawVarDecomp_meanFreqG), aes(level, variance, group = signal)) + geom_line()

# combine genomic contribution variance decomposition and variance magnitude decomposition
indivFreqAllVarDecompG <- merge(propVarDecomp_indivFreqG, rawVarDecomp_indivFreqG)


# 5.3. ---- combine mean ancestry and individual-level ancestry tables, genetic units ----

allVarDecompG <- rbind(indivFreqAllVarDecompG, meanFreqAllVarDecompG)

#ggplot(allVarDecompG, aes(level, variance, group = signal, color = signal)) + geom_point() + geom_line()




# 5.4 ----- mean ancestry, physical units -----

# wavelet variance estimates from modwt
meanWavVarP_Chrs <- gnomP[ID==ID[1], waveletVarianceModwt(meanFreq, allcolsP, na.condition="na"), 
                          by =  chr] %>%
  melt(., id.vars = "chr", variable.name = "level", value.name = "variance")

# weighted average variance magnitudes per scale across chromosomes (exclude chromosomes without a certain scale present)
chrWaveletWeightsP <- meanModwtP[!is.na(coefficient), .(n.coeff=nrow(.SD)), by = .(chr, level)]
chrWaveletWeightsP[, wavelet.weight:= n.coeff/sum(n.coeff), by = level][, n.coeff:=NULL][]
meanWavVarP_Chrs <- merge(meanWavVarP_Chrs, chrWaveletWeightsP, all=T) # all=T important so that NA values are not removed

# weighted average over chromosomes
meanWavVarP_magnitude <- meanWavVarP_Chrs[, .(variance = weighted.mean(variance, wavelet.weight, na.rm=T)), by = level]

# proportion of genomic variance by scale per chromosome 
meanWavVarP_Chrs[is.na(variance), variance := 0] # reassign NA values to zero
meanWavVarP_Chrs[, propVar := variance/sum(variance), by = chr]

# weight chromosomes by total genetic length
meanWavVarP_Chrs <- merge(chrLenP, meanWavVarP_Chrs, by = "chr")

# average proportion over chromosomes, weighting by length
meanWavVarP <- meanWavVarP_Chrs[, .(propVar = weighted.mean(propVar, weight)), by = level]
#ggplot(meanWavVarP, aes(x = level, y = propVar)) + geom_point()

# adjust proportions to account for additional chromosome-level variance
meanWavVarP[, propVar := propVar*(1-chrPropVar_meanFreqP)]

# total genomic proportion of ancestry variance
propVarDecomp_meanFreqP <- rbind(meanWavVarP, data.table(level = "chr", propVar = chrPropVar_meanFreqP))

# should sum to 1
#propVarDecomp_meanFreqP[, sum(propVar)]
#ggplot(propVarDecomp_meanFreqP, aes(level, propVar))+geom_point()

# 5.4.1 ------ combine tables for both raw variance decomposition and proportion of genomic variance decomosition ----- 

rawVarDecomp_meanFreqP <- rbind(meanWavVarP_magnitude, data.table(level = "chr", variance = chrVar_meanFreqP))
rawVarDecomp_meanFreqP[, signal := "mean"]
#ggplot(rawVarDecomp_meanFreqP, aes(level, variance)) + geom_point()

# combine genomic contribution variance decomposition and variance magnitude decomposition
meanFreqAllVarDecompP <- merge(propVarDecomp_meanFreqP, rawVarDecomp_meanFreqP)


# 5.5. ----- individual-level ancestry, physical units -----

# wavelet variance estimates from modwt
indivWavVarP_Chrs <- gnomP[, waveletVarianceModwt(indivFreq, allcolsP, na.condition="na"), 
                           by =  .(ID,chr)] %>%
  melt(., id.vars = c("ID","chr"), variable.name = "level", value.name = "variance")

# average over chromosomes within individual using weights of how many wavelets per chromosome
# (exclude chromosomes without a certain scale present)
indivWavVarP_Chrs <- merge(indivWavVarP_Chrs, chrWaveletWeightsP)
indivWavVarP_magnitude <- indivWavVarP_Chrs[, .(variance= weighted.mean(variance, wavelet.weight, na.rm=T)), by= .(ID,level)]

# average over individuals
indivWavVarP_magnitude <- indivWavVarP_magnitude[, .(variance = mean(variance)), by = level]

# proportion of variance by scale per chromosome per individual
indivWavVarP_Chrs[is.na(variance), variance := 0] # reassign NA values to zero
indivWavVarP_Chrs[, propVar := variance/sum(variance), by = .(ID,chr)]

# weight chromosomes by total genetic length
indivWavVarP_Chrs <- merge(chrLenP, indivWavVarP_Chrs, by=  "chr")

# average proportion over chromosomes within individuals
indivWavVarP <- indivWavVarP_Chrs[, .(propVar = weighted.mean(propVar, weight)), by = .(ID,level)]
#ggplot(meanWavVarP, aes(x = level, y = propVar)) + geom_point()


# adjust proportions to account for additional chromosome-level variance
indivWavVarP <- merge(chrPropVar_indivFreqP, indivWavVarP)
indivWavVarP[, propVar := propVar*(1-chrPropVar), by = ID]
indivWavVarP[, chrPropVar := NULL]

chrPropVar_indivFreqP[, level := "chr"]
setnames(chrPropVar_indivFreqP, "chrPropVar", "propVar")

# total genomic proportion of ancestry variance
propVarDecomp_indivFreqP <- rbind(indivWavVarP, chrPropVar_indivFreqP)

# should sum to 1 within each individual
#propVarDecomp_indivFreqP[, sum(propVar), by = ID]

# average over individuals
propVarDecomp_indivFreqP <- propVarDecomp_indivFreqP[, .(propVar = mean(propVar)), by = level]
propVarDecomp_indivFreqP[, signal := "individual"]

# 5.5.1 ------ combine tables for both raw variance decomposition and proportion of genomic variance decomosition ----- 

rawVarDecomp_indivFreqP <- rbind(indivWavVarP_magnitude, data.table(level = "chr", variance = chrVar_indivFreqP[, mean(chrVar)]))
rawVarDecomp_indivFreqP[, signal := "individual"]

#ggplot(rbind(rawVarDecomp_indivFreqP, rawVarDecomp_meanFreqP), aes(level, variance, group = signal)) + geom_line()

# combine genomic contribution variance decomposition and variance magnitude decomposition
indivFreqAllVarDecompP <- merge(propVarDecomp_indivFreqP, rawVarDecomp_indivFreqP)


# 5.6. ---- combine mean ancestry and individual-level ancestry tables, physical units ----

allVarDecompP <- rbind(indivFreqAllVarDecompP, meanFreqAllVarDecompP)

#ggplot(allVarDecompP, aes(level, propVar, group = signal, color = signal)) + geom_point() + geom_line()





#  ----- recombination rate, physical units =-----

# wavelet variance estimates from modwt
cmTrWavVarP_Chrs <- gnomP[ID==ID[1], waveletVarianceModwt(cmTr, allcolsP, na.condition="na"), 
                          by =  chr] %>%
  melt(., id.vars = "chr", variable.name = "level", value.name = "variance")

# weighted average variance magnitudes per scale across chromosomes (exclude chromosomes without a certain scale present)

cmTrWavVarP_Chrs <- merge(cmTrWavVarP_Chrs, chrWaveletWeightsP, all=T) # all=T important so that NA values are not removed

# weighted average over chromosomes
cmTrWavVarP_magnitude <- cmTrWavVarP_Chrs[, .(variance = weighted.mean(variance, wavelet.weight, na.rm=T)), by = level]

# proportion of genomic variance by scale per chromosome 
cmTrWavVarP_Chrs[is.na(variance), variance := 0] # reassign NA values to zero
cmTrWavVarP_Chrs[, propVar := variance/sum(variance), by = chr]

# weight chromosomes by total genetic length
cmTrWavVarP_Chrs <- merge(chrLenP, cmTrWavVarP_Chrs, by = "chr")

# average proportion over chromosomes, weighting by length
cmTrWavVarP <- cmTrWavVarP_Chrs[, .(propVar = weighted.mean(propVar, weight)), by = level]
#ggplot(cmTrWavVarP, aes(x = level, y = propVar)) + geom_point()

# adjust proportions to account for additional chromosome-level variance
cmTrWavVarP[, propVar := propVar*(1-chrPropVar_cmTrP)]

# total genomic proportion of ancestry variance
propVarDecomp_cmTrP <- rbind(cmTrWavVarP, data.table(level = "chr", propVar = chrPropVar_cmTrP))

# should sum to 1
#propVarDecomp_cmTrP[, sum(propVar)]
#ggplot(propVarDecomp_cmTrP, aes(level, propVar))+geom_point()


# 5.4.1 ------ combine tables for both raw variance decomposition and proportion of genomic variance decomosition ----- 

rawVarDecomp_cmTrP <- rbind(cmTrWavVarP_magnitude, data.table(level = "chr", variance = chrVar_cmTrP))
rawVarDecomp_cmTrP[, signal := "mean"]
#ggplot(rawVarDecomp_meanFreqP, aes(level, variance)) + geom_point()

# combine genomic contribution variance decomposition and variance magnitude decomposition
cmTrAllVarDecompP <- merge(propVarDecomp_cmTrP, rawVarDecomp_cmTrP)
cmTrAllVarDecompP[, units := "physical"]
cmTrAllVarDecompP[, variable := "log10cM"]





#  ----- cds rate, physical units =-----

# wavelet variance estimates from modwt
cdsWavVarP_Chrs <- gnomP[ID==ID[1], waveletVarianceModwt(coding_bp, allcolsP, na.condition="na"), 
                          by =  chr] %>%
  melt(., id.vars = "chr", variable.name = "level", value.name = "variance")

# weighted average variance magnitudes per scale across chromosomes (exclude chromosomes without a certain scale present)
cdsWavVarP_Chrs <- merge(cdsWavVarP_Chrs, chrWaveletWeightsP, all=T) # all=T important so that NA values are not removed

# weighted average over chromosomes
cdsWavVarP_magnitude <- cdsWavVarP_Chrs[, .(variance = weighted.mean(variance, wavelet.weight, na.rm=T)), by = level]

# proportion of genomic variance by scale per chromosome 
cdsWavVarP_Chrs[is.na(variance), variance := 0] # reassign NA values to zero
cdsWavVarP_Chrs[, propVar := variance/sum(variance), by = chr]

# weight chromosomes by total genetic length
cdsWavVarP_Chrs <- merge(chrLenP, cdsWavVarP_Chrs, by = "chr")

# average proportion over chromosomes, weighting by length
cdsWavVarP <- cdsWavVarP_Chrs[, .(propVar = weighted.mean(propVar, weight)), by = level]
#ggplot(cmTrWavVarP, aes(x = level, y = propVar)) + geom_point()

# adjust proportions to account for additional chromosome-level variance
cdsWavVarP[, propVar := propVar*(1-chrPropVar_cdsP)]

# total genomic proportion of ancestry variance
propVarDecomp_cdsP <- rbind(cdsWavVarP, data.table(level = "chr", propVar = chrPropVar_cdsP))

# should sum to 1
#propVarDecomp_cdsP[, sum(propVar)]
#ggplot(propVarDecomp_cdsP, aes(level, propVar))+geom_point()




# 5.4.1 ------ combine tables for both raw variance decomposition and proportion of genomic variance decomosition ----- 

rawVarDecomp_cdsP <- rbind(cdsWavVarP_magnitude, data.table(level = "chr", variance = chrVar_cdsP))
rawVarDecomp_cdsP[, signal := "mean"]

# combine genomic contribution variance decomposition and variance magnitude decomposition
cdsAllVarDecompP <- merge(propVarDecomp_cdsP, rawVarDecomp_cdsP)
cdsAllVarDecompP[, units := "physical"]
cdsAllVarDecompP[, variable := "CDS_density"]


# 5.7 ----- combine all ancestry variance decompositions -----=
allVarDecompG[, units := "genetic"]
allVarDecompP[, units := "physical"]

allVarDecomp <- rbind(allVarDecompG, allVarDecompP)
allVarDecomp[, variable := "ancestry"]

allVarDecomp_anc_rec_cds <- rbind(allVarDecomp, cmTrAllVarDecompP, cdsAllVarDecompP)

save(allVarDecomp_anc_rec_cds, file =  paste0("ACUA_",year,"/allVarDecomp_anc_rec_cds.RData"))




# ===== Correlation Analysis ==== 
meanModwtP2 <- meanModwtP
meanModwtP2[, ancestry_coeff := coefficient]
cmTrModwtP2 <- cmTrModwtP
cmTrModwtP2[, rec_coeff := coefficient]
cdsModwtP2 <- cdsModwtP
cdsModwtP2[, cds_coeff := coefficient]



allModwtP <- merge(meanModwtP2[, .(chr,level,position,ancestry_coeff)], 
      merge(cmTrModwtP2[,.(chr,level,position,rec_coeff)],
            cdsModwtP2[,.(chr,level,position,cds_coeff)]),
      by = c("chr","level","position"))





# 8.2. --------- Wavelet covariances -----------

# unweighted wavelet covariances at each scale separately for each chromosome
recAncCov <- allModwtP[, .(rec_anc_cov = sum(rec_coeff*ancestry_coeff, na.rm=T)/length(ancestry_coeff[!is.na(ancestry_coeff)])), by = .(chr, level)]

cdsAncCov <- allModwtP[, .(cds_anc_cov = sum(cds_coeff*ancestry_coeff, na.rm=T)/length(ancestry_coeff[!is.na(ancestry_coeff)])), by = .(chr, level)]

cdsRecCov <- allModwtP[, .(cds_rec_cov = sum(cds_coeff*rec_coeff, na.rm=T)/length(ancestry_coeff[!is.na(ancestry_coeff)])), by = .(chr, level)]



# combine covariance values
allVarsCov <- merge(merge(recAncCov, cdsAncCov), cdsRecCov)

# Weighted average of covariances over chromosomes
allVarsCov <- merge(allVarsCov, chrWaveletWeightsP)[, lapply(.SD, weighted.mean, w=wavelet.weight),
                               .SDcols = c("rec_anc_cov","cds_anc_cov",
                                           "cds_rec_cov"),
                               by = level]

# merge covariance table with wavelet variance table

allVarDecompP[, anc_var := variance][, variance := NULL]
cmTrAllVarDecompP[, rec_var := variance][, variance := NULL]
cdsAllVarDecompP[, cds_var := variance][, variance := NULL]

allVarsWavVarCov <- merge(allVarsCov, 
                          merge(merge(allVarDecompP[signal == "mean", .(level,anc_var)],
                     cmTrAllVarDecompP[,.(level,rec_var)],by="level"),
               cdsAllVarDecompP[,.(level,cds_var)], by = "level"), by = "level")

 
# 8.3. ----------- Compute Wavelet Correlations between signal pairs ------------

recAncWavCor <- allVarsWavVarCov[, .(rec_anc_cor = rec_anc_cov/(sqrt(anc_var)*sqrt(rec_var))), by = level]

cdsAncWavCor <- allVarsWavVarCov[, .(cds_anc_cor = cds_anc_cov/(sqrt(cds_var)*sqrt(anc_var))), by = level]

cdsRecWavCor <- allVarsWavVarCov[, .(cds_rec_cor = cds_rec_cov/(sqrt(cds_var)*sqrt(rec_var))), by = level]

# combine wavelet correlation tables
allVarsWavCor <-  merge(merge(recAncWavCor, cdsAncWavCor), cdsRecWavCor, by = "level")

allCorDecomp_anc_rec_cds <- 
  rbind(allVarsWavCor, 
      data.table(level = "chr",
                 rec_anc_cor = cor(chrMeans_cmTrP$meanCmTr, chrMeans_meanFreqP$meanFreq),
                 cds_anc_cor = cor(chrMeans_cdsP$meanCds, chrMeans_meanFreqP$meanFreq),
                 cds_rec_cor = cor(chrMeans_cmTrP$meanCmTr, chrMeans_cdsP$meanCds)))

save(allCorDecomp_anc_rec_cds, file =  paste0("ACUA_",year,"/allCorDecomp_anc_rec_cds.RData"))

#ggplot(allVarsWavCor[grepl("d",level,fixed=T)], aes(x = level, y = rec_anc_cor)) + geom_point()








# 
# 
# 
# 
# dwt.nondyadic <- function(x){
#   # discrete wavelet transform for signal that is not length power of 2
#   # (allows for maximum number of coefficients at each scale w/o extending past end of chrom)
#   x <- x-mean(x)
#   M <- length(x)
#   N <- ceiling(log(M, 2)) # next highest power of 2
#   J <- floor(log2(length(x))) # max power of 2
#   xx <- c(x, rep(0, 2^N - M)) # append zeros to end of sequence to reach length N
#   y <- dwt(xx, wf = "haar", n.levels = N) # a list with coefficients for each scale
#   
#   for(j in 1:(length(y))){
#     if(j <= J){
#       y[[j]] <- y[[j]][1:floor(M/2^j)] # truncates to only those coefficients not affected by appending zeros to end
#     } 
#   }
#   y[(J+1):length(y)] <- NULL # also changed to simply remove other coefficients
# 
#   return(y)
# }
# 
# # To get number of non-boundary coefficients
# numCoeff <- function(x){length(x[!is.na(x)])} 
# 
# # Unbiased estimator of wavelet variance for MODWT with brick wall boundary condition
# wav_var <- function(x){sum(x^2,na.rm=TRUE)/(length(x[!is.na(x)]))} 
# 
# # To compute dwt on variable in a data table, for a single group
# dwtAllScales <- function(x,variable,allcols){
#   y <- x[, dwt.nondyadic(variable)]
#   dt <- as.data.table(stri_list2matrix(y))
#   setnames(dt, names(y))
#   
#   # add empty detail columns for higher scales not present on particular chromosome
#   naCols <- setdiff(allcols, names(dt))
#   if(length(naCols) > 0){ dt[,(naCols) := as.numeric(NA)] }
#   return(dt[,lapply(.SD, as.double)])
# }
# 
# # 4. ============ DWT: Ancestry, recombination, and coding bp ==========
# # define levels based on longest chromosome
# maxLevelP <- max(gnomP[ID==ID[1],floor(log2(length(position))),by=chr][,2])
# allColsP <- paste0("d",1:maxLevelP)
# 
# dwtAnc <- gnomP[ID==ID[1], dwtAllScales(.SD, variable = meanFreq, allcols = allColsP), by = chr]
# dwtRec <- gnomP[ID==ID[1], dwtAllScales(.SD, variable = cmTr, allcols = allColsP), by = chr]
# dwtCds <- gnomP[ID==ID[1], dwtAllScales(.SD, variable = coding_bp,  allcols = allColsP), by = chr]
# 
# # combine tables
# lapply(list(dwtAnc, dwtRec, dwtCds), function(x){ x[,position := seq_len(.N), by = chr] })
# 
# dwtAnc <- melt(dwtAnc, id.vars = c("chr","position"), 
#                variable.name = "scale",
#                value.name = "ancestry_coeff")
# 
# dwtRec <- melt(dwtRec, id.vars = c("chr","position"), 
#                variable.name = "scale",
#                value.name = "rec_coeff")
# dwtCds <- melt(dwtCds, id.vars = c("chr","position"), 
#                variable.name = "scale",
#                value.name = "cds_coeff")
# 
# allDWT <- Reduce(function(...) merge(..., all=TRUE), list(dwtAnc, dwtRec, dwtCds))
# save(allDWT, file = paste0("ACUA_",year,"/anc_rec_cds_DWT.RData"))
# 

# # 7.2. --------- MODWT on Mean Ancestry ----------
# popModwtP <- gnomP[ID==ID[1], modwtAllScales(.SD,variable=meanFreq,lenCol=position,allcols=allColsP), by = chr]
# popWavVarP <- popModwtP[,lapply(.SD,wav_var),by=chr]
# setnames(popWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
# popWavVarP <- melt(popWavVarP, measure.vars = as.character(1:maxLevelP),
#                          variable.name = "scale", value.name = "anc_variance")
# 
# # 5.3. ----- MODWT on Recomb Rate (cM length of 1kb windows) ----------
# recModwtP <- gnomP[ID==ID[1],modwtAllScales(.SD,variable=cmTr,lenCol=position,allcols=allColsP), by = chr]
# 
# # wavelet variance for rec rate
# recWavVarP <- recModwtP[,lapply(.SD,wav_var),by=chr]
# setnames(recWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
# recWavVarP <- melt(recWavVarP, measure.vars = as.character(1:maxLevelP),
#                    variable.name = "scale", value.name = "rec_variance")
# 
# # 5.4. ----- MODWT on # of coding bp per 1kb window ----------
# cdsModwtP <- gnomP[ID==ID[1],modwtAllScales(.SD,variable=coding_bp,lenCol=position,allcols=allColsP), by = chr]
# 
# # wavelet variance for CDS
# cdsWavVarP <- cdsModwtP[,lapply(.SD,wav_var),by=chr]
# setnames(cdsWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
# cdsWavVarP <- melt(cdsWavVarP, measure.vars = as.character(1:maxLevelP),
#                    variable.name = "scale", value.name = "cds_variance")
# 
# # 5.5. ----- MODWT on CDS per cM coding bp --------------
# cdsPerCmModwtP <- gnomP[ID==ID[1],modwtAllScales(.SD,variable=coding_bp/cmTr,lenCol=position,allcols=allColsP), by = chr]
# 
# # wavelet variance for cds density
# cdsPerCmWavVarP <- cdsPerCmModwtP[,lapply(.SD,wav_var),by=chr]
# setnames(cdsPerCmWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
# cdsPerCmWavVarP <- melt(cdsPerCmWavVarP, measure.vars = as.character(1:maxLevelP),
#                         variable.name = "scale", value.name = "cdsPerCm_variance")
# 
# # 5.6. ------- Combine tables with wavelet variances for all signals -------
# # save object with chromosome-specific variances of signals
# allVarsWavVar_Chrs <- merge(Reduce(merge, list(popWavVarP, recWavVarP, cdsWavVarP, cdsPerCmWavVarP)),
#                             chrWeightsP)
# allVarsWavVar_Chrs[,`:=`(distance = "physical", 
#                          decomp = "pop_mean")]
# 
# 
# # weighted average over chromosomes
# signalVars <- c("anc_variance", "rec_variance", "cds_variance", "cdsPerCm_variance")
# allVarsWavVar <- allVarsWavVar_Chrs[, lapply(.SD, weighted.mean, w=weight), .SDcols = signalVars, by = .(scale,distance,decomp)]
# 
# 
# # 6. ========== Chromosome-Level Analysis: Physical Scale ==========
# 
# # total physical length as weights for chromosomes 
# chrLenP <- gnomP[, max(position), by = .(chr)]
# setnames(chrLenP, "V1", "len")
# chrLenP[, weight := len/sum(len)]
# 
# # 6.1. ---------- Chromosome-level variance of mean ancestry, rec, cds ----------
# signalCols <- c("meanFreq", "cmTr", "coding_bp")
# 
# # computed chromosome means
# chrSignalMeans <- gnomP[, lapply(.SD, mean), .SDcols = signalCols, by = .(chr)]
# chrSignalMeans[, weight := chrLenP[,weight]]
# 
# # ----- save chromosome means -----
# save(chrSignalMeans, file = paste0("ACUA_",year,"/chrSignalMeans.RData"))
# 
# 
# # compute weighted variance 
# chrVarP <- chrSignalMeans[, lapply(.SD, function(x){sum(weight*(x - weighted.mean(x, weight))^2)}),
#                                           .SDcols = signalCols]
# setnames(chrVarP, c("anc_variance", "rec_variance", "cds_variance"))
# chrVarP[, `:=`(decomp = "pop_mean", 
#                distance = "physical",
#                scale = "chr")]
# 
# # 6.2. ---------- Ancestry: Individual-Level ----------
# 
# # average minor parent ancestry for individuals for each chrom
# chrMeansIndP <- gnomP[, mean(indivFreq), by = .(ID,chr)]
# setnames(chrMeansIndP, "V1", "avg_frq")
# 
# # calculate genome-wide mean for individuals by weighting chromosomes
# chrMeansIndP <- merge(chrMeansIndP, chrLenP, by = "chr") # add chrom weights
# chrMeansIndP[, gnomWideMean := weighted.mean(avg_frq,weight), by = ID]
# 
# # good place to look at distribution of genome-wide admixture proportions
# #boxplot(chrMeansIndP$gnomWideMean)
# 
# # chromosome-level weighted variance by individual, then average over individuals
# chrVarIndP <- chrMeansIndP[, sum(weight*(avg_frq - gnomWideMean)^2), by = ID][,mean(V1)]
# 
# 
# # 6.3 ------ Combine Chr-level Variances for individual and population ----------
# chrVarP <- rbind(list(anc_variance = chrVarIndP, 
#                       rec_variance = NA, 
#                       cds_variance = NA, 
#                       decomp = "mean_individual", 
#                       distance = "physical", 
#                       scale = "chr"), chrVarP)
# 
# 
# # 7. ========== Output Wavelet and Chrom-Level Variance Data ==========
# 
# # Combine Wavelet variance tables
# wvFinalG <- rbind(indMeanWavVarG, popMeanWavVarG)
# wvFinalG[, distance := "genetic"] 
# wvFinalP <- merge(indMeanWavVarP, allVarsWavVar, by = c("scale", "decomp", "distance", "anc_variance"), all=T)
# wvFinalAll <- merge(wvFinalP, wvFinalG, by = c("scale", "decomp", "distance", "anc_variance"), all =T)
# 
# # combine wavelet variances by chromosome (using mean ancestry)
# wvChrsAll <- merge(popWavVarG_Chrs, allVarsWavVar_Chrs, by = c("chr", "scale", "distance", "anc_variance"), all = T)
# 
# # combine chromosome-level variances
# chrVarAll <- merge(chrVarP, chrVarG, by = c("scale","distance", "anc_variance", "decomp"),all=T)
# 
# save(wvFinalAll, chrVarAll, wvChrsAll, file = paste0("ACUA_",year,"/anc_rec_cds_varDecomp.RData"))
# 
# 
# # 8. ========== Correlation Analysis ==========
# 
# # 8.1. ---------- Reformat MODWT tables ----------
# popModwtP[, position := seq_len(.N), by = chr]
# setnames(popModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
# popModwtP <- melt(popModwtP, measure.vars = as.character(1:maxLevelP),
#                   variable.name = "scale", value.name = "ancestry_coeff")
# 
# recModwtP[, position := seq_len(.N), by = chr]
# setnames(recModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
# recModwtP <- melt(recModwtP, measure.vars = as.character(1:maxLevelP),
#                    variable.name = "scale", value.name = "rec_coeff")
# 
# cdsModwtP[, position := seq_len(.N), by = chr]
# setnames(cdsModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
# cdsModwtP <- melt(cdsModwtP, measure.vars = as.character(1:maxLevelP),
#                   variable.name = "scale", value.name = "cds_coeff")
# 
# cdsPerCmModwtP[, position := seq_len(.N), by = chr]
# setnames(cdsPerCmModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
# cdsPerCmModwtP <- melt(cdsPerCmModwtP, measure.vars = as.character(1:maxLevelP),
#                   variable.name = "scale", value.name = "cdsPerCm_coeff")
# # combine
# allVarsModwtP <- Reduce(merge, list(popModwtP, recModwtP, cdsModwtP, cdsPerCmModwtP))
#  
# # save modwt coefficients for linear model analysis
# save(allVarsModwtP, file = paste0("ACUA_",year,"/anc_rec_cds_MODWT.RData"))
# 
# # 8.2. --------- Wavelet covariances -----------
# 
# # unweighted wavelet covariances at each scale separately for each chromosome
# recAncCov <- allVarsModwtP[, sum(rec_coeff*ancestry_coeff, na.rm=T)/length(ancestry_coeff[!is.na(ancestry_coeff)]), by = .(chr, scale)]
# setnames(recAncCov, "V1", "rec_anc_covariance")
# 
# cdsAncCov <- allVarsModwtP[, sum(cds_coeff*ancestry_coeff, na.rm=T)/length(ancestry_coeff[!is.na(ancestry_coeff)]), by = .(chr, scale)]
# setnames(cdsAncCov, "V1", "cds_anc_covariance")
# 
# cdsRecCov <- allVarsModwtP[, sum(cds_coeff*rec_coeff, na.rm=T)/length(ancestry_coeff[!is.na(ancestry_coeff)]), by = .(chr, scale)]
# setnames(cdsRecCov, "V1", "cds_rec_covariance")
# 
# cdsPerCmAncCov <- allVarsModwtP[, sum(cdsPerCm_coeff*ancestry_coeff, na.rm=T)/length(cds_coeff[!is.na(cdsPerCm_coeff)]), by = .(chr, scale)]
# setnames(cdsPerCmAncCov, "V1", "cdsPerCm_anc_covariance")
# 
# # combine covariance values
# allVarsCov <- Reduce(merge, list(recAncCov, cdsAncCov, cdsRecCov, cdsPerCmAncCov))
# 
# # Weighted average of covariances over chromosomes
# allVarsCov <- merge(allVarsCov, chrWeightsP)[, lapply(.SD, weighted.mean, w=weight), 
#                                .SDcols = c("rec_anc_covariance","cds_anc_covariance",
#                                            "cds_rec_covariance","cdsPerCm_anc_covariance"),
#                                by = scale]
# 
# # merge covariance table with wavelet variance table
# allVarsWavVarCov <- merge(allVarsWavVar, allVarsCov)
# 
# 
# # 8.3. ----------- Compute Wavelet Correlations between signal pairs ------------
# recAncWavCor <- allVarsWavVarCov[, rec_anc_covariance/(sqrt(rec_variance)*sqrt(anc_variance)), by = scale]
# setnames(recAncWavCor, "V1", "rec_anc_cor")
# #plot(y = recAncWavCor$rec_anc_cor, x=1:14)
# 
# cdsAncWavCor <- allVarsWavVarCov[, cds_anc_covariance/(sqrt(cds_variance)*sqrt(anc_variance)), by = scale]
# setnames(cdsAncWavCor, "V1", "cds_anc_cor")
# #plot(y = cdsAncWavCor$cds_anc_cor, x=1:14)
# 
# cdsRecWavCor <- allVarsWavVarCov[, cds_rec_covariance/(sqrt(cds_variance)*sqrt(rec_variance)), by = scale]
# setnames(cdsRecWavCor, "V1", "cds_rec_cor")
# #plot(y = cdsRecWavCor$cds_rec_cor, x=1:14)
# 
# cdsPerCmAncWavCor <- allVarsWavVarCov[, cdsPerCm_anc_covariance/(sqrt(cdsPerCm_variance)*sqrt(anc_variance)), by = scale]
# setnames(cdsPerCmAncWavCor, "V1", "cdsPerCm_anc_cor")
# #plot(y = cdsPerCmAncWavCorFinal$cdsPerCm_anc_cor, x=1:14)
# 
# # combine wavelet correlation tables
# allVarsWavCor <- Reduce(merge, list(recAncWavCor, cdsAncWavCor, cdsRecWavCor, cdsPerCmAncWavCor))
# 
# 
# # 8.4. ---------- Chromosome-Scale correlations ----------
# 
# 
# recAncChrCor <- chrSignalMeans[, sum( (meanFreq-weighted.mean(meanFreq, w=weight))*(cmTr-weighted.mean(cmTr, w=weight)) )/
#                  sqrt( sum((meanFreq - weighted.mean(meanFreq,w=weight))^2) * sum((cmTr - weighted.mean(cmTr,w=weight))^2))]
# 
# cdsAncChrCor <- chrSignalMeans[, sum( (meanFreq-weighted.mean(meanFreq, w=weight))*(coding_bp-weighted.mean(coding_bp, w=weight)) )/
#                                  sqrt( sum((meanFreq - weighted.mean(meanFreq,w=weight))^2) * sum((coding_bp - weighted.mean(coding_bp,w=weight))^2))]
# 
# cdsRecChrCor <- chrSignalMeans[, sum( (cmTr-weighted.mean(cmTr, w=weight))*(coding_bp-weighted.mean(coding_bp, w=weight)) )/
#                                  sqrt( sum((cmTr - weighted.mean(cmTr,w=weight))^2) * sum((coding_bp - weighted.mean(coding_bp,w=weight))^2))]
# 
# cdsPerCmChrCor <- chrSignalMeans[, sum( (meanFreq-weighted.mean(meanFreq, w=weight))*((coding_bp/cmTr)-weighted.mean((coding_bp/cmTr), w=weight)) )/
#                                  sqrt( sum((meanFreq - weighted.mean(meanFreq,w=weight))^2) * sum(((coding_bp/cmTr) - weighted.mean((coding_bp/cmTr),w=weight))^2))]
# 
# allVarsCorDecomp <- rbind(allVarsWavCor,
#                           data.table(scale = "chr",
#                                      rec_anc_cor = recAncChrCor,
#                                      cds_anc_cor = cdsAncChrCor,
#                                      cds_rec_cor = cdsRecChrCor,
#                                      cdsPerCm_anc_cor = cdsPerCmChrCor))
# 
# # ========== Save Correlation Analysis Outputs ==========
# save(allVarsCorDecomp, file = paste0("ACUA_",year,"/anc_rec_cds_corDecomp.RData"))
# 



