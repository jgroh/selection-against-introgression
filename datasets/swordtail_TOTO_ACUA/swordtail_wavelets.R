# 1. ========== Load Dependencies and Data ==========
library(data.table)
library(tools)
library(waveslim)

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]

scaffFiles <- dir(path=paste0("ACUA_",year),pattern="ScyDAA6*",full.names=T)
scaffs <- basename(file_path_sans_ext(scaffFiles))

# read cds density file
cdsCM <- fread("allChr_cds_and_cm_per_1kb.txt")

# each file has the same object name for the scaffold 
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
cdsCM[, position := start + 500]
gnomP <- merge(gnomP, cdsCM, by = c("chr", "position"))


# 2. ========== MODWT Functions ==========

# To compute modwt on variable in a data table:
modwtAllScales <- function(x,variable,lenCol,allcols){
  # can work for wavelet decomp on genetic and physical scale
  dt <- setDT( 
    x[,brick.wall(wf="haar", x=modwt(variable,"haar",n.levels=floor(log2(length(unique(lenCol))))))] 
    )
  # remove smooth coeff column 
  smooth.col <- grep("s",names(dt))
  dt <- dt[,-smooth.col,with=F] 
  
  # add empty detail columns for higher scales not present on particular chromosome
  naCols <- setdiff(allcols, names(dt))
  if(length(naCols) > 0){ dt[,(naCols) := as.numeric(NA)] }
  return(dt)
}

# To get number of non-boundary coefficients
numCoeff <- function(x){length(x[!is.na(x)])} 

# Unbiased estimator of wavelet variance for MODWT with brick wall boundary condition
wav_var <- function(x){sum(x^2,na.rm=TRUE)/(length(x[!is.na(x)]))} 


# 3. ========== MODWT: Genetic Scale ==========
# define levels based on longest chromosome
maxLevelG <- max(gnomG[,floor(log2(length(unique(Morgan)))),by=chr][,2])
allColsG <- paste0("d",1:maxLevelG)

# 3.1. ---------- MODWT on individuals ----------
indModwtG <- gnomG[, modwtAllScales(.SD,variable=indivFreq,lenCol=Morgan,allcols=allColsG), by = .(chr,ID)]

# number of wavelet coefficients per scale on each chromosome
chrWeightsG <- indModwtG[ID==ID[1],lapply(.SD,numCoeff), by = .(chr,ID)]
setnames(chrWeightsG, old = paste0("d",1:maxLevelG), new = as.character(1:maxLevelG))
chrWeightsG <- melt(chrWeightsG[,-"ID"], id.vars = "chr", measure.vars = as.character(1:maxLevelG),
                         variable.name = "scale", value.name = "numCoeffs")

# obtain final weights for each chromosome by scale
chrWeightsG[,weight := numCoeffs/sum(numCoeffs), by=.(scale)]

# wavelet variance for individual x chromosome
indWavVarG <- indModwtG[,lapply(.SD,wav_var), by = .(chr,ID)]
setnames(indWavVarG, old = paste0("d",1:maxLevelG), new = as.character(1:maxLevelG))
indWavVarG <- melt(indWavVarG, id.vars = c("ID","chr"), 
                    variable.name = "scale",
                    value.name = "variance")

# average over chromosomes within an individual, then average over individuals
indWavVarG <- merge(indWavVarG, chrWeightsG) 
indMeanWavVarG <- indWavVarG[, weighted.mean(variance, weight), by = .(ID,scale)][,mean(V1), by = scale]
setnames(indMeanWavVarG, "V1", "anc_variance")
indMeanWavVarG[,decomp :="mean_individual"]


# 3.2. ---------- MODWT on population mean ----------

# run modwt
popMeanModwtG <- gnomG[ID==ID[1],modwtAllScales(.SD,variable=meanFreq,lenCol=Morgan,allcols=allColsG), by = chr]

# wavelet variance by chromosome
popMeanWavVarG <- popMeanModwtG[,lapply(.SD,wav_var),by=chr]
setnames(popMeanWavVarG, paste0("d",1:maxLevelG), as.character(1:maxLevelG))
popMeanWavVarG <- melt(popMeanWavVarG, measure.vars = as.character(1:maxLevelG),
                         variable.name = "scale", value.name = "anc_variance")
# save the chromosome-level version before averaging for looking at chrs separately
popWavVarG_Chrs <- popMeanWavVarG
popWavVarG_Chrs[, distance := "genetic"]

# weighted average over chromosomes
popMeanWavVarG <- merge(popMeanWavVarG, chrWeightsG) 
popMeanWavVarG <- popMeanWavVarG[, weighted.mean(anc_variance, weight), by = scale]
setnames(popMeanWavVarG, "V1", "anc_variance")
popMeanWavVarG[,decomp :="pop_mean"]


# 4. ========== Chromosome-Level analysis: Genetic Scale ==========

# 4.1. ---------- Population Mean ----------

# total genetic length as weights for chromosomes 
chrLenG <- gnomG[, max(Morgan), by = .(chr)]
setnames(chrLenG, "V1", "len")
chrLenG[, weight := len/sum(len)]

# chromosome means of population mean minor parent ancestry
chrMeansPopMeanG <- gnomG[, mean(meanFreq), by = .(chr)]
setnames(chrMeansPopMeanG, "V1", "avg_frq")
chrMeansPopMeanG <- merge(chrMeansPopMeanG, chrLenG)

# weighted average of population-mean ancestry across chromosomes 
weightedMeanPopMeanG <- chrMeansPopMeanG[, weighted.mean(avg_frq, weight)]

# chromosome-level weighted variance of pop mean ancestry
chrVarPopG <- chrMeansPopMeanG[, sum(weight*(avg_frq - weightedMeanPopMeanG)^2)]


# 4.2. ---------- Individual-Level ----------

# chromosome-level average minor parent ancestry for individuals
chrMeansIndG <- gnomG[, mean(indivFreq), by = .(ID,chr)]
setnames(chrMeansIndG, "V1", "avg_frq")

# calculate genome-wide mean for individuals by weighting chromosomes
chrMeansIndG <- merge(chrMeansIndG, chrLenG) # add chrom weights
chrMeansIndG[, gnomWideMean := weighted.mean(avg_frq,weight), by = ID]

# chromosome-level weighted variance by individual, then average over individuals
chrVarIndG <- chrMeansIndG[, sum(weight*(avg_frq - gnomWideMean)^2), by = ID][,mean(V1)]

# output weighted chromosome-level variance for avg. individual-level signal 
# and population mean 
chrVarG <- data.table("anc_variance" = c(chrVarIndG, chrVarPopG), 
           "decomp" = c("mean_individual", "pop_mean"), 
           "scale" = "chr", "distance" = "genetic")


# 5. ========== MODWT of ancestry on Physical Scale ==========

# define levels based on longest chromosome
maxLevelP <- max(gnomP[,floor(log2(length(unique(position)))),by=chr][,2])
allColsP <- paste0("d",1:maxLevelP)

# 5.1. ---------- MODWT on Individuals ----------
indModwtP <- gnomP[, modwtAllScales(.SD,variable=indivFreq,lenCol=position,allcols=allColsP), by = .(chr,ID)]

# obtain chromosome weights for each scale based on number of non-boundary coefficients on the chromosome
chrWeightsP <- indModwtP[ID==ID[1],lapply(.SD,numCoeff), by = .(chr,ID)]
setnames(chrWeightsP, old = paste0("d",1:maxLevelP), new = as.character(1:maxLevelP))
chrWeightsP <- melt(chrWeightsP[,-"ID"], id.vars = "chr", measure.vars = as.character(1:maxLevelP),
                   variable.name = "scale", value.name = "numCoeffs")
chrWeightsP[,weight := numCoeffs/sum(numCoeffs), by=.(scale)]

# compute wavelet variance for individual x chromosome
indWavVarP <- indModwtP[,lapply(.SD,wav_var), by = .(chr,ID)]
setnames(indWavVarP, old = paste0("d",1:maxLevelP), new = as.character(1:maxLevelP))
indWavVarP <- melt(indWavVarP, id.vars = c("ID","chr"), 
                    variable.name = "scale",
                    value.name = "variance")

# average over chromosomes within an individual, then average over individuals
indWavVarP <- merge(indWavVarP, chrWeightsP) 
indMeanWavVarP <- indWavVarP[, weighted.mean(variance, weight), by = .(ID,scale)][,mean(V1), by = scale]
setnames(indMeanWavVarP, "V1", "anc_variance")
indMeanWavVarP[, decomp := "mean_individual"]


# 5.2. --------- MODWT on Population Mean ----------

# run wavelet decomp on population mean for each chrom separately
popModwtP <- gnomP[ID==ID[1], modwtAllScales(.SD,variable=meanFreq,lenCol=position,allcols=allColsP), by = chr]
popWavVarP <- popModwtP[,lapply(.SD,wav_var),by=chr]
setnames(popWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
popWavVarP <- melt(popWavVarP, measure.vars = as.character(1:maxLevelP),
                         variable.name = "scale", value.name = "anc_variance")
# save the chromosome-level version before averaging for looking at chrs separately
popWavVarP_Chrs <- popWavVarP
popWavVarP_Chrs[, distance := "physical"]

# weighted average over chromosomes
popWavVarP <- merge(popWavVarP, chrWeightsP) 
popMeanWavVarP <- popWavVarP[, weighted.mean(anc_variance, weight), by = scale]
setnames(popMeanWavVarP, "V1", "anc_variance")
popMeanWavVarP[, decomp := "pop_mean"]



# 6. ========== Chromosome-Level Analysis: Physical Scale ==========

# 6.1. ---------- Population Mean ----------

# total physical length as weights for chromosomes 
chrLenP <- gnomP[, max(position), by = .(chr)]
setnames(chrLenP, "V1", "len")
chrLenP[, weight := len/sum(len)]

# chromosome means of population mean minor parent ancestry
chrMeansPopMeanP <- gnomP[, mean(meanFreq), by = .(chr)]
setnames(chrMeansPopMeanP, "V1", "avg_frq")
chrMeansPopMeanP[, weight := chrLenP[,weight]]

# weighted average of population-mean ancestry across chromosomes 
weightedMeanPopMeanP <- chrMeansPopMeanP[, weighted.mean(avg_frq, weight)]

# chromosome-level weighted variance of pop mean ancestry
chrVarPopP <- chrMeansPopMeanP[, sum(weight*(avg_frq - weightedMeanPopMeanP)^2)]


# 6.2. ---------- Individual-Level ----------

# chromosome-level average minor parent ancestry for individuals
chrMeansIndP <- gnomP[, mean(indivFreq), by = .(ID,chr)]
setnames(chrMeansIndP, "V1", "avg_frq")

# calculate genome-wide mean for individuals by weighting chromosomes
chrMeansIndP <- merge(chrMeansIndP, chrLenP, by = "chr") # add chrom weights
chrMeansIndP[, gnomWideMean := weighted.mean(avg_frq,weight), by = ID]

# chromosome-level weighted variance by individual, then average over individuals
chrVarIndP <- chrMeansIndP[, sum(weight*(avg_frq - gnomWideMean)^2), by = ID][,mean(V1)]

# output weighted chromosome-level variance for avg. individual-level signal 
# and population mean 
chrVarP <- data.table("anc_variance" = c(chrVarIndP, chrVarPopP), 
                      "decomp" = c("mean_individual", "pop_mean"), 
                      "scale" = "chr", "distance" = "physical")

# 7. ========== Output Wavelet and Chrom-Level Variance Data ==========

# Combine Wavelet variance tables
wvFinalG <- rbind(indMeanWavVarG, popMeanWavVarG)
wvFinalG[, distance := "genetic"] 
wvFinalP <- rbind(indMeanWavVarP, popMeanWavVarP)
wvFinalP[, distance := "physical"]
wvFinalAll <- rbind(wvFinalP, wvFinalG)

# combine wavelet variances by chromosome
wvChrsAll <- rbind(popWavVarG_Chrs, popWavVarP_Chrs)

# combine chromosome-level variances
chrVarAll <- rbind(chrVarP, chrVarG)

save(wvFinalAll, chrVarAll, wvChrsAll, file = paste0("ACUA_",year,"/ancVarDecomp.RData"))


# 8. ========== Correlation Analysis ==========

# 8.1. ---------- Wavelet Correlation Analysis ----------

# run modwt on recombination rate (cM length of 1kb window)
recModwtP <- gnomP[ID==ID[1],modwtAllScales(.SD,variable=cM,lenCol=position,allcols=allColsP), by = chr]

# wavelet variance for rec rate
recWavVarP <- recModwtP[,lapply(.SD,wav_var),by=chr]
setnames(recWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
recWavVarP <- melt(recWavVarP, measure.vars = as.character(1:maxLevelP),
                       variable.name = "scale", value.name = "rec_variance")

# run modwt on coding bp per 1kb window
cdsModwtP <- gnomP[ID==ID[1],modwtAllScales(.SD,variable=coding_bp,lenCol=position,allcols=allColsP), by = chr]

# wavelet variance for cds density
cdsWavVarP <- cdsModwtP[,lapply(.SD,wav_var),by=chr]
setnames(cdsWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
cdsWavVarP <- melt(cdsWavVarP, measure.vars = as.character(1:maxLevelP),
                 variable.name = "scale", value.name = "cds_variance")

# run modwt on cM per kb coding bp
cdsPerCmModwtP <- gnomP[ID==ID[1],modwtAllScales(.SD,variable=coding_bp/cM,lenCol=position,allcols=allColsP), by = chr]

# wavelet variance for cds density
cdsPerCmWavVarP <- cdsPerCmModwtP[,lapply(.SD,wav_var),by=chr]
setnames(cdsPerCmWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
cdsPerCmWavVarP <- melt(cdsPerCmWavVarP, measure.vars = as.character(1:maxLevelP),
                   variable.name = "scale", value.name = "cdsPerCm_variance")

# combine tables with wavelet variances for all signals
allVarsWavVar <- Reduce(merge, list(popWavVarP, recWavVarP, cdsWavVarP, cdsPerCmWavVarP))

# get modwt data tables into more convenient shape (only after calculating wav var)
popModwtP[, position := seq_len(.N), by = chr]
setnames(popModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
popModwtP <- melt(popModwtP, measure.vars = as.character(1:maxLevelP),
                  variable.name = "scale", value.name = "anc_coeff")

recModwtP[, position := seq_len(.N), by = chr]
setnames(recModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
recModwtP <- melt(recModwtP, measure.vars = as.character(1:maxLevelP),
                   variable.name = "scale", value.name = "rec_coeff")

cdsModwtP[, position := seq_len(.N), by = chr]
setnames(cdsModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
cdsModwtP <- melt(cdsModwtP, measure.vars = as.character(1:maxLevelP),
                  variable.name = "scale", value.name = "cds_coeff")

cdsPerCmModwtP[, position := seq_len(.N), by = chr]
setnames(cdsPerCmModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
cdsPerCmModwtP <- melt(cdsPerCmModwtP, measure.vars = as.character(1:maxLevelP),
                  variable.name = "scale", value.name = "cdsPerCm_coeff")


allVarsModwtP <- Reduce(merge, list(popModwtP, recModwtP, cdsModwtP, cdsPerCmModwtP))
 

# unweighted wavelet covariances at each scale separately for each chromosome
recAncCov <- allVarsModwtP[, sum(rec_coeff*anc_coeff, na.rm=T)/length(rec_coeff[!is.na(rec_coeff)]), by = .(chr, scale)]
setnames(recAncCov, "V1", "rec_anc_covariance")

cdsAncCov <- allVarsModwtP[, sum(cds_coeff*anc_coeff, na.rm=T)/length(cds_coeff[!is.na(cds_coeff)]), by = .(chr, scale)]
setnames(cdsAncCov, "V1", "cds_anc_covariance")

cdsRecCov <- allVarsModwtP[, sum(cds_coeff*rec_coeff, na.rm=T)/length(cds_coeff[!is.na(cds_coeff)]), by = .(chr, scale)]
setnames(cdsRecCov, "V1", "cds_rec_covariance")

cdsPerCmAncCov <- allVarsModwtP[, sum(cdsPerCm_coeff*anc_coeff, na.rm=T)/length(cds_coeff[!is.na(cdsPerCm_coeff)]), by = .(chr, scale)]
setnames(cdsPerCmAncCov, "V1", "cdsPerCm_anc_covariance")

# combine covariance values
allVarsCov <- Reduce(merge, list(recAncCov, cdsAncCov, cdsRecCov, cdsPerCmAncCov))

# merge covariance table with wavelet variance table
allVarsWavVarCov <- merge(allVarsWavVar, allVarsCov)

# Weighted mean covariances and variances over chromosomes, weighting by number of wavelets on chromosome at each scale
weightedWavVarCov <- allVarsWavVarCov[, lapply(.SD, weighted.mean, w=weight), by = scale, 
                             .SDcols = c("anc_variance","rec_variance","cds_variance","cdsPerCm_variance",
                                         "rec_anc_covariance","cds_anc_covariance","cds_rec_covariance","cdsPerCm_anc_covariance")]

# Compute final correlation using weights
cdsPerCmAncWavCorFinal <- weightedWavVarCov[, cdsPerCm_anc_covariance/(sqrt(cdsPerCm_variance)*sqrt(anc_variance)), by = scale]
setnames(cdsPerCmAncWavCorFinal, "V1", "cdsPerCm_anc_cor")
#plot(y = cdsPerCmAncWavCorFinal$V1, x=1:14)


# 8.2. ---------- Chromosome-Scale correlation ----------

chrMeansCdsPerCm <- gnomP[, mean(coding_bp/cM), by = chr]
setnames(chrMeansCdsPerCm, "V1", "cdsPerCm")

chrMeansCdsPerCmAnc <- merge(chrMeansPopMeanP, chrMeansCdsPerCm)

# compute weighted means of signals overall all chromosomes
gnomMeansCdsPerCmAnc <- chrMeansCdsPerCmAnc[, lapply(.SD,weighted.mean, weight=weight), .SDcols = c("avg_frq", "cdsPerCm")]

# numerator is a weighted covariance
chrCorNum <- chrMeansCdsPerCmAnc[, sum(weight*(avg_frq-gnomMeansCdsPerCmAnc$avg_frq)*(cdsPerCm-gnomMeansCdsPerCmAnc$cdsPerCm))]

# denom is product of weighted standard deviations of two signals
chrCorDenom <- chrMeansCdsPerCmAnc[, sqrt(
  (sum(weight*(avg_frq - gnomMeansCdsPerCmAnc$avg_frq)^2))*
    (sum(weight*(cdsPerCm - gnomMeansCdsPerCmAnc$cdsPerCm)^2)) )
  ]

cdsPerCmAncChrCor <- chrCorNum/chrCorDenom

# ========== Save Correlation Analysis Outputs ==========
save(cdsPerCmAncWavCorFinal, cdsPerCmAncChrCor, file = paste0("ACUA_",year,"cdsPerCm_x_anc_corDecomp.RData"))
