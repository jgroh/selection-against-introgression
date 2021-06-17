# 1. ========== Load Dependencies and Data ==========
library(data.table)
library(tools)
library(waveslim)
library(stringi)

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]

scaffFiles <- dir(path=paste0("ACUA_",year),pattern="ScyDAA6*",full.names=T)
scaffs <- basename(file_path_sans_ext(scaffFiles))

# read cds density file
cdsCM <- fread("xbir_1kb_CDS_and_cM.txt")

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
cdsCM[, position := start + 500] # this is the midpoint of the 1kb intervals where I interpolate ancestry
gnomP <- merge(gnomP, cdsCM, by = c("chr", "position"))


# 2. ===== Data Transformation ===== 
# logit transform ancestry proportion
epsilon <- gnomP[indivFreq > 0, min(indivFreq)]/2
gnomP[, indivFreqTr := log(indivFreq/(1-indivFreq))]
gnomG[, indivFreqTr := log(indivFreq/(1-indivFreq))]

gnomP[indivFreq == 0, indivFreqTr := log(epsilon/(1-epsilon))]
gnomG[indivFreq == 0, indivFreqTr := log(epsilon/(1-epsilon))]

gnomP[indivFreq == 1, indivFreqTr := log((1-epsilon)/epsilon)]
gnomG[indivFreq == 1, indivFreqTr := log((1-epsilon)/epsilon)]

# mean transformed ancestry
gnomP[,meanFreqTr := mean(indivFreqTr), by = .(position,chr)]

gnomG[, position := seq_len(.N), by = .(ID,chr)]
gnomG[, meanFreqTr := mean(indivFreqTr), by = .(position,chr)]

# log transform recombination rates
gnomP[, cmTr := log(cM)]

#ggplot(gnomP[sample(1:8e6,1000)], aes(x = log(coding_bp+1), y = meanFreqTr, colour = chr)) + 
#  geom_point() + geom_smooth(method = "lm")

#ggplot(gnomP[chr == chr[1] & position < 1000000], aes(x= position, y = coding_bp)) + geom_point()

# 3. ========== Wavelet Functions ==========

# To compute modwt on variable in a data table:
modwtAllScales <- function(x,variable,lenCol,allcols){
  # can work for wavelet decomp on genetic and physical scale
  dt <- setDT( 
    x[,brick.wall(wf="haar", x=modwt(variable-mean(variable),"haar",n.levels=floor(log2(length(unique(lenCol))))))] 
    )
  # remove smooth coeff column 
  smooth.col <- grep("s",names(dt))
  dt <- dt[,-smooth.col,with=F] 
  
  # add empty detail columns for higher scales not present on particular chromosome
  naCols <- setdiff(allcols, names(dt))
  if(length(naCols) > 0){ dt[,(naCols) := as.numeric(NA)] }
  return(dt)
}

dwt.nondyadic <- function(x){
  # discrete wavelet transform for signal that is not length power of 2
  # (allows for maximum number of coefficients at each scale w/o extending past end of chrom)
  x <- x-mean(x)
  M <- length(x)
  N <- ceiling(log(M, 2)) # next highest power of 2
  J <- floor(log2(length(x))) # max power of 2
  xx <- c(x, rep(0, 2^N - M)) # append zeros to end of sequence to reach length N
  y <- dwt(xx, wf = "haar", n.levels = N) # a list with coefficients for each scale
  
  for(j in 1:(length(y))){
    if(j <= J){
      y[[j]] <- y[[j]][1:floor(M/2^j)] # truncates to only those coefficients not affected by appending zeros to end
    } 
  }
  y[(J+1):length(y)] <- NULL # also changed to simply remove other coefficients

  return(y)
}

# To get number of non-boundary coefficients
numCoeff <- function(x){length(x[!is.na(x)])} 

# Unbiased estimator of wavelet variance for MODWT with brick wall boundary condition
wav_var <- function(x){sum(x^2,na.rm=TRUE)/(length(x[!is.na(x)]))} 

# To compute dwt on variable in a data table, for a single group
dwtAllScales <- function(x,variable,allcols){
  y <- x[, dwt.nondyadic(variable)]
  dt <- as.data.table(stri_list2matrix(y))
  setnames(dt, names(y))
  
  # add empty detail columns for higher scales not present on particular chromosome
  naCols <- setdiff(allcols, names(dt))
  if(length(naCols) > 0){ dt[,(naCols) := as.numeric(NA)] }
  return(dt[,lapply(.SD, as.double)])
}

# 4. ============ DWT: Ancestry, recombination, and coding bp ==========
# define levels based on longest chromosome
maxLevelP <- max(gnomP[ID==ID[1],floor(log2(length(position))),by=chr][,2])
allColsP <- paste0("d",1:maxLevelP)

dwtAnc <- gnomP[ID==ID[1], dwtAllScales(.SD, variable = meanFreqTr, allcols = allColsP), by = chr]
dwtRec <- gnomP[ID==ID[1], dwtAllScales(.SD, variable = cmTr, allcols = allColsP), by = chr]
dwtCds <- gnomP[ID==ID[1], dwtAllScales(.SD, variable = coding_bp,  allcols = allColsP), by = chr]

# combine tables
lapply(list(dwtAnc, dwtRec, dwtCds), function(x){ x[,position := seq_len(.N), by = chr] })

dwtAnc <- melt(dwtAnc, id.vars = c("chr","position"), 
               variable.name = "scale",
               value.name = "ancestry_coeff")

dwtRec <- melt(dwtRec, id.vars = c("chr","position"), 
               variable.name = "scale",
               value.name = "rec_coeff")
dwtCds <- melt(dwtCds, id.vars = c("chr","position"), 
               variable.name = "scale",
               value.name = "cds_coeff")

allDWT <- Reduce(function(...) merge(..., all=TRUE), list(dwtAnc, dwtRec, dwtCds))
save(allDWT, file = paste0("ACUA_",year,"/anc_rec_cds_DWT.RData"))


# 5. ========== MODWT: Genetic Scale ==========
# define levels based on longest chromosome
maxLevelG <- max(gnomG[,floor(log2(length(unique(Morgan)))),by=chr][,2])
allColsG <- paste0("d",1:maxLevelG)

# 5.1. ---------- MODWT on individuals ----------
indModwtG <- gnomG[, modwtAllScales(.SD,variable=indivFreqTr,lenCol=Morgan,allcols=allColsG), by = .(chr,ID)]

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


# 5.2. ---------- MODWT on population mean ----------

# run modwt
popMeanModwtG <- gnomG[ID==ID[1],modwtAllScales(.SD,variable=meanFreqTr,lenCol=Morgan,allcols=allColsG), by = chr]

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


# 6. ========== Chromosome-Level analysis: Genetic Scale ==========

# 6.1. ---------- Population Mean ----------

# total genetic length as weights for chromosomes 
chrLenG <- gnomG[, max(Morgan), by = .(chr)]
setnames(chrLenG, "V1", "len")
chrLenG[, weight := len/sum(len)]

# chromosome means of population mean minor parent ancestry
chrMeansPopMeanG <- gnomG[, mean(meanFreqTr), by = .(chr)]
setnames(chrMeansPopMeanG, "V1", "avg_frq")
chrMeansPopMeanG <- merge(chrMeansPopMeanG, chrLenG)

# weighted average of population-mean ancestry across chromosomes 
weightedMeanPopMeanG <- chrMeansPopMeanG[, weighted.mean(avg_frq, weight)]

# chromosome-level weighted variance of pop mean ancestry
chrVarPopG <- chrMeansPopMeanG[, sum(weight*(avg_frq - weightedMeanPopMeanG)^2)]


# 6.2. ---------- Individual-Level ----------

# chromosome-level average minor parent ancestry for individuals
chrMeansIndG <- gnomG[, mean(indivFreqTr), by = .(ID,chr)]
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


# 7. ========== MODWT of ancestry on Physical Scale ==========

# 7.1. ---------- MODWT on Individuals ----------
indModwtP <- gnomP[, modwtAllScales(.SD,variable=indivFreqTr,lenCol=position,allcols=allColsP), by = .(chr,ID)]

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
indMeanWavVarP[, distance := "physical"]

# 7.2. --------- MODWT on Mean Ancestry ----------
popModwtP <- gnomP[ID==ID[1], modwtAllScales(.SD,variable=meanFreqTr,lenCol=position,allcols=allColsP), by = chr]
popWavVarP <- popModwtP[,lapply(.SD,wav_var),by=chr]
setnames(popWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
popWavVarP <- melt(popWavVarP, measure.vars = as.character(1:maxLevelP),
                         variable.name = "scale", value.name = "anc_variance")

# 5.3. ----- MODWT on Recomb Rate (cM length of 1kb windows) ----------
recModwtP <- gnomP[ID==ID[1],modwtAllScales(.SD,variable=cmTr,lenCol=position,allcols=allColsP), by = chr]

# wavelet variance for rec rate
recWavVarP <- recModwtP[,lapply(.SD,wav_var),by=chr]
setnames(recWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
recWavVarP <- melt(recWavVarP, measure.vars = as.character(1:maxLevelP),
                   variable.name = "scale", value.name = "rec_variance")

# 5.4. ----- MODWT on # of coding bp per 1kb window ----------
cdsModwtP <- gnomP[ID==ID[1],modwtAllScales(.SD,variable=coding_bp,lenCol=position,allcols=allColsP), by = chr]

# wavelet variance for CDS
cdsWavVarP <- cdsModwtP[,lapply(.SD,wav_var),by=chr]
setnames(cdsWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
cdsWavVarP <- melt(cdsWavVarP, measure.vars = as.character(1:maxLevelP),
                   variable.name = "scale", value.name = "cds_variance")

# 5.5. ----- MODWT on CDS per cM coding bp --------------
cdsPerCmModwtP <- gnomP[ID==ID[1],modwtAllScales(.SD,variable=coding_bp/cmTr,lenCol=position,allcols=allColsP), by = chr]

# wavelet variance for cds density
cdsPerCmWavVarP <- cdsPerCmModwtP[,lapply(.SD,wav_var),by=chr]
setnames(cdsPerCmWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
cdsPerCmWavVarP <- melt(cdsPerCmWavVarP, measure.vars = as.character(1:maxLevelP),
                        variable.name = "scale", value.name = "cdsPerCm_variance")

# 5.6. ------- Combine tables with wavelet variances for all signals -------
# save object with chromosome-specific variances of signals
allVarsWavVar_Chrs <- merge(Reduce(merge, list(popWavVarP, recWavVarP, cdsWavVarP, cdsPerCmWavVarP)),
                            chrWeightsP)
allVarsWavVar_Chrs[,`:=`(distance = "physical", 
                         decomp = "pop_mean")]


# weighted average over chromosomes
signalVars <- c("anc_variance", "rec_variance", "cds_variance", "cdsPerCm_variance")
allVarsWavVar <- allVarsWavVar_Chrs[, lapply(.SD, weighted.mean, w=weight), .SDcols = signalVars, by = .(scale,distance,decomp)]


# 6. ========== Chromosome-Level Analysis: Physical Scale ==========

# total physical length as weights for chromosomes 
chrLenP <- gnomP[, max(position), by = .(chr)]
setnames(chrLenP, "V1", "len")
chrLenP[, weight := len/sum(len)]

# 6.1. ---------- Chromosome-level variance of mean ancestry, rec, cds ----------
signalCols <- c("meanFreqTr", "cmTr", "coding_bp")

# computed chromosome means
chrSignalMeans <- gnomP[, lapply(.SD, mean), .SDcols = signalCols, by = .(chr)]
chrSignalMeans[, weight := chrLenP[,weight]]

# compute weighted variance 
chrVarP <- chrSignalMeans[, lapply(.SD, function(x){sum(weight*(x - weighted.mean(x, weight))^2)}),
                                          .SDcols = signalCols]
setnames(chrVarP, c("anc_variance", "rec_variance", "cds_variance"))
chrVarP[, `:=`(decomp = "pop_mean", 
               distance = "physical",
               scale = "chr")]

# 6.2. ---------- Ancestry: Individual-Level ----------

# average minor parent ancestry for individuals for each chrom
chrMeansIndP <- gnomP[, mean(indivFreqTr), by = .(ID,chr)]
setnames(chrMeansIndP, "V1", "avg_frq")

# calculate genome-wide mean for individuals by weighting chromosomes
chrMeansIndP <- merge(chrMeansIndP, chrLenP, by = "chr") # add chrom weights
chrMeansIndP[, gnomWideMean := weighted.mean(avg_frq,weight), by = ID]

# good place to look at distribution of genome-wide admixture proportions
#boxplot(chrMeansIndP$gnomWideMean)

# chromosome-level weighted variance by individual, then average over individuals
chrVarIndP <- chrMeansIndP[, sum(weight*(avg_frq - gnomWideMean)^2), by = ID][,mean(V1)]


# 6.3 ------ Combine Chr-level Variances for individual and population ----------
chrVarP <- rbind(list(anc_variance = chrVarIndP, 
                      rec_variance = NA, 
                      cds_variance = NA, 
                      decomp = "mean_individual", 
                      distance = "physical", 
                      scale = "chr"), chrVarP)


# 7. ========== Output Wavelet and Chrom-Level Variance Data ==========

# Combine Wavelet variance tables
wvFinalG <- rbind(indMeanWavVarG, popMeanWavVarG)
wvFinalG[, distance := "genetic"] 
wvFinalP <- merge(indMeanWavVarP, allVarsWavVar, by = c("scale", "decomp", "distance", "anc_variance"), all=T)
wvFinalAll <- merge(wvFinalP, wvFinalG, by = c("scale", "decomp", "distance", "anc_variance"), all =T)

# combine wavelet variances by chromosome (using mean ancestry)
wvChrsAll <- merge(popWavVarG_Chrs, allVarsWavVar_Chrs, by = c("chr", "scale", "distance", "anc_variance"), all = T)

# combine chromosome-level variances
chrVarAll <- merge(chrVarP, chrVarG, by = c("scale","distance", "anc_variance", "decomp"),all=T)

save(wvFinalAll, chrVarAll, wvChrsAll, file = paste0("ACUA_",year,"/anc_rec_cds_varDecomp.RData"))


# 8. ========== Correlation Analysis ==========

# 8.1. ---------- Reformat MODWT tables ----------
popModwtP[, position := seq_len(.N), by = chr]
setnames(popModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
popModwtP <- melt(popModwtP, measure.vars = as.character(1:maxLevelP),
                  variable.name = "scale", value.name = "ancestry_coeff")

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
# combine
allVarsModwtP <- Reduce(merge, list(popModwtP, recModwtP, cdsModwtP, cdsPerCmModwtP))
 
# save modwt coefficients for linear model analysis
save(allVarsModwtP, file = paste0("ACUA_",year,"/anc_rec_cds_MODWT.RData"))

# 8.2. --------- Wavelet covariances -----------

# unweighted wavelet covariances at each scale separately for each chromosome
recAncCov <- allVarsModwtP[, sum(rec_coeff*ancestry_coeff, na.rm=T)/length(ancestry_coeff[!is.na(ancestry_coeff)]), by = .(chr, scale)]
setnames(recAncCov, "V1", "rec_anc_covariance")

cdsAncCov <- allVarsModwtP[, sum(cds_coeff*ancestry_coeff, na.rm=T)/length(ancestry_coeff[!is.na(ancestry_coeff)]), by = .(chr, scale)]
setnames(cdsAncCov, "V1", "cds_anc_covariance")

cdsRecCov <- allVarsModwtP[, sum(cds_coeff*rec_coeff, na.rm=T)/length(ancestry_coeff[!is.na(ancestry_coeff)]), by = .(chr, scale)]
setnames(cdsRecCov, "V1", "cds_rec_covariance")

cdsPerCmAncCov <- allVarsModwtP[, sum(cdsPerCm_coeff*ancestry_coeff, na.rm=T)/length(cds_coeff[!is.na(cdsPerCm_coeff)]), by = .(chr, scale)]
setnames(cdsPerCmAncCov, "V1", "cdsPerCm_anc_covariance")

# combine covariance values
allVarsCov <- Reduce(merge, list(recAncCov, cdsAncCov, cdsRecCov, cdsPerCmAncCov))

# Weighted average of covariances over chromosomes
allVarsCov <- merge(allVarsCov, chrWeightsP)[, lapply(.SD, weighted.mean, w=weight), 
                               .SDcols = c("rec_anc_covariance","cds_anc_covariance",
                                           "cds_rec_covariance","cdsPerCm_anc_covariance"),
                               by = scale]

# merge covariance table with wavelet variance table
allVarsWavVarCov <- merge(allVarsWavVar, allVarsCov)


# 8.3. ----------- Compute Wavelet Correlations between signal pairs ------------
recAncWavCor <- allVarsWavVarCov[, rec_anc_covariance/(sqrt(rec_variance)*sqrt(anc_variance)), by = scale]
setnames(recAncWavCor, "V1", "rec_anc_cor")
#plot(y = recAncWavCor$rec_anc_cor, x=1:14)

cdsAncWavCor <- allVarsWavVarCov[, cds_anc_covariance/(sqrt(cds_variance)*sqrt(anc_variance)), by = scale]
setnames(cdsAncWavCor, "V1", "cds_anc_cor")
#plot(y = cdsAncWavCor$cds_anc_cor, x=1:14)

cdsRecWavCor <- allVarsWavVarCov[, cds_rec_covariance/(sqrt(cds_variance)*sqrt(rec_variance)), by = scale]
setnames(cdsRecWavCor, "V1", "cds_rec_cor")
#plot(y = cdsRecWavCor$cds_rec_cor, x=1:14)

cdsPerCmAncWavCor <- allVarsWavVarCov[, cdsPerCm_anc_covariance/(sqrt(cdsPerCm_variance)*sqrt(anc_variance)), by = scale]
setnames(cdsPerCmAncWavCor, "V1", "cdsPerCm_anc_cor")
#plot(y = cdsPerCmAncWavCorFinal$cdsPerCm_anc_cor, x=1:14)

# combine wavelet correlation tables
allVarsWavCor <- Reduce(merge, list(recAncWavCor, cdsAncWavCor, cdsRecWavCor, cdsPerCmAncWavCor))


# 8.4. ---------- Chromosome-Scale correlations ----------


recAncChrCor <- chrSignalMeans[, sum( (meanFreqTr-weighted.mean(meanFreqTr, w=weight))*(cmTr-weighted.mean(cmTr, w=weight)) )/
                 sqrt( sum((meanFreqTr - weighted.mean(meanFreqTr,w=weight))^2) * sum((cmTr - weighted.mean(cmTr,w=weight))^2))]

cdsAncChrCor <- chrSignalMeans[, sum( (meanFreqTr-weighted.mean(meanFreqTr, w=weight))*(coding_bp-weighted.mean(coding_bp, w=weight)) )/
                                 sqrt( sum((meanFreqTr - weighted.mean(meanFreqTr,w=weight))^2) * sum((coding_bp - weighted.mean(coding_bp,w=weight))^2))]

cdsRecChrCor <- chrSignalMeans[, sum( (cmTr-weighted.mean(cmTr, w=weight))*(coding_bp-weighted.mean(coding_bp, w=weight)) )/
                                 sqrt( sum((cmTr - weighted.mean(cmTr,w=weight))^2) * sum((coding_bp - weighted.mean(coding_bp,w=weight))^2))]

cdsPerCmChrCor <- chrSignalMeans[, sum( (meanFreqTr-weighted.mean(meanFreqTr, w=weight))*((coding_bp/cmTr)-weighted.mean((coding_bp/cmTr), w=weight)) )/
                                 sqrt( sum((meanFreqTr - weighted.mean(meanFreqTr,w=weight))^2) * sum(((coding_bp/cmTr) - weighted.mean((coding_bp/cmTr),w=weight))^2))]

allVarsCorDecomp <- rbind(allVarsWavCor,
                          data.table(scale = "chr",
                                     rec_anc_cor = recAncChrCor,
                                     cds_anc_cor = cdsAncChrCor,
                                     cds_rec_cor = cdsRecChrCor,
                                     cdsPerCm_anc_cor = cdsPerCmChrCor))

# ========== Save Correlation Analysis Outputs ==========
save(allVarsCorDecomp, file = paste0("ACUA_",year,"/anc_rec_cds_corDecomp.RData"))




