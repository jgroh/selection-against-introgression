# 1. ========== Load Dependencies and Data ==========
library(data.table)
library(tools)
library(waveslim)

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]

scaffFiles <- dir(path=paste0("ACUA_",year),pattern="ScyDAA6*",full.names=T)
scaffs <- basename(file_path_sans_ext(scaffFiles))

# each file has the same object name for the scaffold 
# so we load into separate environments
loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 

# combine data into one file for the year
gnomG <- rbindlist(lapply(as.list(scaffFiles), 
                  function(x){loadFrom(x, "chrom1.interp.gen")}))
gnomP <- rbindlist(lapply(as.list(scaffFiles), 
                             function(x){loadFrom(x, "chrom1.interp.phys")}))

gnomG[,ind_minor_parent := 1 - ind_frq_interp]
gnomG[,pop_mean_minor_parent := 1 - pop_mean_interp]

gnomP[,ind_minor_parent := 1 - ind_frq_interp]
gnomP[,pop_mean_minor_parent := 1 - pop_mean_interp]

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
indModwtG <- gnomG[, modwtAllScales(.SD,variable=ind_minor_parent,lenCol=Morgan,allcols=allColsG), by = .(chr,ID)]

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
setnames(indMeanWavVarG, "V1", "variance")
indMeanWavVarG[,decomp :="mean_individual"]


# 3.2. ---------- MODWT on population mean ----------

# run modwt
popMeanModwtG <- gnomG[ID==ID[1],modwtAllScales(.SD,variable=pop_mean_minor_parent,lenCol=Morgan,allcols=allColsG), by = chr]

# wavelet variance by chromosome
popMeanWavVarG <- popMeanModwtG[,lapply(.SD,wav_var),by=chr]
setnames(popMeanWavVarG, paste0("d",1:maxLevelG), as.character(1:maxLevelG))
popMeanWavVarG <- melt(popMeanWavVarG, measure.vars = as.character(1:maxLevelG),
                         variable.name = "scale", value.name = "variance")
# save the chromosome-level version before averaging for looking at chrs separately
popWavVarG_Chrs <- popMeanWavVarG

# weighted average over chromosomes
popMeanWavVarG <- merge(popMeanWavVarG, chrWeightsG) 
popMeanWavVarG <- popMeanWavVarG[, weighted.mean(variance, weight), by = scale]
setnames(popMeanWavVarG, "V1", "variance")
popMeanWavVarG[,decomp :="pop_mean"]


# 4. ========== Chromosome-Level analysis: Genetic Scale ==========

# 4.1. ---------- Population Mean ----------

# total genetic length as weights for chromosomes 
chrLenG <- gnomG[, max(Morgan), by = .(chr)]
setnames(chrLenG, "V1", "len")
chrLenG[, weight := len/sum(len)]

# chromosome means of population mean minor parent ancestry
chrMeansPopMeanG <- gnomG[, mean(pop_mean_minor_parent), by = .(chr)]
setnames(chrMeansPopMeanG, "V1", "avg_frq")
chrMeansPopMeanG <- merge(chrMeansPopMeanG, chrLenG)

# weighted average of population-mean ancestry across chromosomes 
weightedMeanPopMeanG <- chrMeansPopMeanG[, weighted.mean(avg_frq, weight)]

# chromosome-level weighted variance of pop mean ancestry
chrVarPopG <- chrMeansPopMeanG[, sum(weight*(avg_frq - weightedMeanPopMeanG)^2)]


# 4.2. ---------- Individual-Level ----------

# chromosome-level average minor parent ancestry for individuals
chrMeansIndG <- gnomG[, mean(ind_minor_parent), by = .(ID,chr)]
setnames(chrMeansIndG, "V1", "avg_frq")

# calculate genome-wide mean for individuals by weighting chromosomes
chrMeansIndG <- merge(chrMeansIndG, chrLenG) # add chrom weights
chrMeansIndG[, gnomWideMean := weighted.mean(avg_frq,weight), by = ID]

# chromosome-level weighted variance by individual, then average over individuals
chrVarIndG <- chrMeansIndG[, sum(weight*(avg_frq - gnomWideMean)^2), by = ID][,mean(V1)]

# output weighted chromosome-level variance for avg. individual-level signal 
# and population mean 
chrVarG <- data.table("variance" = c(chrVarIndG, chrVarPopG), 
           "decomp" = c("mean_individual", "pop_mean"), 
           "scale" = "chr", "distance" = "genetic")


# 5. ========== MODWT on Physical Scale ==========

# define levels based on longest chromosome
maxLevelP <- max(gnomP[,floor(log2(length(unique(position)))),by=chr][,2])
allColsP <- paste0("d",1:maxLevelP)

# 5.1. ---------- MODWT on Individuals ----------
indModwtP <- gnomP[, modwtAllScales(.SD,variable=ind_minor_parent,lenCol=position,allcols=allColsP), by = .(chr,ID)]

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
setnames(indMeanWavVarP, "V1", "variance")
indMeanWavVarP[, decomp := "mean_individual"]


# 5.2. --------- MODWT on Population Mean ----------

# run wavelet decomp on population mean for each chrom separately
popModwtP <- gnomP[ID==ID[1], modwtAllScales(.SD,variable=pop_mean_minor_parent,lenCol=position,allcols=allColsP), by = chr]
popWavVarP <- popModwtP[,lapply(.SD,wav_var),by=chr]
setnames(popWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
popWavVarP <- melt(popWavVarP, measure.vars = as.character(1:maxLevelP),
                         variable.name = "scale", value.name = "variance")
# save the chromosome-level version before averaging for looking at chrs separately
popWavVarP_Chrs <- popWavVarP

# weighted average over chromosomes
popWavVarP <- merge(popWavVarP, chrWeightsP) 
popMeanWavVarP <- popWavVarP[, weighted.mean(variance, weight), by = scale]
setnames(popMeanWavVarP, "V1", "variance")
popMeanWavVarP[, decomp := "pop_mean"]



# 6. ========== Chromosome-Level Analysis: Physical Scale ==========

# 6.1. ---------- Population Mean ----------

# total physical length as weights for chromosomes 
chrLenP <- gnomP[, max(position), by = .(chr)]
setnames(chrLenP, "V1", "len")
chrLenP[, weight := len/sum(len)]

# chromosome means of population mean minor parent ancestry
chrMeansPopMeanP <- gnomP[, mean(pop_mean_minor_parent), by = .(chr)]
setnames(chrMeansPopMeanP, "V1", "avg_frq")
chrMeansPopMeanP[, weight := chrLenP[,weight]]

# weighted average of population-mean ancestry across chromosomes 
weightedMeanPopMeanP <- chrMeansPopMeanP[, weighted.mean(avg_frq, weight)]

# chromosome-level weighted variance of pop mean ancestry
chrVarPopP <- chrMeansPopMeanP[, sum(weight*(avg_frq - weightedMeanPopMeanP)^2)]


# 6.2. ---------- Individual-Level ----------

# chromosome-level average minor parent ancestry for individuals
chrMeansIndP <- gnomP[, mean(ind_minor_parent), by = .(ID,chr)]
setnames(chrMeansIndP, "V1", "avg_frq")

# calculate genome-wide mean for individuals by weighting chromosomes
chrMeansIndP <- merge(chrMeansIndP, chrLenP) # add chrom weights
chrMeansIndP[, gnomWideMean := weighted.mean(avg_frq,weight), by = ID]

# chromosome-level weighted variance by individual, then average over individuals
chrVarIndP <- chrMeansIndP[, sum(weight*(avg_frq - gnomWideMean)^2), by = ID][,mean(V1)]

# output weighted chromosome-level variance for avg. individual-level signal 
# and population mean 
chrVarP <- data.table("variance" = c(chrVarIndP, chrVarPopP), 
                      "decomp" = c("mean_individual", "pop_mean"), 
                      "scale" = "chr", "distance" = "physical")

# 7. ========== Output Wavelet and Chrom-Level Variance Data ==========

# Combine Wavelet variance tables
wvFinalG <- rbind(indMeanWavVarG, popMeanWavVarG)
wvFinalG[, distance := "genetic"] 
wvFinalP <- rbind(indMeanWavVarP, popMeanWavVarP)
wvFinalP[, distance := "physical"]
wvFinalAll <- rbind(wvFinalP, wvFinalG)

# combine chromosome-level variances
chrVarAll <- rbind(chrVarP, chrVarG)

save(wvFinalAll, chrVarAll, popWavVarG_Chrs, popWavVarP_Chrs, file = paste0("ACUA_",year,"/varDecompAll.RData"))


# 8. ========== Correlation Analysis ==========

# 8.1. ---------- Wavelet Correlation Analysis ----------

# run modwt on recombination rate
rModwtP <- gnomP[ID==ID[1],modwtAllScales(.SD,variable=r_interp,lenCol=position,allcols=allColsP), by = chr]

# wavelet variance for rec rate
rWavVarP <- rModwtP[,lapply(.SD,wav_var),by=chr]
setnames(rWavVarP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
rWavVarP <- melt(rWavVarP, measure.vars = as.character(1:maxLevelP),
                       variable.name = "scale", value.name = "r_variance")

# combine tables with wavelet variances for both signals
rAncWavVar <- merge(rWavVarP, popWavVarP, by = c("chr", "scale"))

# get modwt data tables into more convenient shape (only after calculating wav var)
rModwtP[, position := seq_len(.N), by = chr]
setnames(rModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
rModwtP <- melt(rModwtP, measure.vars = as.character(1:maxLevelP),
                   variable.name = "scale", value.name = "r_coeff")

popModwtP[, position := seq_len(.N), by = chr]
setnames(popModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
popModwtP <- melt(popModwtP, measure.vars = as.character(1:maxLevelP),
                  variable.name = "scale", value.name = "anc_coeff")

rAncModwtP <- merge(popModwtP, rModwtP)

# unweighted wavelet covariance at each scale separately for each chromosome
rAncCov <- rAncModwtP[, sum(r_coeff*anc_coeff, na.rm=T)/length(r_coeff[!is.na(r_coeff)]), by = .(chr, scale)]
setnames(rAncCov, "V1", "covariance")

# merge covariance table with wavelet variance table
rAncCovWavVar <- merge(rAncCov, rAncWavVar, by = c("chr", "scale"))

# Weight covariances and variances
weightedVarCov <- rAncCovWavVar[, lapply(.SD, weighted.mean, w=weight), by = scale, 
                             .SDcols = c("covariance", "variance", "r_variance")]

# Compute final correlation using weights
wvCorFinal <- weightedVarCov[, covariance/(sqrt(variance)*sqrt(r_variance)), by = scale]
setnames(wvCorFinal, "V1", "cor")
#plot(y =wvCorAll$cor, x=1:15)

# 8.2. ---------- Chromosome-Scale correlation ----------

chrMeansRecP <- gnomP[, mean(r_interp), by = chr]
setnames(chrMeansRecP, "V1", "rec")

chrMeansRecAnc <- merge(chrMeansPopMeanP, chrMeansRecP)

# compute weighted means of signals
chrW8MeansRecAnc <- chrMeansRecAnc[, lapply(.SD,weighted.mean), .SDcols = c("avg_frq", "rec")]

# numerator
chrCorNum <- chrMeansRecAnc[, sum( weight*(avg_frq-chrW8MeansRecAnc$avg_frq)*(rec-chrW8MeansRecAnc$rec))/length(avg_frq)]

# weighted st. dev. of two signals
chrCorDenom <- denominator <- chrMeansRecAnc[, sqrt(
  (
    sum(weight*(avg_frq - chrW8MeansRecAnc$avg_frq)^2)/length(avg_frq))*
    (
      sum(weight*(rec - chrW8MeansRecAnc$rec)^2)/length(rec))
  )]

chrCor <- chrCorNum/chrCorDenom

# ========== Save Correlation Analysis Outputs ==========
save(wvCorFinal, chrCor, file = paste0("ACUA_",year,"/rAncCorDecompAll.RData"))
