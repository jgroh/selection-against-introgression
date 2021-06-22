# ===== Inputs =====

library(data.table)
library(waveslim)

# arguments from command line
args <- commandArgs(trailingOnly = TRUE)
metaFile <- args[1]
pop <- args[2]
ancPath <- args[3]
outPath <- args[4]

# run this block only if running locally
#pop <- "RIMMA0366"
#ancPath <- "hmm_anc_interp/genetic/"
#metaFile <- "HILO_MAIZE55_PARV50_meta.txt"
#

# Read files
meta <- fread(metaFile)

# read whichever files present from individuals in focal population
indFiles <- as.list(paste0(ancPath, 
                           intersect(list.files(ancPath), 
                                     paste0(meta[RI_ACCESSION==pop,ID],".txt"))))
# combine individuals from population into one table
gnoms <- rbindlist(lapply(indFiles, fread))

# ===== Define Functions =====

# maximum overlap discrete wavelet transform
modwtAllScales <- function(x,variable,lenCol,allcols){
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

# Unbiased estimator of wavelet variance for MODWT with brick wall boundary condition
wav_var <- function(x){sum(x^2,na.rm=TRUE)/(length(x[!is.na(x)]))} 


# ===== Data transformations =====

# logit transform ancestry (?)
gnoms[freqMex > 1, freqMex := 1] # strangely some values slightly > 1, precision error?

# replace values of zero or 1 by small deviation so logit works
epsilon <- gnoms[freqMex > 0, min(freqMex)]/2
gnoms[, freqMexTr := log(freqMex/(1-freqMex))]
gnoms[freqMex == 0, freqMexTr := log(epsilon/(1-epsilon))]
gnoms[freqMex == 1, freqMexTr := log((1-epsilon)/epsilon)]

# take mean over individuals
meanAnc <- gnoms[, lapply(.SD,mean),.SDcols = c('freqMex','freqMexTr'),by=.(Morgan,chr)]

# mean ancestry in population

totalMeanAnc <- meanAnc[, mean(freqMex)]
# ===== Wavelet Transform =====

# define levels based on longest chromosome
maxLevel <- max(gnoms[,floor(log2(length(unique(Morgan)))),by=chr][,2])
allCols <- paste0("d",1:maxLevel)

# MODWT on individuals
indAncModwt <- gnoms[, modwtAllScales(.SD,variable=freqMexTr,lenCol=Morgan,allcols=allCols),by=.(ID,chr)]

# MODWT on mean ancestry
meanAncModwt <- meanAnc[,modwtAllScales(.SD,variable=freqMexTr,lenCol=Morgan,allcols=allCols), by = chr]

# ===== Wavelet Variance =====

# number of wavelet coefficients per scale on each chromosome
chrWeights <- meanAncModwt[,lapply(.SD,function(x)length(x[!is.na(x)]) ), by = chr]
setnames(chrWeights, old = paste0("d",1:maxLevel), new = as.character(1:maxLevel))
chrWeights <- melt(chrWeights, id.vars = "chr", measure.vars = as.character(1:maxLevel),
                    variable.name = "scale", value.name = "numCoeffs")

# rescale weights to percentages
chrWeights[,weight := numCoeffs/sum(numCoeffs), by=.(scale)]

# ----- wavelet variance: chromosome x individual ----- 
indAncWavVarChrs <- indAncModwt[,lapply(.SD,wav_var),by=.(chr,ID)]

# ----- wavelet variance: individual -----
# weighted average over chromosomes within each individual
setnames(indAncWavVarChrs, paste0("d",1:maxLevel), as.character(1:maxLevel))
indAncWavVarChrs <- melt(indAncWavVarChrs, measure.vars = as.character(1:maxLevel),
                          variable.name = "scale", value.name = "anc_variance")
# weighted average over chromosomes
indAncWavVar <- merge(indAncWavVarChrs, chrWeights) 
indAncWavVar <- indAncWavVar[, weighted.mean(anc_variance, weight), by=.(ID,scale)]
setnames(indAncWavVar, "V1", "anc_variance")
# average over individuals
indMeanAncWavVar <- indAncWavVar[, mean(anc_variance), by = scale]
setnames(indMeanAncWavVar, "V1", "anc_variance")
indMeanAncWavVar[,decomp :="individual"]


# ----- wavelet variance: mean -----
meanAncWavVarChrs <- meanAncModwt[,lapply(.SD,wav_var),by=chr]
setnames(meanAncWavVarChrs, paste0("d",1:maxLevel), as.character(1:maxLevel))
meanAncWavVarChrs <- melt(meanAncWavVarChrs, measure.vars = as.character(1:maxLevel),
                       variable.name = "scale", value.name = "anc_variance")
# weighted average over chromosomes
meanAncWavVar <- merge(meanAncWavVarChrs, chrWeights) 
meanAncWavVar <- meanAncWavVar[, weighted.mean(anc_variance, weight), by = scale]
setnames(meanAncWavVar, "V1", "anc_variance")
meanAncWavVar[,decomp :="mean"]


# ----- combine wav var tables -----
allAncWavVar <- rbind(meanAncWavVar, indMeanAncWavVar)



# ===== Chromosome Variance =====
# new weights: total genetic length

chrMorgans <- gnoms[ID==ID[1], max(Morgan), by=chr]
setnames(chrMorgans, "V1","MorganLength")

# ----- individual-level -----
indChrMeanAnc <- gnoms[, mean(freqMex), by = .(chr,ID)]
setnames(indChrMeanAnc,"V1","meanFreqMex")
indChrMeanAnc <- merge(indChrMeanAnc,chrMorgans)
indChrMeanAnc[, weight := MorganLength/sum(MorganLength), by = ID]

# weighted total mean for individuals
indChrMeanAnc[, indMean := weighted.mean(meanFreqMex,w=MorganLength),by=ID]
# chrom-level weighted variance
indChrVarAnc <- data.table(anc_variance=indChrMeanAnc[, sum(weight*(meanFreqMex-indMean)^2),by=ID][,mean(V1)],
                           decomp = "individual",
                           scale = "chrom")

# ----- mean ancestry ------
meanChrMeanAnc <- meanAnc[, mean(freqMex),by=chr]
setnames(meanChrMeanAnc,"V1","meanFreqMex")
meanChrMeanAnc <- merge(meanChrMeanAnc,chrMorgans)
meanChrMeanAnc[, weight := MorganLength/sum(MorganLength)]

#chrom-level weighted variance
meanChrVarAnc <- data.table(anc_variance=meanChrMeanAnc[, sum(weight*(meanFreqMex-totalMeanAnc)^2)],
                            decomp = "mean",
                            scale = "chrom")

# ===== Final Output =====
allAncVarDecomp <- rbind(indChrVarAnc, meanChrVarAnc, allAncWavVar)

# Add metadata to table for output
allAncVarDecomp[,popID := pop]
allAncVarDecomp[,species := meta[RI_ACCESSION==pop,zea][1]]
allAncVarDecomp[,locality := meta[RI_ACCESSION==pop,LOCALITY][1]]
allAncVarDecomp[,totalMean := totalMeanAnc]


save(allAncVarDecomp, file=outPath)

