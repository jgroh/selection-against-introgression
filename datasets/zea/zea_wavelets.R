library(data.table)
library(waveslim)

args <- commandArgs(trailingOnly = TRUE)
metaFile <- args[1]
pop <- args[2]
ancPath <- args[3]
outPath <- args[4]

meta <- fread(metaFile)

# read whichever files present from individuals in focal population
indFiles <- as.list(paste0(ancPath, 
                           intersect(list.files(ancPath), 
                                     paste0(meta[RI_ACCESSION==pop,ID],".txt"))))
# combine individuals
gnoms <- rbindlist(lapply(indFiles, fread))

# logit transform ancestry
gnoms[freqMex > 1, freqMex := 1] # strangely some values slightly > 1, precision error?

# replace values of zero or 1 by small deviation so logit works
epsilon <- gnoms[freqMex > 0, min(freqMex)]/2
gnoms[, freqMexTr := log(freqMex/(1-freqMex))]
gnoms[freqMex == 0, freqMexTr := log(epsilon/(1-epsilon))]
gnoms[freqMex == 1, freqMexTr := log((1-epsilon)/epsilon)]

# take mean over individuals
meanAnc <- gnoms[, lapply(.SD,mean),.SDcols = c('freqMex','freqMexTr'),by=.(Morgan,chr)]

# To compute modwt on variable in a data table:
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

# define levels based on longest chromosome
maxLevel <- max(gnoms[,floor(log2(length(unique(Morgan)))),by=chr][,2])
allCols <- paste0("d",1:maxLevel)

# MODWT on mean ancestry
meanAncModwt <- meanAnc[,modwtAllScales(.SD,variable=freqMexTr,lenCol=Morgan,allcols=allCols), by = chr]

# number of wavelet coefficients per scale on each chromosome
chrWeights <- meanAncModwt[,lapply(.SD,function(x)length(x[!is.na(x)]) ), by = chr]
setnames(chrWeights, old = paste0("d",1:maxLevel), new = as.character(1:maxLevel))
chrWeights <- melt(chrWeights, id.vars = "chr", measure.vars = as.character(1:maxLevel),
                    variable.name = "scale", value.name = "numCoeffs")

# rescale weights to percentages
chrWeights[,weight := numCoeffs/sum(numCoeffs), by=.(scale)]

# Unbiased estimator of wavelet variance for MODWT with brick wall boundary condition
wav_var <- function(x){sum(x^2,na.rm=TRUE)/(length(x[!is.na(x)]))} 

# wavelet variance by chromosome
meanAncWavVarChrs <- meanAncModwt[,lapply(.SD,wav_var),by=chr]
setnames(meanAncWavVarChrs, paste0("d",1:maxLevel), as.character(1:maxLevel))
meanAncWavVarChrs <- melt(meanAncWavVarChrs, measure.vars = as.character(1:maxLevel),
                       variable.name = "scale", value.name = "anc_variance")

# weighted average of wavelet variance over chromosomes
meanAncWavVar <- merge(meanAncWavVarChrs, chrWeights) 
meanAncWavVar <- meanAncWavVar[, weighted.mean(anc_variance, weight), by = scale]
setnames(meanAncWavVar, "V1", "anc_variance")
meanAncWavVar[,decomp :="mean"]
meanAncWavVar[,popID := pop]
meanAncWavVar[,species := meta[RI_ACCESSION==pop,zea][1]]
meanAncWavVar[,locality := meta[RI_ACCESSION==pop,LOCALITY][1]]


save(meanAncWavVar, file=outPath)

