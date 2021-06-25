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
#ancPathG <- "hmm_anc_interp/genetic/"
#ancPathP <- "hmm_anc_interp/physical/"
#metaFile <- "HILO_MAIZE55_PARV50_meta.txt"
#recFile <- "zea_1kb_rec.bed"

# Read files
meta <- fread(metaFile)
rec <- fread(recFile, col.names = c('chr','start','end','cM'))

# ----- read gnom files in genetic and physical distance ----
d <- data.table(x=c("ancPathG","ancPathP"), y=c("gnomsG","gnomsP"))

for(i in 1:2){
  # read whichever files present from individuals in focal population
  indFiles <- as.list(paste0(get(d[i,x]),
                             intersect(list.files(get(d[i,x])), 
                                       paste0(meta[RI_ACCESSION==pop,ID],".txt"))))
  # combine individuals from population into one table
  assign(d[i,y], rbindlist(lapply(indFiles, fread)))
}

# merge ancestry data and genomic features data
rec[, position := start + 500] # this is the midpoint of the 1kb intervals where I interpolate ancestry
gnomsP <- merge(gnomsP, rec, by = c("chr", "position"))


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

# ----- log transform recombination -----
rec[, cmTr := log10(cM)]

# ----- logit transform ancestry -----
lapply(list(gnomsG,gnomsP),function(x){
  x[freqMex > 1, freqMex := 1] # strangely some values that are 1 evaluate to > 1
})

# replace values of zero or 1 by small deviation so logit works
lapply(list(gnomsG,gnomsP),function(x){
  epsilon <- x[freqMex > 0, min(freqMex)]/2
  x[, freqMexTr := log(freqMex/(1-freqMex))]
  x[freqMex == 0, freqMexTr := log(epsilon/(1-epsilon))]
  x[freqMex == 1, freqMexTr := log((1-epsilon)/epsilon)]
})

# ----- mean ancestry -----
# mean over individuals
meanAncG <- gnomsG[, lapply(.SD,mean),.SDcols = c('freqMex','freqMexTr'),by=.(Morgan,chr)]
meanAncP <- gnomsP[, lapply(.SD,mean),.SDcols = c('freqMex','freqMexTr'),by=.(position,chr)]

# total mean ancestry in population
totalMeanAncG <- meanAncG[, mean(freqMex)]
totalMeanAncP <- meanAncP[, mean(freqMex)]


# ===== Wavelet Transforms =====

# assign max level and column names for scales
maxLevelG <- max( gnomsG[ID==ID[1], 
            .(level = floor(log2(length(freqMex)))),
            by=chr][,level])
maxLevelP <- max( gnomsP[ID==ID[1], 
                         .(level = floor(log2(length(freqMex)))),
                         by=chr][,level])
allColsG <- paste0("d",1:maxLevelG) # these are the names waveslim gives to scales
allColsP <- paste0("d",1:maxLevelP)


# MODWT on individuals
indAncModwtG <- gnomsG[, modwtAllScales(.SD,variable=freqMexTr,lenCol=Morgan,allcols=allColsG),by=.(ID,chr)]
indAncModwtP <- gnomsP[, modwtAllScales(.SD,variable=freqMexTr,lenCol=position,allcols=allColsP),by=.(ID,chr)]

# MODWT on mean ancestry
meanAncModwtG <- meanAncG[,modwtAllScales(.SD,variable=freqMexTr,lenCol=Morgan,allcols=allColsG), by = chr]
meanAncModwtP <- meanAncP[,modwtAllScales(.SD,variable=freqMexTr,lenCol=position,allcols=allColsP), by = chr]

# MODWT of recombination rate
recModwtP <- rec[, modwtAllScales(.SD,variable=cmTr,lenCol=position,allcols=allColsP), by=chr]


# ===== Wavelet Variances =====

# number of wavelet coefficients per scale on each chromosome
d <- data.table(A=c("chrWeightsG","chrWeightsP"),
                B=c("meanAncModwtG","meanAncModwtP"),
                C=c("allColsG","allColsP"))

for(i in 1:2){
  assign(d[i,A], 
         get(d[i,B])[,
                     lapply(.SD,function(x)length(x[!is.na(x)])),
                     by = chr] %>%
           melt(., id.vars = "chr", measure.vars = get(d[i,C]),
                variable.name = "scale", value.name = "numCoeffs"))
  get(d[i,A])[, scale := gsub("d","",scale)]
  get(d[i,A])[, weight := numCoeffs/sum(numCoeffs), by = scale]
}


# ----- wavelet variance: chromosome x individual ----- 
indAncWavVarChrs <- indAncModwt[,lapply(.SD,wav_var),by=.(chr,ID)]

# ----- wavelet variance: individual -----
setnames(indAncWavVarChrs, paste0("d",1:maxLevel), as.character(1:maxLevel))
indAncWavVarChrs <- melt(indAncWavVarChrs, measure.vars = as.character(1:maxLevel),
                          variable.name = "scale", value.name = "anc_variance")
# weighted average over chromosomes
indAncWavVar <- merge(indAncWavVarChrs, chrWeights) 
indAncWavVar <- indAncWavVar[, weighted.mean(anc_variance, weight), by=.(ID,scale)]
setnames(indAncWavVar, "V1", "anc_variance")
# now average over individuals
indMeanAncWavVar <- indAncWavVar[, mean(anc_variance), by = scale]
setnames(indMeanAncWavVar, "V1", "anc_variance")
indMeanAncWavVar[,decomp :="individual"]


# ----- wavelet variance: mean ancestry -----
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
# chrom-level weighted variance, then averaged over individuals
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

