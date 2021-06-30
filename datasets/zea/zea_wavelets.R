# ===== Inputs =====

library(data.table)
library(waveslim)
library(magrittr)
library(psych)

# arguments from command line
args <- commandArgs(trailingOnly = TRUE)
metaFile <- args[1]
pop <- args[2]
ancPathG <- args[3]
ancPathP <- args[4]
recFile <- args[5]
outPath <- args[6]

# run this block only if running locally
# metaFile <- "HILO_MAIZE55_PARV50_meta.txt"
# pop <- "RIMME0035" # mexicana 
# pop <- "RIMMA0366" # maize
# ancPathG <- "hmm_anc_interp/genetic/"
# ancPathP <- "hmm_anc_interp/physical/"
# recFile <- "zea_1kb_rec.bed"

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
#qqnorm(rec[chr==1,cmTr])


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
vars <- c("gnoms","meanAnc","totalMeanAnc")
d <- setDT(lapply(vars, function(x)paste0(x,c("G","P"))))
setnames(d,vars)

for(i in 1:2){
  # mean over individuals
  assign(d[i,meanAnc], 
         get(d[i,gnoms])[, 
                         lapply(.SD,mean),
                         .SDcols = c('freqMex','freqMexTr'),
                         by=.(position,chr)])
  
  # total mean ancestry in population
  assign( d[i,totalMeanAnc], get(d[i,meanAnc])[, lapply(.SD,mean),
                                               .SDcols = c('freqMex','freqMexTr')] )
}

#qqnorm(meanAncG$freqMexTr)

# ===== Wavelet Transforms =====
vars <- c("gnoms","maxLevel","allCols","indAncModwt","meanAnc","meanAncModwt")
d <- setDT(lapply(vars, function(x)paste0(x,c("G","P"))));setnames(d,vars)

for(i in 1:2){
  # maximum level of wavelet decomposition
  assign(d[i,maxLevel],
         max( get(d[i,gnoms]) [ID==ID[1],
                              .(level = floor(log2(length(position)))),
                              by=chr][,level] ) )
  # these are the names waveslim gives to scales
  assign(d[i,allCols], paste0("d",1:get(d[i,maxLevel])))
  
  # MODWT on individuals
  assign(d[i,indAncModwt], 
         get(d[i,gnoms])[, modwtAllScales(.SD,
                                          variable=freqMexTr,
                                          lenCol=position,
                                          allcols=get(d[i,allCols])), by=.(ID,chr)])
  
  # MODWT on mean ancestry
  assign(d[i,meanAncModwt],
         get(d[i,meanAnc])[, modwtAllScales(.SD,
                                            variable=freqMexTr,
                                            lenCol=position,
                                            allcols=get(d[i,allCols])), by = chr])
}

# MODWT of recombination rate
recModwtP <- rec[, modwtAllScales(.SD,variable=cmTr,lenCol=position,allcols=allColsP), by=chr]


# ===== Wavelet Variances =====

# perform calculations for both genetic and physical scales
vars <- c("gnoms","indAncTotalVar","meanAncTotalVar","chrWeights","allCols","meanAncModwt",
          "indAncModwt","indAncWVChrs","indAncWV",
          "indMeanAncWV","meanAncWVChrs","meanAncWV","allAncWV")
d <- setDT(lapply(vars, function(x)paste0(x,c("G","P"))))
setnames(d,vars)

for(i in 1:2){ 
  
  # ----- total ancestry variance -----

  # within individuals
  assign(d[i,indAncTotalVar], gnomsG[, .(totalVar = var(freqMexTr)), by=ID])
  
  # mean ancestry
  assign(d[i,meanAncTotalVar], meanAncG[, var(freqMexTr)])
  
  # ----- chromosome weights by scale -----
  assign(d[i,chrWeights], 
         get(d[i,meanAncModwt])[,
                     lapply(.SD,function(x)length(x[!is.na(x)])),
                     by = chr] %>%
           melt(., id.vars = "chr", measure.vars = get(d[i,allCols]),
                variable.name = "scale", value.name = "numCoeffs"))
  get(d[i,chrWeights])[, scale := gsub("d","",scale)]
  get(d[i,chrWeights])[, weight := numCoeffs/sum(numCoeffs), by = scale]
  
  # ----- wavelet variance: chromosome x individual ----- 
  assign(d[i,indAncWVChrs],
         get(d[i,indAncModwt]) [,lapply(.SD,wav_var),by=.(chr,ID)] %>%
           melt(., measure.vars=get(d[i,allCols]),
                variable.name="scale", value.name="anc_variance"))
  get(d[i,indAncWVChrs]) [, scale:=gsub("d","",scale)]
  
  # ----- wavelet variance: individual -----
  
  # convention: if scale not present on chromosome, in order to get proportion of total variance, 
  # set the wavelet variance for that scale to zero
  get(d[i,indAncWVChrs])[is.nan(anc_variance), anc_variance := 0]
 
  #  weighted average over chromosomes 
  assign(d[i,indAncWV], 
         merge(get(d[i,indAncWVChrs]), get(d[i,chrWeights]), by = c("chr","scale"))  %>%
           .[, .(anc_variance = weighted.mean(anc_variance, weight)), by=.(ID,scale)]
         )
  
  # merge with total ancestry variance per individual
  assign(d[i,indAncWV], merge(get(d[i,indAncWV]), get(d[i,indAncTotalVar])))
  
  
  # compute average proportion of variance and average variance
  get(d[i,indAncWV])[, ancPropVar := anc_variance/totalVar]
  a <- get(d[i,indAncWV])[, .(ancPropVar = mean(ancPropVar)), by = scale]
  b <- get(d[i,indAncWV])[, .(anc_variance = mean(anc_variance)), by = scale]
  
  assign(d[i,indMeanAncWV], merge(a,b))
  
  get(d[i,indMeanAncWV])[,decomp :="individual"]
  
  
  # ----- wavelet variance: mean ancestry -----
  assign(d[i,meanAncWVChrs],
         get(d[i,meanAncModwt])[,lapply(.SD,wav_var),by=chr] %>%
           melt(., measure.vars = get(d[i,allCols]),
                variable.name = "scale", value.name = "anc_variance"))
  get(d[i,meanAncWVChrs])[, scale:=gsub("d","",scale)]
  
  # convention: if scale not present on chromosome, in order to get proportion of total variance, 
  # set the wavelet variance for that scale to zero
  get(d[i,meanAncWVChrs])[is.nan(anc_variance), anc_variance := 0]
  
  # ----- weighted average of wavelet variances over chromosome
  assign(d[i,meanAncWV],
         merge(get(d[i,meanAncWVChrs]), get(d[i,chrWeights])) %>%
           .[, .(anc_variance=weighted.mean(anc_variance, weight)), by = scale])
  get(d[i,meanAncWV])[, decomp := "mean"]
  
  # ----- proportion of total variance of mean ancestry along sequence -----
  get(d[i,meanAncWV])[, ancPropVar := anc_variance/get(d[i,meanAncTotalVar])]
  
  # ----- combine ind and mean wav var tables -----
  assign(d[i, allAncWV], 
         rbind(get(d[i,meanAncWV]), get(d[i,indMeanAncWV])))
  
}

allAncWVG[,distance := 'genetic']
allAncWVP[,distance := 'physical']

allAncWV <- rbind(allAncWVG,allAncWVP)


# ----- Wavelet variance of recombination rate -----

recWVChrs <- recModwtP[,lapply(.SD,wav_var),by=chr] %>%
         melt(., measure.vars = allColsP,
              variable.name = "scale", value.name = "rec_variance")
recWVChrs[, scale := gsub("d","",scale)]


# -- weighted average over chromosomes -- 
recWV <- merge(recWVChrs, chrWeightsP) %>%
  .[, .(rec_variance=weighted.mean(rec_variance, weight)), by = scale]
recWV[, decomp := "mean"]
recWV[, distance := "physical"]

allVarsWV <- merge(allAncWV, recWV, all=T)


# ===== Chromosome-Level Variances =====

vars <- c("gnoms","chrLen","indChrMeanAnc",
          "indAncChrVar", "meanAnc","meanAncChrMean",
          "meanAncChrVar","indAncTotalVar","meanAncTotalVar","totalMeanAnc")
d <- setDT(lapply(vars, function(x)paste0(x,c("G","P"))))
setnames(d,vars)

for(i in 1:2){
  # chromosome weights by total length (Morgans or base pairs)
  assign(d[i,chrLen], get(d[i,gnoms])[ID==ID[1],
                                      .(length = max(position)),
                                      by=chr])
  get(d[i,chrLen])[, weight := length/sum(length)]
  
  # ----- chromosome-level variance in individual ancestry -----
  # chromosome x individual ancestry means
  assign(d[i,indChrMeanAnc], 
         merge( get(d[i,gnoms])[, .(meanFreqMexTr = mean(freqMexTr)),by=.(chr,ID)],
                get(d[i,chrLen]),
                by="chr"))
  
  # weighted total ancestry mean for individuals
  get(d[i,indChrMeanAnc])[, indMean := weighted.mean(meanFreqMexTr,w=weight),by=ID]
  
  
  # chrom-level weighted variance within individuals
  assign(d[i,indAncChrVar], 
         merge(
           get(d[i,indChrMeanAnc])[, .(anc_variance = sum(weight*(meanFreqMexTr-indMean)^2)),by=ID],
           get(d[i,indAncTotalVar])))
  
  # get average variance and average proportion of variance across individuals
  assign(d[i,indAncChrVar], get(d[i,indAncChrVar])[,ancPropVar:=anc_variance/totalVar][, lapply(.SD,mean), 
                                                                                         .SDcols = c("anc_variance","ancPropVar")])
  get(d[i,indAncChrVar])[,`:=`(scale="chrom",decomp="individual")]
  
  # ----- chromosome-level variance in mean ancestry -----
  
  # chromosome-level ancestry means
  assign(d[i,meanAncChrMean], 
         merge(
           get(d[i,meanAnc])[, .(meanFreqMexTr = mean(freqMexTr)),by=chr],
           get(d[i,chrLen]),by="chr"))
  
  # chrom-level weighted variance
  assign(d[i,meanAncChrVar],
         data.table(anc_variance=get(d[i,meanAncChrMean])[,sum(weight*(meanFreqMexTr-get(d[i,totalMeanAnc])[,freqMexTr])^2)],
                    decomp = "mean",
                    scale = "chrom"))
  
  # calculate proportion of variance
  get(d[i,meanAncChrVar])[, ancPropVar := anc_variance/get(d[i,meanAncTotalVar])]
}

# add distance and combine
lapply(list(indAncChrVarG,meanAncChrVarG), function(x){
  x[,distance := 'genetic']
  })
lapply(list(indAncChrVarP,meanAncChrVarP), function(x){
  x[,distance := 'physical']
})

# ----- chromosome-level variance in recombination -----
meanRecChrs <- merge(rec[, .(meanCmTr = mean(cmTr)),by=chr], chrLenP, by = 'chr')

totalMeanRec <- meanRecChrs[, weighted.mean(meanCmTr,weight)]

recVarChrs <- data.table(
  rec_variance=meanRecChrs[,sum(weight*(meanCmTr-totalMeanRec)^2)],
  decomp = "mean",
  scale = "chrom",
  distance = "physical")

# ===== Variance Output =====

allVarDecomp <- rbind(merge(recVarChrs, meanAncChrVarP),
                      meanAncChrVarG[, rec_variance := NA],
                      allVarsWV)
        
# Add metadata to table for output
allVarDecomp[,popID := pop]
allVarDecomp[,species := meta[RI_ACCESSION==pop,zea][1]]
allVarDecomp[,locality := meta[RI_ACCESSION==pop,LOCALITY][1]]
allVarDecomp[,totalMean := totalMeanAncP[,freqMex]]


# ===== Correlations =====

# reformat modwt tables
meanAncModwtP[, position := seq_len(.N), by = chr]
setnames(meanAncModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
meanAncModwtP2 <- melt(meanAncModwtP, measure.vars = as.character(1:maxLevelP),
                   variable.name = "scale", value.name = "ancestry_coeff")

recModwtP[, position := seq_len(.N), by = chr]
setnames(recModwtP, paste0("d",1:maxLevelP), as.character(1:maxLevelP))
recModwtP2 <- melt(recModwtP, measure.vars = as.character(1:maxLevelP),
                  variable.name = "scale", value.name = "rec_coeff")

# combine tables
allVarsModwtP <- merge(meanAncModwtP2, recModwtP2)


# correlations:
# 1. compute for scale x chromosome
# 2. Fisher z transform
# 3. average over chromosomes
# 4. back-transform

recAncCor <- allVarsModwtP[, .(z = fisherz(cor(ancestry_coeff,
                                               rec_coeff,
                                               use="pairwise.complete.obs"))),
                           by=.(scale,chr)][,
                                            .(z=mean(z,na.rm=T)),
                                            by=scale][, 
                                                      rec_anc_cor:=fisherz2r(z)][,z:=NULL][]

# ===== Chromosome-scale correlation =====
chrSignalMeans <- merge(meanAncChrMeanP, meanRecChrs, by = c("chr","weight"))


m1 <- totalMeanAncP[,freqMexTr]
m2 <- totalMeanRec
v1 <- allVarDecomp[scale=="chrom"&distance=="physical",anc_variance]
v2 <- allVarDecomp[scale=="chrom"&distance=="physical",rec_variance]

recAncCor <- rbind(data.table(scale = "chrom",
                              rec_anc_cor=chrSignalMeans[, sum(weight*(meanFreqMexTr-m1)*(meanCmTr-m2))/
                                                           sqrt(v1*v2)
                                                         ]),
                   recAncCor)



# ===== Final Output =====
save(allVarDecomp, recAncCor, file = outPath)






