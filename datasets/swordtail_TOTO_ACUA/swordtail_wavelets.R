# script combines scaffold data from each year and runs wavelets
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

# MODWT ----------------
# run wavelet decomposition on individuals, separately for each chromosome
# since chromosomes have different levels, write a function
max.level <- max(gnomG[,floor(log2(length(unique(Morgan)))),by=chr][,2])
all.cols <- paste0("d",1:max.level)

modwtAllScales <- function(x,variable){
  dt <- setDT( 
    x[,brick.wall(wf="haar", x=modwt(variable,"haar",n.levels=floor(log2(length(unique(Morgan))))))] 
    )
  # remove smooth coeff column 
  smooth.col <- grep("s",names(dt))
  dt <- dt[,-smooth.col,with=F] 
  
  # add empty detail columns for higher scales not present
  naCols <- setdiff(all.cols,names(dt))
  if(length(naCols) > 0){ dt[,(naCols) := as.numeric(NA)] }
  return(dt)
}
ind.modwt <- gnomG[, modwtAllScales(.SD,variable=ind_frq_interp), by = .(chr,ID)]

# get number of wavelet coefficients per scale on each chromosome
numCoeff <- function(x){length(x[!is.na(x)])} # exclude NAs in length as these are boundary coefficients
chrWeights <- ind.modwt[ID==ID[1],lapply(.SD,numCoeff), by = .(chr,ID)]
setnames(chrWeights, old = paste0("d",1:max.level), new = as.character(1:max.level))
chrWeights <- melt(chrWeights[,-"ID"], id.vars = "chr", measure.vars = as.character(1:max.level),
                         variable.name = "scale", value.name = "numCoeffs")
# obtain final weights for each chromosome by scale
chrWeights[,weight := numCoeffs/sum(numCoeffs), by=.(scale)]

# compute wavelet variance
wav_var <- function(x){sum(x^2,na.rm=TRUE)/(length(x[!is.na(x)]))} # exclude NAs in length as these are boundary coefficients
ind.wav.var <- ind.modwt[,lapply(.SD,wav_var), by = .(chr,ID)]
setnames(ind.wav.var, old = paste0("d",1:max.level), new = as.character(1:max.level))
ind.wav.var <- melt(ind.wav.var, id.vars = c("ID","chr"), 
                    variable.name = "scale",
                    value.name = "variance")

# compute mean across individuals for each chromosome separately
ind.mean.wav.var <- ind.wav.var[, mean(variance,na.rm=T), by = .(chr,scale)]
setnames(ind.mean.wav.var, "V1", "variance")

# run wavelet decomp on population mean for each chrom separately
pop.mean.modwt <- gnomG[,modwtAllScales(.SD,variable=pop_mean_interp), by = chr]
pop.mean.wav.var <- pop.mean.modwt[,lapply(.SD,wav_var),by=chr]
setnames(pop.mean.wav.var, paste0("d",1:max.level), as.character(1:max.level))
pop.mean.wav.var <- melt(pop.mean.wav.var, measure.vars = as.character(1:max.level),
                         variable.name = "scale", value.name = "variance")
pop.mean.wav.var[,scale:=as.numeric(scale)]

# combine mean of individuals and pop mean wav var tables
ind.mean.wav.var[,decomp :="mean_individual"]
pop.mean.wav.var[,decomp :="pop_mean"]

wvFinal <- rbind(ind.mean.wav.var, pop.mean.wav.var)
wvFinal <- merge(wvFinal, chrWeights)[,-"numCoeffs"][]

# output chromosome means of pop mean (to account for total variance)
chrMeans <- gnomP[, mean(1-pop_mean_interp), by = .(chr)]
setnames(chrMeans, "V1", "avg_frq")
chrVarPop <- chrMeans[, mean(avg_frq^2)-mean(avg_frq)^2]

# chromosome level variance for individuals
chrMeansInd <- gnomP[, mean(1-ind_frq_interp), by = .(ID,chr)]
setnames(chrMeansInd, "V1", "avg_frq")
# Take mean of chromosome-level variance over individuals 
chrVarInd <- chrMeansInd[, mean(avg_frq^2) - mean(avg_frq)^2, by = ID][, mean(V1)]

chrVar <- data.table("chrVar" = c(chrVarInd, chrVarPop), 
           "decomp" = c("individual", "pop_mean"))

save(wvFinal, chrVar, file = paste0("ACUA_",year,"/wvFinal.RData"))


# Correlation Analysis ---------

max.level <- max(gnomP[,floor(log2(length(unique(position)))),by=chr][,2])
all.levels <- paste0("d",1:max.level)


wavcor <- function(data,xvar,yvar){
  x <- unlist(data[,xvar, with=F])
  y <- unlist(data[,yvar, with=F])
  modwt.x <- brick.wall(wf="haar",x=modwt(x,"haar",n.levels=floor(log2(length(x)))))
  modwt.y <- brick.wall(wf="haar",x=modwt(y,"haar",n.levels=floor(log2(length(y)))))

  wc <- as.data.table(
    t( wave.correlation(modwt.x, modwt.y, N=length(unique(data$position)))[,"wavecor"] ) )

  # remove smooth coeff col
  smooth.col <- grep("s",names(wc))
  wc <- wc[,-smooth.col,with=F] 
  
  # add empty detail columns for higher scales not present
  naCols <- setdiff(all.levels,names(wc))
  if(length(naCols) > 0){ wc[,(naCols) := as.numeric(NA)] }
  return(wc)
}

gnomP[,ind_minor_parent := 1 - ind_frq_interp]
gnomP[,pop_mean_minor_parent := 1 - pop_mean_interp]

# run wavelet correlation on individuals
ind_wavcor <- gnomP[, wavcor(.SD, xvar="r_interp", yvar="ind_minor_parent"), by=.(ID,chr)]

# ind_wavcor
# note individuals either have strong negative or strong positive correlations at largest scale
# would log be more informative? 

# compute mean across individuals for each chromosome separately
ind_mean_wavcor <- ind_wavcor[, lapply(.SD,mean), by = chr, .SDcols = all.cols]

# run wavelet correlation on population mean
pop_mean_wavcor <- gnomP[ID==ID[1],wavcor(.SD,xvar="r_interp",yvar="pop_mean_minor_parent"), by = chr]

# combine mean of individuals and pop mean wav cor
ind_mean_wavcor[,signal :="individual"]
pop_mean_wavcor[,signal :="pop_mean"]

wavcorFinal <- rbind(ind_mean_wavcor, pop_mean_wavcor)
setnames(wavcorFinal, paste0("d",1:max.level), as.character(1:max.level))
wavcorFinal <- melt(wavcorFinal, measure.vars = as.character(1:max.level),
                         variable.name = "scale", value.name = "variance")

# and add chromosome weights to output
wavcorFinal <- merge(wavcorFinal, chrWeights)[,-"numCoeffs"][]
save(wavcorFinal, file = paste0("ACUA_",year,"/wavcorFinal.RData"))



# ggplot(ind.wav.var, aes(x=as.numeric(scale), y=variance)) + 
#   theme(legend.position = "none") + 
#   geom_line(aes(group = ID, color = ID, alpha = 0.8)) +
#   geom_line(data = pop.mean.wav.var, size = 1) + 
#   geom_line(data = ind_mean_wav_var, color = "red", size = 1) + 
#   facet_wrap(~chr)+
#   #geom_ribbon(data = pop.mean.wv, aes(ymin=lower,ymax=upper)) +
#   labs(x = "Scale log2 (N x Morgans)", y = "Wavelet variance")#   +
#  #geom_line(data = pop.mean.roll.wav.var, color = "blue", size = 1) 
