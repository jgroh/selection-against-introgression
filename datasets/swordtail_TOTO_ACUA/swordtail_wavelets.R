# script combines scaffold data from each year and runs wavelets
library(data.table)
library(tools)
library(waveslim)

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]

scaffFiles <- dir(paste0("ACUA_",year),full.names=T)
scaffs <- basename(file_path_sans_ext(scaffFiles))

# each file has the same object name for the scaffold 
# so we load into separate environments
loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 

# combine data into one file for the year
gnom <- rbindlist(lapply(as.list(scaffFiles), 
                  function(x){loadFrom(x, "chrom1.interp")}))


# MODWT ----------------
# run wavelet decomposition on individuals, separately for each chromosome
# since chromosomes have different levels, write a function
max.level <- max(gnom[,floor(log2(length(unique(gen_pos)))),by=chr][,2])
all.cols <- paste0("d",1:max.level)

modwtAllScales <- function(x,variable){
  dt <- setDT( 
    x[,modwt(variable,"haar",n.levels=floor(log2(length(unique(gen_pos)))))] 
    )
  # remove smooth coeff column 
  smooth.col <- grep("s",names(dt))
  dt <- dt[,-smooth.col,with=F] 
  
  # add empty detail columns for higher scales not present
  naCols <- setdiff(all.cols,names(dt))
  if(length(naCols) > 0){ dt[,(naCols) := as.numeric(NA)] }
  return(dt)
}
ind.modwt <- gnom[, modwtAllScales(.SD,variable=ind_frq_interp), by = .(chr,ID)]
      

# compute wavelet variance
wav_var <- function(x){sum(x^2)/(length(x))}
ind.wav.var <- ind.modwt[,lapply(.SD,wav_var), by = .(chr,ID)]
setnames(ind.wav.var, old = paste0("d",1:max.level), new = as.character(1:max.level))
ind.wav.var <- melt(ind.wav.var, id.vars = c("ID","chr"), 
                    variable.name = "scale",
                    value.name = "variance")

# compute mean across individuals for each chromosome separately
ind.mean.wav.var <- ind.wav.var[, mean(variance), by = .(chr,scale)]
setnames(ind.mean.wav.var, "V1", "variance")

# run wavelet decomp on population mean for each chrom separately
pop.mean.modwt <- gnom[,modwtAllScales(.SD,variable=pop_mean_interp), by = chr]
pop.mean.wav.var <- pop.mean.modwt[,lapply(.SD,wav_var),by=chr]
setnames(pop.mean.wav.var, paste0("d",1:max.level), as.character(1:max.level))
pop.mean.wav.var <- melt(pop.mean.wav.var, measure.vars = as.character(1:max.level),
                         variable.name = "scale", value.name = "variance")
pop.mean.wav.var[,scale:=as.numeric(scale)]

# combine mean of individuals and pop mean wav var tables
ind.mean.wav.var[,decomp :="mean_individual"]
pop.mean.wav.var[,decomp :="pop_mean"]

wvFinal <- rbind(ind.mean.wav.var,pop.mean.wav.var)
save(wvFinal, file = paste0("ACUA_",year,"/wvFinal.RData"))

# ggplot(ind.wav.var, aes(x=as.numeric(scale), y=variance)) + 
#   theme(legend.position = "none") + 
#   geom_line(aes(group = ID, color = ID, alpha = 0.8)) +
#   geom_line(data = pop.mean.wav.var, size = 1) + 
#   geom_line(data = ind_mean_wav_var, color = "red", size = 1) + 
#   facet_wrap(~chr)+
#   #geom_ribbon(data = pop.mean.wv, aes(ymin=lower,ymax=upper)) +
#   labs(x = "Scale log2 (N x Morgans)", y = "Wavelet variance")#   +
#  #geom_line(data = pop.mean.roll.wav.var, color = "blue", size = 1) 
