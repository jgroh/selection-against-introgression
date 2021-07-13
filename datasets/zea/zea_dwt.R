# ===== Inputs =====

library(data.table)
library(waveslim)
library(magrittr)
library(psych)
source("../wavelet_functions.R")

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
# pop <- "RIMMA0366" # maize
# ancPathG <- "hmm_anc_interp/genetic/"
# ancPathP <- "hmm_anc_interp/physical/"
# recFile <- "zea_5kb_rec.bed"


# Read files
meta <- fread(metaFile)
rec <- fread(recFile, col.names = c('chr','start','end','cM'))


d <- data.table(units = c("genetic", "physical"), 
                gnoms = c("gnomsG", "gnomsP"),
                path = c(ancPathG, ancPathP))


# ===== Ancestry Variance Decomposition =====

for(i in 1:2){
  # repeat for genetic and physical
  
  unit <- d[i, units]
  path <- d[i,path]
  
  # list of individual files 
  indFiles <- as.list(paste0(path,
                             intersect(list.files(path), 
                                       paste0(meta[RI_ACCESSION==pop,ID],".txt"))))
  # combine individuals from population into one table (save these objects independently)
  assign(d[i,gnoms], rbindlist(lapply(indFiles, fread)))
  
  gnoms <- copy(get(d[i,gnoms]))
  
  # ----- mean ancestry -----
  # mean over individuals
  meanAnc <- gnoms[, .(freq_sp2 = mean(freq_sp2)), by=.(position,chr)]
  
  # total mean ancestry 
  totalMeanAnc <- meanAnc[, mean(freq_sp2)]
  
  # ===== Ancestry Variance Decomposition =====
  
  # max level of wavelet decomp
  maxLevel <- meanAnc[, .(level=ceiling(log2(length(position)))), 
                      by = chr][, max(level)]
  
  # ----- mean ancestry -----
  mean_chr_var_decomp <- meanAnc[, haar_dwt_nondyadic_var(.SD, 
                                                          max.level=maxLevel, 
                                                          variable = "freq_sp2"),
                                 by=chr]
  
  # weighted average over chromosomes by scale
  mean_var_decomp <- mean_chr_var_decomp[, .(variance=weighted.mean(variance,w=n.wavelets)),
                                         by=level]
  mean_var_decomp[,signal := "mean"]
  mean_var_decomp[, units := unit]
  
  # ----- individual -----
  ind_chr_var_decomp <- gnoms[, haar_dwt_nondyadic_var(.SD, 
                                                       variable="freq_sp2",
                                                       max.level=maxLevel), 
                              by=.(ID,chr)]
  
  # average over chromosomes within individuals, then average over individuals
  ind_var_decomp <- ind_chr_var_decomp[, .(variance=weighted.mean(variance,w=n.wavelets)),
                                       by=.(ID,level)][, .(variance=mean(variance)), by = level]
  ind_var_decomp[, signal := "individual"]
  ind_var_decomp[, units := unit]
  
  # combine
  if(i==1){
    ancestry_var_decomp <- rbind(ind_var_decomp,mean_var_decomp)
  } else if(i==2){
    ancestry_var_decomp <- rbindlist(list(ancestry_var_decomp,ind_var_decomp,mean_var_decomp),fill=T)
  }
  
  
  # ----- Chromosome-level -----
  
  chr_means <- meanAnc[, .(len = max(position), freq_sp2 = mean(freq_sp2)), by = chr]
  chr_means[, weight := len/sum(len)]
  chr_means[, total_mean := weighted.mean(freq_sp2,weight)]
  
  # add chromosome-level variance to table
  ancestry_var_decomp <- rbindlist(list(ancestry_var_decomp, 
                                        list(level="chrom",
                                             variance=chr_means[, sum(weight*(freq_sp2-total_mean)^2)],
                                             signal="mean",
                                             units=unit)),
                                   fill=T)
  
  # ===== Variance Output =====
  ancestry_var_decomp[units==unit, `:=`(popID = pop,
                                        species = meta[RI_ACCESSION==pop,zea][1],
                                        locality = meta[RI_ACCESSION==pop,LOCALITY][1],
                                        totalMean = totalMeanAnc)][]
}


# ====== Recombination Variance Decomposition =====
# (did physical scale 2nd in loop so maxLevel the same)

rec[, cm_tr := log10(cM)]
rec_var_decomp <- rec[, haar_dwt_nondyadic_var(.SD, variable="cm_tr", max.level = maxLevel), by = chr]

# average over chromosomes
rec_var_decomp <- rec_var_decomp[, .(rec_var = weighted.mean(variance, n.wavelets)),by=level]


# chromosome-level
chr_means_rec <- rec[, .(len = max(position), cm_tr = mean(cm_tr)), by = chr]
chr_means_rec[, weight := len/sum(len)]
chr_means_rec[, total_mean := weighted.mean(cm_tr,weight)]

# add chromosome-level variance to table
rec_var_decomp <- rbindlist(list(rec_var_decomp, 
                                      list(level="chrom",
                                           rec_var=chr_means_rec[, sum(weight*(cm_tr-total_mean)^2)])))

all_var_decomp <- merge(ancestry_var_decomp, rec_var_decomp, all=T)


# ===== Correlation Analysis =====




# ===== Final Output =====
save(all_var_decomp, file = outPath)






