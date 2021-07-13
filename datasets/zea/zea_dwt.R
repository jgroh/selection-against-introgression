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
ancPath <- args[3]
#recFile <- args[5]
outPath <- args[4]

# run this block only if running locally
# metaFile <- "HILO_MAIZE55_PARV50_meta.txt"
# pop <- "RIMME0035" # mexicana
# pop <- "RIMMA0366" # maize
# ancPath <- "hmm_anc_interp/genetic/"
# recFile <- "zea_1kb_rec.bed"


# Read files
meta <- fread(metaFile)
rec <- fread(recFile, col.names = c('chr','start','end','cM'))

# list of individual files 
indFiles <- as.list(paste0(ancPath,
                           intersect(list.files(ancPath), 
                                     paste0(meta[RI_ACCESSION==pop,ID],".txt"))))
# combine individuals from population into one table
gnoms <- rbindlist(lapply(indFiles, fread))


# merge ancestry data and genomic features data
# rec[, position := start + 500] # this is the midpoint of the 1kb intervals where I interpolate ancestry
# gnoms <- merge(gnoms, rec, by = c("chr", "position"), all=T)


# ===== Data transformations =====

# ----- log transform recombination -----
#gnoms[, cmTr := log10(cM)]
#qqnorm(rec[chr==1,cmTr])


# ----- mean ancestry -----
# mean over individuals
meanAnc <- gnoms[, .(freq_sp2 = mean(freq_sp2)), by=.(position,chr)]

# total mean ancestry 
totalMeanAnc <- meanAnc[, .(freq_sp2 = mean(freq_sp2))]


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

# ----- individual -----
ind_chr_var_decomp <- gnoms[, haar_dwt_nondyadic_var(.SD, 
                                                     variable="freq_sp2",
                                                     max.level=maxLevel), 
                            by=.(ID,chr)]

# average over chromosomes within individuals, then average over individuals
ind_var_decomp <- ind_chr_var_decomp[, .(variance=weighted.mean(variance,w=n.wavelets)),
                   by=.(ID,level)][, .(variance=mean(variance)), by = level]
ind_var_decomp[, signal := "individual"]

# combine
ancestry_var_decomp <- rbind(ind_var_decomp,mean_var_decomp)


# ----- Chromosome-level -----

chr_means <- meanAnc[, .(len = max(position), freq_sp2 = mean(freq_sp2)), by = chr]
chr_means[, weight := len/sum(len)]
chr_means[, total_mean := weighted.mean(freq_sp2,weight)]

# add chromosome-level variance to table
ancestry_var_decomp <- rbindlist(list(ancestry_var_decomp, 
                                      list(level="chrom",
                                           variance=chr_means[, sum(weight*(freq_sp2-total_mean)^2)],
                                           signal="mean")))
     
# ===== Variance Output =====
ancestry_var_decomp[,popID := pop]
ancestry_var_decomp[,species := meta[RI_ACCESSION==pop,zea][1]]
ancestry_var_decomp[,locality := meta[RI_ACCESSION==pop,LOCALITY][1]]
ancestry_var_decomp[,totalMean := totalMeanAnc]


# ===== Final Output =====
save(ancestry_var_decomp, file = outPath)






