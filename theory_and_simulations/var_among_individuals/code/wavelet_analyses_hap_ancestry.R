library(data.table)
library(cubature)
library(ggplot2)

# ===== Read and format simulation data =====

if(Sys.getenv("RSTUDIO") == "1"){
  
  source("~/workspace/gnomwav/R/multi_modwts.R")
  source("~/workspace/gnomwav/R/variance_decomp.R")
  source("~/workspace/gnomwav/R/theory.R")
  source("~/workspace/gnomwav/R/correlation_decomp.R")
  
  setwd("/Users/Jeff/workspace/selection-against-introgression/theory_and_simulations/var_among_individuals/results/neutral_sims/equilibrium/")
  n.sample <- 100
  ancestry <- fread("alpha0.5_replicate0_hap_ancestry.txt", col.names = c("gen", "alpha", "rep", "id", "left", "right", "source"))
  
} else {
  
 # source("/Users/brogroh/gnomwav/R/multi_modwts.R")
 # source("/Users/brogroh/gnomwav/R/variance_decomp.R")
 # source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
  #source("/Users/brogroh/gnomwav/R/theory.R")
  source("~/workspace/gnomwav/R/multi_modwts.R")
  source("~/workspace/gnomwav/R/variance_decomp.R")
  source("~/workspace/gnomwav/R/theory.R")
  source("~/workspace/gnomwav/R/correlation_decomp.R")
  

  args <- commandArgs(trailingOnly = TRUE)
  n.sample <- as.numeric(args[1])
  ancestry <- fread(args[2],col.names = c("gen", "alpha", "rep", "id", "left", "right", "source"))
  
}
ancestry[gen != 0, gen := gen-1]
ancestry[, alpha_empirical := mean(source), by = .(rep, gen)]

# ===== Wavelet Variance of True Ancestry =====

xout = seq(2^-16, 1, by = 2^-16)

# ----- Method 1: directly measure ancestry at fine scales, most accurate comparison to theory (but less comparable to snp stat results due to differences in interpolation)

ancestry[, leftMorgan := left/2^26]
ancestry_grid1 <- ancestry[, approx(x = leftMorgan, y = source, xout = xout, method = 'constant', rule = 2), 
                           by = .(id,gen,alpha,alpha_empirical,rep)] 
setnames(ancestry_grid1, c('x','y'), c('Morgan', 'source'))

# wavelet variance per haplotype
hap_true_ancestry_wv <- ancestry_grid1[, gnom_var_decomp(.SD, chromosome = NA, signals = "source"), 
                                       by = .(alpha,alpha_empirical,rep,gen,id)]

# average over haplotypes
hap_true_ancestry_wv <- hap_true_ancestry_wv[, .(single_hap = mean(variance.source)), 
                                             by = .(rep, alpha, alpha_empirical, gen, level)]

# ====== Wavelet Variance of sample mean of true ancestry =====

# reformat so that haplotypes are in separate columns to take mean

#ancestry_grid1[, pos := seq_len(.N), by = .(rep, alpha, alpha_empirical, gen, id)]
#keycols <- c('gen', 'id', 'Morgan')
#setkeyv(ancestry_grid1, cols = keycols)
#ancestry_grid1[, id2 := rep(rep(paste0('id',1:n.sample), each = length(xout)), length(unique(gen)))]
#ancestry_grid1_wide <- dcast(ancestry_grid1, rep + alpha + alpha_empirical + gen + Morgan ~ id2, value.var = 'source')

#smpl_mean_true_ancestry <- ancestry_grid1_wide[, .(smpl_mean = rowSums(.SD)/n.sample), .SDcols = paste0('id',1:n.sample), 
#                                               by = .(rep, alpha, alpha_empirical, gen, Morgan)]

#smpl_mean_true_ancestry_wv <- smpl_mean_true_ancestry[, gnom_var_decomp(.SD, signals = 'smpl_mean', chromosome = NA), 
#                                                      by = .(rep, alpha, alpha_empirical, gen)]
#setnames(smpl_mean_true_ancestry_wv, 'variance.smpl_mean',  'smpl_mean')

# combine results for output
#true_ancestry_allWV <- rbind(melt(hap_true_ancestry_wv, variable.name = 'signal', 
#     id.vars = c('alpha', 'alpha_empirical', 'rep','gen','level'), value.name = 'variance'),
#     melt(smpl_mean_true_ancestry_wv, variable.name = 'signal',
#     id.vars = c('alpha', 'alpha_empirical', 'rep','gen','level'), value.name = 'variance'))


true_ancestry_allWV <- melt(hap_true_ancestry_wv, variable.name = 'signal', 
            id.vars = c('alpha', 'alpha_empirical', 'rep','gen','level'), value.name = 'variance')

# ===== Return individual haplotype's mean ancestry ==== 

ind_gnom_avg <- ancestry_grid1[, .(gnom_avg = mean(source)), by = .(gen, alpha, alpha_empirical, rep, id)] 



save(true_ancestry_allWV, ind_gnom_avg, file = gsub('_hap_ancestry.txt','_wavelet_results.RData', args[2]))
