library(data.table)
library(tidyverse)
library(boot)

args <- commandArgs(trailingOnly = T)
dwt_file <- args[1]
modwt_file <- args[2]
n.cores <- args[3]

# 1. ===== Load Data =====

loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 

# chromosome info 
#chrs <- fread("match_chromnames_xbir-10x_chroms.txt", header = F)
#chrLen <- fread("xbir10x_chrlengths.txt", col.names = c("chr", "length"))

#acua2018chrMeans <- loadFrom("ACUA_2018/chrSignalMeans.RData", "chrSignalMeans")

# dwt data
acua2018dwt <- loadFrom(dwt_file, "allDWT")

# modwt data
acua2018modwt <- loadFrom(modwt_file, "allVarsModwtP")


# 2. ----- r squared Functions ------

# From linear model
r.squared_lm <- function(data, z, formula1, formula2){
  # assumes data has coefficients for only a single scale
  fitData <- data[!is.na(z)]
  
  # fit models with intercept and either one or two predictors
  fit1 <- lm(formula1, fitData)
  r.squared_1 <- summary(fit1)$adj.r.squared
  anova1 <- anova(fit1)
  p1 <- anova1$`Pr(>F)`[1]
  
  fit2 <- lm(formula2, fitData)
  r.squared_2 <- summary(fit2)$adj.r.squared
  anova2 <- anova(fit2)
  p2 <- anova2$`Pr(>F)`[1:2]
  
  return(list(model = c("x","xy"),
              r.squared_adj = c(r.squared_1, r.squared_2),
              p.val_rec = c(p1, p2[1]),
              p.val_cds = c(NA, p2[2])))
}


# calculated adjusted r squared estimates
r.squared_dwt <- acua2018dwt[, r.squared_lm(.SD, 
                                                z=ancestry_coeff, 
                                                formula1=ancestry_coeff~rec_coeff,
                                                formula2=ancestry_coeff~rec_coeff+cds_coeff), 
                                      by=scale]


# ===== bootstrap confidence intervals =====

r.sqrd_boot_single_level <- function(data, indices, formula1, formula2){
  d <- data[indices,] # allows boot to select sample
  fit1 <- lm(formula1, data=d)
  fit2 <- lm(formula2, data=d)
  return(c(summary(fit1)$adj.r.squared, summary(fit2)$adj.r.squared))
}

r.sqrd_boot_all <- function(data){
  # runs bootstrap sampling on single level and outputs
  # to be applied to data table grouped by level for output per level 
  
  d <- data[!is.na(ancestry_coeff)]
  
  result <- boot(d, statistic=r.sqrd_boot_single_level, 
                 formula1=ancestry_coeff~rec_coeff,
                 formula2=ancestry_coeff~rec_coeff+cds_coeff,
                 R=10000,
                 parallel="multicore",
                 ncpus=n.cores)
  
  x_CI <- boot.ci(result, type="bca", index=1) # recomb
  xy_CI <- boot.ci(result, type="bca", index=2) # recomb + gene density
  
  return(data.table(model = c("x","xy"),
                    CI_l = c(x_CI$bca[4],xy_CI$bca[4]),
                    CI_u = c(x_CI$bca[5],xy_CI$bca[5]))
  )
  
}

r.squared_dwt_CI <- acua2018dwt[, r.sqrd_boot_all(.SD), by = scale]

dwt_lm_plot_data <- merge(r.squared_dwt,r.squared_dwt_CI, all = T)

dwt_lm_plot_data[, scale := as.numeric(gsub("d","",scale))]

write.table(lm_plot_data, file = "dwt_lm_plot_data.txt", quote=F, sep="\t", row.names=F, col.names=T)

#ggplot(r.squared_plot_data, aes(x = scale, y = r.squared_adj, colour = model)) + 
#  geom_point() +
#  geom_errorbar(aes(ymin=CI_l, ymax=CI_u))

