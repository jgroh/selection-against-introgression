library(data.table)
library(boot)
source("lm_functions.R")

args <- commandArgs(trailingOnly = T)
dwt_file <- args[1]
out_file <- args[2]
n.cores <- args[3]

# 1. ===== Load Data =====

# chromosome info 
#chrs <- fread("match_chromnames_xbir-10x_chroms.txt", header = F)
#chrLen <- fread("xbir10x_chrlengths.txt", col.names = c("chr", "length"))

#acua2018chrMeans <- loadFrom("ACUA_2018/chrSignalMeans.RData", "chrSignalMeans")

# dwt data
acua2018dwt <- loadFrom(dwt_file, "allDWT")

# calculated adjusted r squared estimates
r.squared_dwt <- acua2018dwt[, r.squared_lm(.SD, 
                                                z=ancestry_coeff, 
                                                formula1=ancestry_coeff~rec_coeff,
                                                formula2=ancestry_coeff~rec_coeff+cds_coeff), 
                                      by=scale]

# ===== bootstrap confidence intervals =====

r.squared_dwt_CI <- acua2018dwt[, r.sqrd_boot_all(.SD), by = scale]

dwt_lm_plot_data <- merge(r.squared_dwt,r.squared_dwt_CI, all = T)

dwt_lm_plot_data[, scale := as.numeric(gsub("d","",scale))]

write.table(lm_plot_data, file = out_file, quote=F, sep="\t", row.names=F, col.names=T)

#ggplot(r.squared_plot_data, aes(x = scale, y = r.squared_adj, colour = model)) + 
#  geom_point() +
#  geom_errorbar(aes(ymin=CI_l, ymax=CI_u))

