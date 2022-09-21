library(data.table)
library(ggplot2)
source("~/workspace/gnomwav/R/multi_modwts.R")
source("~/workspace/gnomwav/R/variance_decomp.R")

#maize_frqs <- fread("allopatric_maize.mafs.gz")
#mex_frqs <- fread("allopatric_mexicana.mafs.gz")
#hyb_frqs <- fread("pop18.mafs.gz")
#hyb_frqs <- fread("HILO9.mafs.gz")
#rmap <- fread("ogut_2015_rmap_v2_to_v4_EXTENDED.txt")


args <- commandArgs(trailingOnly = T)
maize_frqs <- fread(args[1]) # file path
mex_frqs <- fread(args[2])
hyb_frqs <- fread(args[3])
rmap <- fread(args[4])


d <- merge(merge(maize_frqs[,.(chromo, position, maize_frq = phat)],
            mex_frqs[,.(chromo, position, mex_frq = phat)]),
      hyb_frqs[,.(chromo, position, hyb_frq = phat)])


# calculate admixture stat (estimator of mexicana ancestry)
d[, snp_stat := (hyb_frq-maize_frq)/(mex_frq-maize_frq)]

# output global ancestry and fraction of sites with nonzero coverage
mex_ancestry <- d[, mean(snp_stat)]
nonzero_covrg <- nrow(hyb_frqs)/nrow(maize_frqs)
dt <- data.table(group = args[3], mex_ancestry = mex_ancestry, nonzero_covrg = nonzero_covrg)
fwrite(dt, file = "avgAncestry_and_nonzeroCovrg.txt", append=T, quote=F, sep="\t", col.names=!file.exists("avgAncestry_and_nonzeroCovrg.txt"))

# get genetic position of informative SNPs
snp_cm <- rmap[, approx(x=pos_bp, y=pos_cM, xout=d[chromo == .BY, position]), by = chr]
setnames(snp_cm, c("x","y"),c("position","cm"))

d <- merge(d, snp_cm[, .(chromo=chr,position,cm)])
d[, Morgan := cm/100]

# interpolate snp stat to evenly spaced genetic intervals
# resolution: 2^-14 Morgans
# d[, cm-shift(cm), by = chromo][, mean(V1, na.rm=T)/100]

interp <- d[, approx(x=Morgan, y=snp_stat, xout=seq(min(Morgan), max(Morgan), by = 2^-14) ), by = chromo]
setnames(interp, c("x", "y"), c("Morgan", "snp_stat"))

wv <- interp[, gnom_var_decomp(.SD, chromosome = 'chromo', signals = 'snp_stat')]
wv[, level := factor(level, levels = c(paste0('d', 1:15), 's14', 's15', 'chromo'))]

#ggplot(wv, aes(x = level, y = variance.snp_stat)) + geom_point() + 
#  scale_x_discrete(breaks = c(paste0("d",1:15),  's14', 's15', 'chromo'), labels = c(as.character(-14:0), '-1(s)', '0(s)', 'chrom')) +
#  labs(x = expression(Scale: log[2] (Morgan)), y = 'Variance')  + 
#  theme(axis.text.x = element_text(angle=90))

wv[, signal := gsub(".mafs.gz", "", args[3])]
fwrite(wv, file = "snp_stat_wavelet_variances.txt", append=T, quote =F, sep = '\t', row.names=F, col.names=!file.exists("snp_stat_wavelet_variances.txt"))
