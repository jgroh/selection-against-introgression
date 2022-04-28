library(data.table)
library(waveslim)
library(ggplot2)

source("~/workspace/gnomwav/R/multi_modwts.R")
source("~/workspace/gnomwav/R/variance_decomp.R")
source("~/workspace/gnomwav/R/correlation_decomp.R")


chr_files10kb <- dir(path = ".", pattern = "chr.*_frq_10kb_windows.txt")
chr_filesM <- dir(path = ".", pattern = "chr.*_frq_genetic_windows.txt")

gnom10kb <- rbindlist(lapply(chr_files10kb, fread))
#ggplot(gnom10kb, aes(x = rec, y = freq)) + geom_point()  

gnomM <- rbindlist(lapply(chr_filesM, fread))
#ggplot(gnomM, aes(x = rec, y = freq)) + geom_point()  

#gnom10kb[, cor.test(rec,freq, method = "spearman")]


# ===== wavelet variance, physical units =====

wv10kb <- gnom10kb[, gnom_var_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T, avg.over.chroms = T)]
wv10kb[, level := factor(level, levels = c(paste0("d", 1:14), paste0("s", 11:14), "chr"))]

lineData10kb <- wv10kb[grepl("d", level, fixed = T)]
ggplot(wv10kb, 
       #aes(x = level, y = variance.freq, group = 1)) + 
      aes(x = level, y = propvar.freq, group = 1)) + 
  geom_point(size = 2.2) + 
  geom_line(data = lineData10kb) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:14), paste0("s", 11:14), "chr"), 
                   labels = c(as.character(1:14), paste0(11:14, " (scaling var)"), "chromosome"))   + 
  labs(x = expression(Scale: log[2] ("10kb")), 
       #y = "Variance") + 
       y = "Proportion of genome-wide variance") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
      axis.text.x = element_text(angle = 90, size = 12),
      axis.text.y = element_text(size = 12),
      axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 14, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 15, xend = 18, y = -Inf, yend = -Inf))  



# wv10kb_chrs <- gnom10kb[, gnom_var_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T, avg.over.chroms = F)]
# wv10kb_chrs[, level := factor(level, levels = c(paste0("d", 1:14), paste0("s", 11:14), "chr"))]
# 
# totalVar10kb <- gnom10kb[, .(total = var(freq)), by = chr]
# totalWV10kb_s <- wv10kb_chrs[, .(sum = sum(variance.freq)),by = chr]
# totalWV10kb_nos <- wv10kb_chrs[grepl("d", level, fixed=T), .(sum_nos = sum(variance.freq)),by = chr]
# 
# totalVars10kb <- merge(totalWV10kb_nos, merge(totalVar10kb, totalWV10kb_s, by = "chr"), by = "chr")
# 
# with(totalVars10kb, plot(sum ~ total, xlab = "Exact total variance", ylab = "Sum of wavelet variances and scaling variance"))
# abline(0,1)
# with(totalVars10kb, plot(sum_nos ~ total, xlab = "Exact total variance", ylab = "Sum of only wavelet variances"))
# abline(0,1)

# ===== wavelet variance, genetic units =====

wvM <- gnomM[, gnom_var_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T, avg.over.chroms = T)]
wvM[, level := factor(level, levels = c(paste0("d", 1:13), paste0("s", 10:13), "chr"))]

lineDataM <- wvM[grepl("d", level, fixed = T)]
ggplot(wvM, 
       #aes(x = level, y = variance.freq, group = 1)) + 
       aes(x = level, y = propvar.freq, group = 1)) + 
  geom_point(size = 2.2) + 
  geom_line(data = lineDataM) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:13), paste0("s", 10:13), "chr"), 
                   labels = c(as.character(-12:0), paste0(-2:1, " (scaling var)"), "chromosome"))   + 
  labs(x = expression(Scale: log[2] (Morgan)), 
       #y = "Variance") + 
       y = "Proportion of genome-wide variance") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 13, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 14, xend = 17, y = -Inf, yend = -Inf))  

# all chromosomes separately 
# wvM_chr <- gnomM[, gnom_var_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T, avg.over.chroms = F)]
# wvM_chr[, level := factor(level, levels = c(paste0("d", 1:13), paste0("s", 10:13), "chr"))]
# 
# ggplot(wvM_chr[chr == 1], aes(x = level, y = variance.freq, group = as.factor(chr), color = as.factor(chr))) + geom_point() +
#   geom_line()




# ===== Correlations =====

# ---- physical scale 
wc10kb <- gnom10kb[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T)]
wc10kb[, level := factor(level, levels = c(paste0("d", 1:14), paste0("s", 11:14), "chr"))]

ggplot(wc10kb, aes(x = level, y = cor_jack)) + geom_point() + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci)) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:14), paste0("s", 11:14), "chr"), 
                   labels = c(as.character(1:14), paste0(11:14, " (scaling var)"), "chromosome"))   + 
  labs(x = expression(Scale: log[2] ("10kb")), 
       y = "Correlation (archaic allele freq, recomb)") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 14, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 15, xend = 18, y = -Inf, yend = -Inf))  


# ---- genetic scale 
wcM <- gnomM[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T)]
wcM[, level := factor(level, levels = c(paste0("d", 1:13), paste0("s", 10:13), "chr"))]


ggplot(wcM, aes(x = level, y = cor_n)) + geom_point() + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci)) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:13), paste0("s", 10:13), "chr"), 
                   labels = c(as.character(-12:0), paste0(-2:1, " (scaling var)"), "chromosome"))   + 
  labs(x = expression(Scale: log[2] (Morgan)), 
       y = "Variance") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 13, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 14, xend = 17, y = -Inf, yend = -Inf))  



# ===== Contribution to Correlation ====

# ---- physical -----
allWav10kb <- merge(wv10kb, wc10kb)
allWav10kb[, contribution := cor_jack*sqrt(propvar.freq*propvar.rec)]
allWav10kb[, level := factor(level, levels = c(paste0("d", 1:14), paste0("s", 11:14), "chr"))]

ggplot(allWav10kb) + 
  geom_bar(aes(fill = cor_jack, x = level, y = contribution), stat= "identity", color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  scale_x_discrete(breaks = c(paste0("d", 1:14), paste0("s", 11:14), "chr"), 
                   labels = c(as.character(1:14), paste0(11:14, " (scaling)"), "chromosome")) +
 theme_classic() +
  labs(x = expression(Scale: log[2] ("10kb")), 
       y = "Contribution to total correlation", 
       fill = "Correlation") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12))


# ----- genetic ----
allWavM <- merge(wvM, wcM, by = "level")
allWavM[, contribution := cor_jack*sqrt(propvar.freq*propvar.rec)]
allWavM[, level := factor(level, levels = c(paste0("d", 1:13), paste0("s", 10:13), "chr"))]


ggplot(allWavM) + 
  geom_bar(aes(fill = cor_jack, x = level, y = contribution), stat= "identity", color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  scale_x_discrete(breaks = c(paste0("d", 1:13), paste0("s", 10:13), "chr"), 
                   labels = c(as.character(-12:0), paste0(-2:1, " (scaling)"), "chromosome"))   + 
  labs(x = expression(Scale: log[2] (Morgan)), 
       y = "Contribution to total correlation", 
       fill = "Correlation") + 
  theme_classic() +
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12))



