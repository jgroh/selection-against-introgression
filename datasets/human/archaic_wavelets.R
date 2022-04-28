library(data.table)
library(waveslim)
library(ggplot2)

source("~/workspace/gnomwav/R/multi_modwts.R")
source("~/workspace/gnomwav/R/variance_decomp.R")
source("~/workspace/gnomwav/R/correlation_decomp.R")


chr_files1kb <- dir(path = "archaic_freqs/", pattern = "chr.*_frq_1kb_windows.txt", full.names=T)
chr_filesM <- dir(path = "archaic_freqs/", pattern = "chr.*_frq_genetic_windows.txt", full.names=T)

gnom1kb <- rbindlist(lapply(chr_files1kb, fread))
#ggplot(gnom1kb, aes(x = rec, y = freq)) + geom_point()  

gnomM <- rbindlist(lapply(chr_filesM, fread))
#ggplot(gnomM, aes(x = rec, y = freq)) + geom_point()  

#gnom10kb[, cor.test(rec,freq, method = "spearman")]


# ===== wavelet variance, physical units =====

wv1kb <- gnom1kb[, gnom_var_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T, avg.over.chroms = T)]
wv1kb[, level := factor(level, levels = c(paste0("d", 1:17), paste0("s", 15:17), "chr"))]

lineData1kb <- wv1kb[grepl("d", level, fixed = T)]
ggplot(wv1kb, 
       #aes(x = level, y = variance.freq, group = 1)) + 
       aes(x = level, y = propvar.freq, group = 1)) + 
  geom_point(size = 2.2) + 
  geom_line(data = lineData1kb) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = c(as.character(1:17), paste0(15:17, " (scaling var)"), "chromosome"))   + 
  labs(x = expression(Scale: log[2] ("1kb")), 
       #y = "Variance") + 
       y = "Proportion of genome-wide variance") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 18, xend = 20, y = -Inf, yend = -Inf))  



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
wvM[, level := factor(level, levels = c(paste0("d", 1:17), paste0("s", 15:17), "chr"))]

lineDataM <- wvM[grepl("d", level, fixed = T)]
ggplot(wvM, 
       #aes(x = level, y = variance.freq, group = 1)) + 
       aes(x = level, y = propvar.freq, group = 1)) + 
  geom_point(size = 2.2) + 
  geom_line(data = lineDataM) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = c(as.character(-16:0), paste0(-2:0, " (scaling var)"), "chromosome"))   + 
  labs(x = expression(Scale: log[2] (Morgan)), 
       #y = "Variance") + 
       y = "Proportion of genome-wide variance") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 18, xend = 20, y = -Inf, yend = -Inf))  

# all chromosomes separately 
# wvM_chr <- gnomM[, gnom_var_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T, avg.over.chroms = F)]
# wvM_chr[, level := factor(level, levels = c(paste0("d", 1:13), paste0("s", 10:13), "chr"))]
# 
# ggplot(wvM_chr[chr == 1], aes(x = level, y = variance.freq, group = as.factor(chr), color = as.factor(chr))) + geom_point() +
#   geom_line()




# ===== Correlations =====

# ---- physical scale 
wc1kb <- gnom1kb[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T)]
wc1kb[, level := factor(level, levels = c(paste0("d", 1:17), paste0("s", 15:17), "chr"))]

ggplot(wc1kb, aes(x = level, y = cor_jack)) + geom_point() + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci)) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = c(as.character(1:17), paste0(15:17, " (scaling var)"), "chromosome"))   + 
  labs(x = expression(Scale: log[2] ("1kb")), 
       y = "Correlation (archaic allele freq, recomb)") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 18, xend = 20, y = -Inf, yend = -Inf))  


# ---- genetic scale 
wcM <- gnomM[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("freq", "rec"), rm.boundary = T)]
wcM[, level := factor(level, levels = c(paste0("d", 1:17), paste0("s", 15:17), "chr"))]


ggplot(wcM, aes(x = level, y = cor_jack)) + geom_point() + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci)) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = c(as.character(1:17), paste0(15:17, " (scaling var)"), "chromosome"))   +
  labs(x = expression(Scale: log[2] (Morgan)), 
       y = "Correlation (archaic allele freq, recomb)") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 18, xend = 20, y = -Inf, yend = -Inf))  



# ===== Contribution to Correlation ====

# ---- physical -----
allWav1kb <- merge(wv1kb, wc1kb)
allWav1kb[, contribution := cor_jack*sqrt(propvar.freq*propvar.rec)]
allWav1kb[, level := factor(level, levels = c(paste0("d", 1:17), paste0("s", 15:17), "chr"))]

ggplot(allWav1kb) + 
  geom_bar(aes(fill = cor_jack, x = level, y = contribution), stat= "identity", color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = c(as.character(1:17), paste0(15:17, " (scaling)"), "chromosome"))   +
  theme_classic() +
  labs(x = expression(Scale: log[2] ("1kb")), 
       y = "Contribution to total correlation", 
       fill = "Correlation") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12))


# ----- genetic ----
allWavM <- merge(wvM, wcM, by = "level")
allWavM[, contribution := cor_jack*sqrt(propvar.freq*propvar.rec)]
allWavM[, level := factor(level, levels = c(paste0("d", 1:17), paste0("s", 15:17), "chr"))]


ggplot(allWavM) + 
  geom_bar(aes(fill = cor_jack, x = level, y = contribution), stat= "identity", color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = c(as.character(1:17), paste0(15:17, " (scaling)"), "chromosome"))   +
  labs(x = expression(Scale: log[2] (Morgan)), 
       y = "Contribution to total correlation", 
       fill = "Correlation") + 
  theme_classic() +
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12))


# ===== Look at chromosome-level correlation

chrMeans <- gnom10kb[, lapply(.SD, mean), .SDcols = c("rec", "freq"), by = chr]
chrLen <- gnom10kb[, .(nwindow = .N), by = chr]
chrMeans <- merge(chrMeans, chrLen)
ggplot(chrMeans, aes(x = rec, y = freq, size = nwindow)) + geom_point() + geom_text(aes(label = chr), hjust=-.5,vjust=0) +
  theme_classic()

chrMeansM <- gnomM[, lapply(.SD, mean), .SDcols = c("rec", "freq"), by = chr]
chrLenM <- gnomM[, .(nwindow = .N), by = chr]
chrMeansM <- merge(chrMeansM, chrLenM)
ggplot(chrMeansM, aes(x = rec, y = freq, size = nwindow)) + geom_point() + geom_text(aes(label = chr), hjust=-.5,vjust=0) +
  theme_classic()