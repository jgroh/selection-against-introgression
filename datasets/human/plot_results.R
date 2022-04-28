library(data.table)
library(waveslim)
library(ggplot2)
library(gridExtra)
library(ggpubr)



wv_physical <- fread("wavelet_results/wv_physical.txt")
wv_genetic <- fread("wavelet_results/wv_genetic.txt")

wc_freq_rec_physical <- fread("wavelet_results/wc_freq_rec_physical.txt")
wc_freq_rec_genetic <- fread("wavelet_results/wc_freq_rec_genetic.txt")

wc_freq_gd_physical <- fread("wavelet_results/wc_freq_gd_physical.txt")
wc_freq_gd_genetic <- fread("wavelet_results/wc_freq_gd_genetic.txt")

wc_freq_gdr_physical <- fread("wavelet_results/wc_freq_gdr_physical.txt")
wc_freq_gdr_genetic <- fread("wavelet_results/wc_freq_gdr_genetic.txt")

wc_rec_gd_physical <- fread("wavelet_results/wc_rec_gd_physical.txt")
wc_rec_gd_genetic <- fread("wavelet_results/wc_rec_gd_genetic.txt")

lapply(
  list(wv_physical, wv_genetic, wc_freq_rec_physical, wc_freq_rec_genetic, wc_freq_gd_physical, wc_freq_gd_genetic, 
            wc_freq_gdr_physical, wc_freq_gdr_genetic, wc_rec_gd_physical, wc_rec_gd_genetic), 
  function(x){
    x[, level := factor(level, levels = c(paste0("d", 1:17), paste0("s", 15:17), "chr"))]
  }
)

# ===== wavelet variance, physical units =====

lineData_physical <- wv_physical[grepl("d", level, fixed = T)]

ggplot(wv_physical, 
       #aes(x = level, y = variance.freq, group = 1)) + 
       aes(x = level, y = propvar.freq, group = 1)) + 
  geom_point(size = 2.2) + 
  geom_line(data = lineData_physical) + 
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



# ===== wavelet variance, genetic units =====

lineData_genetic <- wv_genetic[grepl("d", level, fixed = T)]

ggplot(wv_genetic, 
       #aes(x = level, y = variance.freq, group = 1)) + 
       aes(x = level, y = propvar.freq, group = 1)) + 
  geom_point(size = 2.2) + 
  geom_line(data = lineData_genetic) + 
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



# ===== Correlations =====

# ---- physical scale 

p1 <- ggplot(wc_freq_rec_physical, aes(x = level, y = cor_jack)) + geom_point() + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci)) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = NULL) + #c(as.character(1:17), paste0(15:17, " (s)"), "chromosome"))   + 
  labs(x="",
    #x = expression(Scale: log[2] ("1kb")), 
       y = "freq, rec") + 
  theme(#aspect.ratio = 1, 
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 18, xend = 20, y = -Inf, yend = -Inf))  

p2 <- ggplot(wc_freq_gd_physical, aes(x = level, y = cor_jack)) + geom_point() + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci)) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = NULL) + #c(as.character(1:17), paste0(15:17, " (s)"), "chromosome"))   + 
  labs(x="",
    #x = expression(Scale: log[2] ("1kb")), 
       y = "freq, gene dnsty") + 
  theme(#aspect.ratio = 1, 
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 18, xend = 20, y = -Inf, yend = -Inf))  

p3 <- ggplot(wc_rec_gd_physical, aes(x = level, y = cor_jack)) + geom_point() + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci)) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = c(as.character(1:17), paste0(15:17, " (s)"), "chromosome"))   + 
  labs(x = "",
    #x = expression(Scale: log[2] ("1kb")), 
       y = "rec, gene dnsty") + 
  theme(#aspect.ratio = 1, 
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 18, xend = 20, y = -Inf, yend = -Inf))  

p4 <- ggplot(wc_freq_gdr_physical, aes(x = level, y = cor_jack)) + geom_point() + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci)) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = c(as.character(1:17), paste0(15:17, " (s)"), "chromosome"))   + 
  labs(x="",
       #x = expression(Scale: log[2] ("1kb")), 
       y = "freq, (gene dnsty / rec)") + 
  theme(#aspect.ratio = 1, 
    axis.title = element_text(size = 13),
    axis.text.x = element_text(angle = 90, size = 12),
    axis.text.y = element_text(size = 12),
    axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 18, xend = 20, y = -Inf, yend = -Inf))  


toprow <- ggarrange(NULL, p1, p2, NULL, labels = c("","",""), ncol = 4, widths = c(0.2,1,1,.3))
#toprow
bottomrow <- ggarrange(NULL,p3, p4, NULL, ncol=4, widths = c(0.2,1,1,.3))
#bottomrow
figure <- ggarrange(NULL, toprow, bottomrow, nrow = 3, labels = c("",""), heights = c(0.2,0.8,1))

annotate_figure(figure,
                bottom =  text_grob(expression(Scale: log[2] ("1kb")), hjust=1, size = 15), 
                left = text_grob("Correlation", rot = 90, size = 15)
)


# ---- genetic scale 

p1g <- ggplot(wc_freq_rec_genetic, aes(x = level, y = cor_jack)) + geom_point() + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci)) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = NULL) + #c(as.character(1:17), paste0(15:17, " (s)"), "chromosome"))   + 
  labs(x="",
       y = "freq, rec") + 
  theme(#aspect.ratio = 1, 
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 90, size = 12),
    axis.text.y = element_text(size = 12),
    axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 18, xend = 20, y = -Inf, yend = -Inf))  

p2g <- ggplot(wc_freq_gd_genetic, aes(x = level, y = cor_jack)) + geom_point() + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci)) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = NULL) + #c(as.character(1:17), paste0(15:17, " (s)"), "chromosome"))   + 
  labs(x="",
       y = "freq, gene dnsty") + 
  theme(#aspect.ratio = 1, 
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 90, size = 12),
    axis.text.y = element_text(size = 12),
    axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 18, xend = 20, y = -Inf, yend = -Inf))  

p3g <- ggplot(wc_rec_gd_genetic, aes(x = level, y = cor_jack)) + geom_point() + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci)) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = c(as.character(-16:0), paste0(-2:0, " (s)"), "chromosome"))   + 
  labs(x = "",
       #x = expression(Scale: log[2] ("1kb")), 
       y = "rec, gene dnsty") + 
  theme(#aspect.ratio = 1, 
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 90, size = 12),
    axis.text.y = element_text(size = 12),
    axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 18, xend = 20, y = -Inf, yend = -Inf))  

p4g <- ggplot(wc_freq_gdr_genetic, aes(x = level, y = cor_jack)) + geom_point() + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci)) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = c(as.character(-16:0), paste0(-2:0, " (s)"), "chromosome"))   + 
  labs(x="",
       #x = expression(Scale: log[2] ("1kb")), 
       y = "freq, (gene dnsty / rec)") + 
  theme(#aspect.ratio = 1, 
    axis.title = element_text(size = 13),
    axis.text.x = element_text(angle = 90, size = 12),
    axis.text.y = element_text(size = 12),
    axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf)) + 
  geom_segment(aes(x = 18, xend = 20, y = -Inf, yend = -Inf))  


toprow_g <- ggarrange(NULL, p1g, p2g, NULL, labels = c("","",""), ncol = 4, widths = c(0.2,1,1,.3))
#toprow
bottomrow_g <- ggarrange(NULL,p3g, p4g, NULL, ncol=4, widths = c(0.2,1,1,.3))
#bottomrow
figure_g <- ggarrange(NULL, toprow_g, bottomrow_g, nrow = 3, labels = c("",""), heights = c(0.2,0.8,1))

annotate_figure(figure_g,
                bottom =  text_grob(expression(Scale: log[2] (Morgan)), hjust=1, size = 15), 
                left = text_grob("Correlation", rot = 90, size = 15)
)




# ===== Contribution to Correlation ====

# ---- physical -----
allWav_physical <- merge(wv_physical, wc_freq_gdr_physical)
allWav_physical[, contribution := cor_jack*sqrt(propvar.freq*propvar.rec)]

ggplot(allWav_physical) + 
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