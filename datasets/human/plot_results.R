library(data.table)
library(waveslim)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(cubature)

setwd("~/workspace/selection-against-introgression/datasets/human/")
source("~/workspace/gnomwav/R/theory.R")


wv_physical <- fread("wavelet_results/wv_physical.txt")
wv_genetic <- fread("wavelet_results/wv_genetic.txt")

wc_freq_rec_physical <- fread("wavelet_results/wc_freq_rec_physical.txt")
wc_freq_rec_genetic <- fread("wavelet_results/wc_freq_rec_genetic.txt")

wc_freq_gd_physical <- fread("wavelet_results/wc_freq_gd_physical.txt")
wc_freq_gd_genetic <- fread("wavelet_results/wc_freq_gd_genetic.txt")

wc_freq_gdr_physical <- fread("wavelet_results/wc_freq_gdr_physical.txt")

wc_rec_gd_physical <- fread("wavelet_results/wc_rec_gd_physical.txt")
wc_rec_gd_genetic <- fread("wavelet_results/wc_rec_gd_genetic.txt")

lm_physical <- fread("wavelet_results/lm_physical.txt")
lm_genetic <- fread("wavelet_results/lm_genetic.txt")

wc_freq_B_physical <- fread("wavelet_results/wc_freq_B_physical.txt")
wc_freq_B_genetic <- fread("wavelet_results/wc_freq_B_genetic.txt")

lapply(
  list(wv_physical, wv_genetic, wc_freq_rec_physical,wc_freq_B_physical, wc_freq_B_genetic, wc_freq_rec_genetic, wc_freq_gd_physical, wc_freq_gd_genetic, 
            wc_freq_gdr_physical,  wc_rec_gd_physical, wc_rec_gd_genetic, lm_physical, lm_genetic), 
  function(x){
    x[, level := factor(level, levels = c(paste0("d", 1:17), paste0("s", 14:17), "chr"))]
  }
)

# ===== wavelet variance, physical units =====

lineData_physical <- wv_physical[grepl("d", level, fixed = T)]
wv_physical_collapsed <- rbind(wv_physical[grepl('s', level, fixed=T), 
                                         .(level = 'scl', variance.freq = sum(variance.freq), 
                                           propvar.freq = sum(propvar.freq))], 
                              wv_physical[!grepl('s', level, fixed=T), 
                                         .(level, variance.freq, propvar.freq)]
)
wv_physical_collapsed[, level := factor(level, levels = c(paste0("d", 1:17), 'scl', 'chr'))] 

ggplot(wv_physical_collapsed, 
       #aes(x = level, y = variance.freq, group = 1)) + 
       aes(x = level, y = variance.freq, group = 1)) + 
  geom_point(size = 2.2) + 
  geom_line(data = lineData_physical) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), 'scl', "chr"), 
                   labels = c(as.character(0:16),"scl", "chrom"))   + 
  labs(x = expression(Scale: log[2] ("1kb")), 
       #y = "Variance") + 
       y = "Variance") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank())  +
  geom_segment(aes(x = 0, xend = 17.05, y = -Inf, yend = -Inf)) 


# ===== wavelet variance, genetic units =====

wvtheory <- wavelet_variance_equilbrium(n.pop=20000, n.sample = 20000, unit.dist = 2^-16, level = 1:17, alpha = 0.005, gen = 2000)
setnames(wvtheory, "variance", "variance.freq")
wvtheory[, propvar.freq := variance.freq/sum(variance.freq)]
wvtheory[, data := "theory"]
wvtheory[, level := paste0("d", level)]
wvtheory[, level := factor(level, levels = c(paste0("d", 1:17), 'scl', 'chr'))] 

lineData_genetic <- wv_genetic[grepl("d", level, fixed = T)]

wv_genetic_collapsed <- rbind(wv_genetic[grepl('s', level, fixed=T), 
                                         .(level = 'scl', variance.freq = sum(variance.freq), 
                                           propvar.freq = sum(propvar.freq))], 
                              wv_genetic[!grepl('s', level, fixed=T), 
                                         .(level, variance.freq, propvar.freq)]
      )

wv_genetic_collapsed[, level := factor(level, levels = c(paste0("d", 1:17), 'scl', 'chr'))] 


ggplot(wv_genetic_collapsed, 
       #aes(x = level, y = variance.freq, group = 1)) + 
       aes(x = level, y = propvar.freq)) + 
  #geom_point(data = wvtheory, aes(x = level, y = propvar.freq), size=3, color = "#c24633") +
  geom_point(size = 3) + 
  
  geom_line(data = wvtheory, aes(x = level, y = propvar.freq), group = 1, color = "darkgrey", size = 1) +
  #geom_line(data = lineData_genetic, group = 1, size = 2, color = 'dodgerblue') + 
  geom_point(size = 3) + 
  
  labs(color = "") +
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), 'scl', 'chr'), 
                   labels = c(as.character(-16:0), 'scl', 'chrom')) + 
  labs(x = expression(Scale: log[2] (Morgan)), 
       #y = "Variance") + 
       y = "Proportion of variance") + 
  geom_segment(aes(x = 0, xend = 17.05, y = -Inf, yend = -Inf))  +
  #geom_line(data = wvtheory, aes(x = level, y = propvar.freq), group = 1, color = 'red') +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10))) 


# ===== Correlations =====

# ---- physical scale 


p1 <- ggplot(wc_freq_rec_physical[!grepl('s',level,fixed=T)], aes(x = level, y = cor_jack)) + 
  geom_abline(slope = 0, intercept = 0, col = 'darkgrey') +
  geom_point(size=2, color = "#eb5ab4") + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci), width = 0, size = 1, color = "#eb5ab4") + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(0:16, 'chrom')) +
  #scale_y_continuous(limits = c(-.45,0.2)) +
  labs(x = expression(Scale: log[2] ("kb")), 
       y = "Correlation",
       title = "") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf))  +
  theme(aspect.ratio = 1,
      text = element_text(size=15),
      axis.ticks.x = element_line(size=1),
      axis.line.x = element_blank(),
      axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
      axis.text.y = element_text(size=12),
      axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
      plot.title = element_text(hjust = -.1)) 


p1


p2 <- ggplot(wc_freq_gd_physical[!grepl('s',level,fixed=T)], aes(x = level, y = cor_jack)) + 
  geom_abline(slope = 0, intercept = 0, col = 'darkgrey') +
  geom_point(size=2, col = '#eda915') + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci), width = 0, size = 1, col = '#eda915') + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(1:17, 'chrom')) +
  labs(x = expression(Scale: log[2] ("kb")), 
       y = "Correlation", 
       title = "B") + 
  scale_y_continuous(limits = c(-.45,0.2)) +
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf))  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1))

p2

# correlation of recomb and gene density in same plot
p1 + geom_point(data = wc_freq_gd_physical[!grepl('s',level,fixed=T)], size=2, col = '#3dc2b9') + 
  geom_errorbar(data = wc_freq_gd_physical[!grepl('s',level,fixed=T)], aes(ymin=lower95ci, ymax = upper95ci), width = 0, size = 1, col = '#3dc2b9') + 
  geom_point(data = wc_rec_gd_physical[!grepl('s',level,fixed=T)], size=2, col = '#e89a1b') + 
  geom_errorbar(data = wc_rec_gd_physical[!grepl('s',level,fixed=T)], aes(ymin=lower95ci, ymax = upper95ci), width = 0, size = 1, col = '#e89a1b')  
  

p3 <- ggplot(wc_rec_gd_physical[!grepl('s',level,fixed=T)], aes(x = level, y = cor_jack)) + 
  geom_abline(slope = 0, intercept = 0, col = 'darkgrey') +
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci), width = 0, size = 1) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(1:17, 'chrom')) +
  labs(x = expression(Scale: log[2] ("kb")), 
       y = "Correlation", 
       title = "A") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf))  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),plot.title = element_text(hjust = -.1)) 

p3
p1 + p2 


p4 <- ggplot(wc_freq_gdr_physical[!grepl('s',level,fixed=T)], aes(x = level, y = cor_jack)) + 
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci), width = 0, size = 1) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(1:17, 'chrom')) +
  labs(x = expression(Scale: log[2] ("kb")), 
       y = "Correlation") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf))  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10))) 

p4


p5 <- ggplot(wc_freq_B_physical[!grepl('s',level,fixed=T)], aes(x = level, y = cor_jack)) + 
  geom_abline(slope = 0, intercept = 0, col = 'darkgrey') +
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci), width = 0, size = 1) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(0:16, 'chrom')) +
  labs(x = expression(Scale: log[2] ("kb")), 
       y = "Correlation", 
       title = "") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf))  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust = -.1),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10))) 

p5
p3 + p5


# ----- genetic units ----
p1.1 <- ggplot(wc_freq_rec_genetic[!grepl('s',level,fixed=T)], aes(x = level, y = cor_jack)) + 
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci), width = 0, size = 1) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), labels = c(-16:0, 'chrom')) +
  labs(x = expression(Scale: log[2] (Morgan)), 
       y = "Correlation") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf))  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10))) 

p1.1


p2.1 <- ggplot(wc_freq_gd_genetic[!grepl('s',level,fixed=T)], aes(x = level, y = cor_jack)) + 
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci), width = 0, size = 1) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), labels = c(-16:0, 'chrom')) +
  labs(x = expression(Scale: log[2] (Morgan)), y = "Correlation") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf))  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10))) 

p2.1


p3.1 <- ggplot(wc_rec_gd_genetic[!grepl('s',level,fixed=T)], aes(x = level, y = cor_jack)) + 
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci), width = 0, size = 1) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), labels = c(-16:0, 'chrom')) +
  labs(x = expression(Scale: log[2] (Morgan)), y = "Correlation") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf))  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10))) 

p1.1 + p2.1 + p3.1


p4.1 <- ggplot(wc_freq_gdr_genetic[!grepl('s',level,fixed=T)], aes(x = level, y = cor_jack)) + 
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci), width = 0, size = 1) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), labels = c(-16:0, 'chrom')) +
  labs(x = expression(Scale: log[2] (Morgan)), y = "Correlation") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf))  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10))) 

p5.1 <- ggplot(wc_freq_B_genetic[!grepl('s',level,fixed=T)], aes(x = level, y = cor_jack)) + 
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=lower95ci, ymax = upper95ci), width = 0, size = 1) + 
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), labels = c(-16:0, 'chrom')) +
  labs(x = expression(Scale: log[2] (Morgan)), y = "Correlation") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf))  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10))) 

p5.1






# toprow <- ggarrange(NULL, p1, p2, NULL, labels = c("","",""), ncol = 4, widths = c(0.2,1,1,.3))
# #toprow
# bottomrow <- ggarrange(NULL,p3, p4, NULL, ncol=4, widths = c(0.2,1,1,.3))
# #bottomrow
# figure <- ggarrange(NULL, toprow, bottomrow, nrow = 3, labels = c("",""), heights = c(0.2,0.8,1))
# 
# annotate_figure(figure,
#                 bottom =  text_grob(expression(Scale: log[2] ("1kb")), hjust=1, size = 15), 
#                 left = text_grob("Correlation", rot = 90, size = 15)
# )


# toprow <- ggarrange(NULL, p1g, NULL,  labels = c("","",""), ncol = 3, widths = c(.3,1,.3))
# middle <- ggarrange(NULL, p2g, NULL,  labels = c("","",""), ncol = 3, widths = c(.3,1,.3))
# bottomrow <- ggarrange(NULL, p3g, NULL,  labels = c("","",""), ncol = 3, widths = c(.3,1,.3))
# 
# 
# #toprow
# #bottomrow
# figure_g <- ggarrange(toprow, middle, bottomrow,  nrow = 3, labels = c("",""), heights = c(0.8,0.8,1))
# 
# annotate_figure(figure_g,
#                 bottom =  text_grob(expression(Scale: log[2] (Morgan)), hjust=1, size = 15), 
#                 left = text_grob("Correlation", rot = 90, size = 15)
# )




# ===== Contribution to Correlation ====

# ---- physical -----
allWav_physical <- merge(wv_physical, wc_freq_rec_physical)
allWav_physical[, contribution := cor_jack*sqrt(propvar.freq*propvar.rec)]

ggplot(allWav_physical[!grepl("s", level,fixed=T)]) + 
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
allWav_genetic <- merge(wv_genetic, wc_freq_gd_genetic, by = "level")
allWav_genetic[, contribution := cor_jack*sqrt(propvar.freq*propvar.gd)]


ggplot(allWav_genetic) + 
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


# ===== lm results =====

ggplot(lm_physical[variable == "rsqrd"], aes(x = level, y = jn_bc_estimate, color = variable)) + 
  geom_abline(slope=0, intercept = 0) + 
  geom_point() + 
  geom_errorbar(aes(ymin=jn_bc_estimate - 1.96*jn_se, ymax=jn_bc_estimate + 1.96*jn_se) ) + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = c(as.character(1:17), paste0(15:17, " (s)"), "chromosome"))   +
  theme_classic() +
  labs(x = expression(Scale: log[2] ("1kb")), 
       y = "Slope estimate") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12)) 



ggplot(lm_physical[variable == "rsqrd" & !grepl("s",level,fixed=T)], aes(x = level, y = jn_bc_estimate)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=jn_bc_estimate - 1.96*jn_se, ymax=jn_bc_estimate + 1.96*jn_se) ) + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), paste0("s", 15:17), "chr"), 
                   labels = c(as.character(1:17), paste0(15:17, " (s)"), "chromosome"))   +
  theme_classic() +
  labs(x = expression(Scale: log[2] ("1kb")), 
       y = "R squared") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12))





