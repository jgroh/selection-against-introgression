library(data.table)
library(waveslim)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(cubature)

setwd("~/workspace/selection-against-introgression/datasets/human/")
source("~/workspace/gnomwav/R/theory.R")

wc_calls_listy <- list()
wc_calls_files <- list.files("wavelet_results", pattern = "wc_calls*", full.names=T)
for(i in 1:length(wc_calls_files)){
  f <- wc_calls_files[i]
  d <- strsplit(gsub(".txt", "", gsub("wavelet_results/wc_calls_", "", f)), "_")
  ff <- fread(f)
  ff[, c("units", "assembly", "thresh") := as.list(d[[1]])]
  wc_calls_listy[[i]] <- ff
}
wc_calls <- rbindlist(wc_calls_listy)

wv_listy <- list()
wv_files <- list.files("wavelet_results", pattern = "wv_*", full.names=T)
for(i in 1:length(wv_files)){
  f <- wv_files[i]
  d <- strsplit(gsub(".txt", "", gsub("wavelet_results/wv_", "", f)), "_")
  ff <- fread(f)
  ff[, c("units", "assembly", "thresh") := as.list(d[[1]])]
  wv_listy[[i]] <- ff
}
wv <- rbindlist(wv_listy)

wc_freq_rec_listy <- list()
wc_freq_rec_files <- list.files("wavelet_results", pattern = "wc_freq_rec_*", full.names=T)
for(i in 1:length(wc_freq_rec_files)){
  f <- wc_freq_rec_files[i]
  d <- strsplit(gsub(".txt", "", gsub("wavelet_results/wc_freq_rec_", "", f)), "_")
  ff <- fread(f)
  ff[, c("units", "assembly", "thresh") := as.list(d[[1]])]
  wc_freq_rec_listy[[i]] <- ff
}
wc_freq_rec <- rbindlist(wc_freq_rec_listy)

wc_freq_log10rec_listy <- list()
wc_freq_log10rec_files <- list.files("wavelet_results", pattern = "wc_freq_log10rec_*", full.names=T)
for(i in 1:length(wc_freq_log10rec_files)){
  f <- wc_freq_log10rec_files[i]
  d <- strsplit(gsub(".txt", "", gsub("wavelet_results/wc_freq_log10rec_", "", f)), "_")
  ff <- fread(f)
  ff[, c("units", "assembly", "thresh") := as.list(d[[1]])]
  wc_freq_log10rec_listy[[i]] <- ff
}
wc_freq_log10rec <- rbindlist(wc_freq_log10rec_listy)



wc_freq_cds_listy <- list()
wc_freq_cds_files <- list.files("wavelet_results", pattern = "wc_freq_cds_[^M]*", full.names=T)
for(i in 1:length(wc_freq_cds_files)){
  f <- wc_freq_cds_files[i]
  d <- strsplit(gsub(".txt", "", gsub("wavelet_results/wc_freq_cds_", "", f)), "_")
  ff <- fread(f)
  ff[, c("units", "assembly", "thresh") := as.list(d[[1]])]
  wc_freq_cds_listy[[i]] <- ff
}
wc_freq_cds <- rbindlist(wc_freq_cds_listy)


wc_freq_cdsM_listy <- list()
wc_freq_cdsM_files <- list.files("wavelet_results", pattern = "wc_freq_cdsM", full.names=T)
for(i in 1:length(wc_freq_cdsM_files)){
  f <- wc_freq_cdsM_files[i]
  d <- strsplit(gsub(".txt", "", gsub("wavelet_results/wc_freq_cdsM_", "", f)), "_")
  ff <- fread(f)
  ff[, c("units", "assembly", "thresh") := as.list(d[[1]])]
  wc_freq_cdsM_listy[[i]] <- ff
}
wc_freq_cdsM <- rbindlist(wc_freq_cdsM_listy)


wc_rec_cds_listy <- list()
wc_rec_cds_files <- list.files("wavelet_results", pattern = "wc_rec_cds", full.names=T)
for(i in 1:length(wc_rec_cds_files)){
  f <- wc_rec_cds_files[i]
  d <- strsplit(gsub(".txt", "", gsub("wavelet_results/wc_rec_cds_", "", f)), "_")
  ff <- fread(f)
  ff[, c("units", "assembly", "thresh") := as.list(d[[1]])]
  wc_rec_cds_listy[[i]] <- ff
}
wc_rec_cds <- rbindlist(wc_rec_cds_listy)

lm_listy <- list()
lm_files <- list.files("wavelet_results", pattern = "lm", full.names=T)
for(i in 1:length(lm_files)){
  f <- lm_files[i]
  d <- strsplit(gsub(".txt", "", gsub("wavelet_results/lm_", "", f)), "_")
  ff <- fread(f)
  ff[, c("units", "assembly", "thresh") := as.list(d[[1]])]
  lm_listy[[i]] <- ff
}
lm_results <- rbindlist(lm_listy)

wc_freq_B_listy <- list()
wc_freq_B_files <- list.files("wavelet_results", pattern = "wc_freq_B", full.names=T)
for(i in 1:length(wc_freq_B_files)){
  f <- wc_freq_B_files[i]
  d <- strsplit(gsub(".txt", "", gsub("wavelet_results/wc_freq_B_", "", f)), "_")
  ff <- fread(f)
  ff[, c("units", "assembly", "thresh") := as.list(d[[1]])]
  wc_freq_B_listy[[i]] <- ff
}
wc_freq_B <- rbindlist(wc_freq_B_listy)

lapply(
  list(wc_calls, wv, wc_freq_rec, wc_freq_log10rec, wc_freq_cds, wc_rec_cds, wc_freq_cdsM, wc_freq_B, lm_results), 
  function(x){
    x[, level := factor(level, levels = c(paste0("d", 1:17), paste0("s", 14:17), "chr"))]
  }
)

# ===== Compare calls =====

# ---- physical (assembly makes little difference)
ggplot(wc_calls[!grepl("s", level, fixed = T) & units == 'physical' & assembly == 'hg19' & thresh == 'thresh'], aes(x = level, y = cor_n, color = study)) + 
  #geom_point(size = 2.2, aes(shape = thresh)) + 
  geom_point(size = 2.2) + 
  
  #geom_line(data = wc_calls[grepl("d", level, fixed = T) & units == 'physical' & assembly == 'hg19'], aes(group = interaction(study, thresh), linetype = thresh)) +
  geom_line(data = wc_calls[grepl("d", level, fixed = T) & units == 'physical' & assembly == 'hg19' & thresh == 'thresh'], aes(group = study)) +
  
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), 'scl', "chr"), 
                   labels = c(as.character(0:16),"scl", "chrom"))   + 
  labs(x = expression(Scale: log[2] ("kb")), 
       color = "Comparison", 
       y = "Correlation") + 
  #scale_color_manual(values = c("#19CEBF", "#EFCA2F", "#E87DE0")) +
  scale_color_brewer(type = 'qual', palette = 'Accent') +
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12)) +
  geom_segment(aes(x = 0, xend = 17.05, y = -Inf, yend = -Inf), color = 'black') 

# ---- genetic
ggplot(wc_calls[!grepl("s", level, fixed = T) & units == 'genetic' & assembly == 'hg19'], aes(x = level, y = cor_n, color = study)) + 
  geom_point(size = 2.2, aes(shape = thresh)) + 
  geom_line(data = wc_calls[grepl("d", level, fixed = T) & units == 'genetic' & assembly == 'hg19'], aes(group = interaction(study, thresh), linetype = thresh)) +
  
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), 'scl', 'chr'), 
                   labels = c(as.character(-16:0), 'scl', 'chrom')) + 
  labs(x = expression(Scale: log[2] (Morgan)), 
       color = "Comparison",
       y = "Correlation") + 
  #scale_color_manual(values = c("#19CEBF", "#EFCA2F", "#E87DE0")) +
  scale_color_brewer(type = 'qual', palette = 'Accent') +
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12)) +
  geom_segment(aes(x = 0, xend = 17.05, y = -Inf, yend = -Inf), color = 'black') 



# ===== wavelet variance  =====

# ----- physical units


wv_physical_collapsed <- rbind(wv[grepl('s', level, fixed=T) & units == 'physical', 
                                         .(level = 'scl', variance.skov_freq = sum(variance.skov_freq), 
                                           variance.sank_freq = sum(variance.sank_freq), 
                                           variance.stein_freq = sum(variance.stein_freq),
                                           propvar.skov_freq = sum(propvar.skov_freq),
                                           propvar.sank_freq = sum(propvar.sank_freq),
                                           propvar.stein_freq = sum(propvar.stein_freq)), by= .(assembly, thresh)], 
                              wv[!grepl('s', level, fixed=T) & units == 'physical' & assembly == 'hg19', 
                                         .(level, variance.skov_freq, propvar.skov_freq, 
                                           variance.sank_freq, propvar.sank_freq, 
                                           variance.stein_freq, propvar.stein_freq), by = .(assembly, thresh)]
)
wv_physical_collapsed[, level := factor(level, levels = c(paste0("d", 1:17), 'scl', 'chr'))] 

wv_physical_collapsed <- melt(wv_physical_collapsed, measure.vars = c("variance.skov_freq", "variance.sank_freq", "variance.stein_freq"), 
     variable.name = "study", value.name = "variance")
wv_physical_collapsed[, study := gsub("_freq", "", gsub("variance.","",study))]

ggplot(wv_physical_collapsed[assembly == 'hg19'], 
       aes(x = level, y = variance, color = study)) + 
  theme_classic() + 
  geom_point(size = 2.2, aes(shape = thresh)) +
  geom_line(data=wv_physical_collapsed[grepl("d", level)], aes(group = interaction(thresh, study), linetype = thresh)) +
  scale_x_discrete(breaks = c(paste0("d", 1:17), 'scl', "chr"), 
                   labels = c(as.character(0:16),"scl", "chrom"))   + 
  labs(x = expression(Scale: log[2] ("kb")), 
       #y = "Variance") + 
       y = "Variance") + 
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_blank())  +
  scale_color_brewer(type = 'qual', palette = 'Set2') +
  geom_segment(aes(x = 0, xend = 17.05, y = -Inf, yend = -Inf), color = 'black') 


# ===== wavelet variance, genetic units =====

wvtheory <- wavelet_variance_equilbrium(n.pop=2000, n.sample = 2000, unit.dist = 2^-16, level = 1:17, alpha = 0.01, gen = 2000)
setnames(wvtheory, "variance", "variance.freq")
wvtheory[, propvar.freq := variance.freq/sum(variance.freq)]
wvtheory[, data := "theory"]
wvtheory[, level := paste0("d", level)]
wvtheory[, level := factor(level, levels = c(paste0("d", 1:17), 'scl', 'chr'))] 

# reformat data for plotting
lineData_genetic_hg19 <- wv_genetic_hg19[grepl("d", level, fixed = T)]


wv_genetic_collapsed <- rbind(wv[grepl('s', level, fixed=T) & units == 'genetic', 
                                  .(level = 'scl', 
                                    variance.skov_freq = sum(variance.skov_freq), 
                                    variance.sank_freq = sum(variance.sank_freq), 
                                    variance.stein_freq = sum(variance.stein_freq),
                                    propvar.skov_freq = sum(propvar.skov_freq),
                                    propvar.sank_freq = sum(propvar.sank_freq),
                                    propvar.stein_freq = sum(propvar.stein_freq),
                                    variance.skov_freq.jack.se = NA, 
                                    variance.sank_freq.jack.se = NA, 
                                    variance.stein_freq.jack.se = NA), 
                                 by= .(assembly, thresh)], 
                               wv[!grepl('s', level, fixed=T) & units == 'genetic' & assembly == 'hg19', 
                                  .(level, variance.skov_freq, propvar.skov_freq, 
                                    variance.sank_freq, propvar.sank_freq, 
                                    variance.stein_freq, propvar.stein_freq, variance.skov_freq.jack.se, 
                                    variance.sank_freq.jack.se, 
                                    variance.stein_freq.jack.se), by = .(assembly, thresh)]
)
wv_genetic_collapsed[, level := factor(level, levels = c(paste0("d", 1:17), 'scl', 'chr'))] 

wv_genetic_collapsed_variances <- melt(wv_genetic_collapsed[, .(assembly, thresh, level, variance.skov_freq, 
                                                                variance.sank_freq, variance.stein_freq)], 
                                       measure.vars = c("variance.skov_freq", "variance.sank_freq", "variance.stein_freq"), 
                              variable.name = "study", value.name = "variance")
wv_genetic_collapsed_variances[, study := gsub("_freq", "", gsub("variance.","",study))]

wv_genetic_collapsed_propvars <- melt(wv_genetic_collapsed[, .(assembly, thresh, level, propvar.skov_freq, propvar.sank_freq, propvar.stein_freq)], 
                                      measure.vars = c("propvar.skov_freq", "propvar.sank_freq", "propvar.stein_freq"), 
                                           variable.name = "study", value.name = "propvar")
wv_genetic_collapsed_propvars[, study := gsub("_freq", "", gsub("propvar.", "", study))]

wv_genetic_collapsed_se <- melt(wv_genetic_collapsed[, .(assembly, thresh, level, variance.skov_freq.jack.se, 
                                                                   variance.sank_freq.jack.se, 
                                                                   variance.stein_freq.jack.se)], 
                                     measure.vars = c("variance.skov_freq.jack.se", 
                                                      "variance.sank_freq.jack.se", 
                                                      "variance.stein_freq.jack.se"), 
                                     variable.name = "study", value.name = "se") 
wv_genetic_collapsed_se[, study := gsub("_freq.jack.se", "", gsub("variance.", "", study))]


wv_genetic_plot_data <- merge(merge(wv_genetic_collapsed_variances, wv_genetic_collapsed_propvars), wv_genetic_collapsed_se)
wv_genetic_plot_data[, se.scld := se/sum(variance), by = study]

ggplot(wv_genetic_plot_data[], 
       #aes(x = level, y = variance.freq, group = 1)) + 
       aes(x = level, y = propvar, color = study)) + 
  geom_point(size = 2, aes(shape = thresh)) + 
  scale_color_manual(values = c("#19CEBF", "#EFCA2F", "#E87DE0")) +  
  geom_line(data = wvtheory, aes(x = level, y = propvar.freq), group = 1, color = "black", size = 1) +
  geom_point(size = 2) +
  geom_line(data=wv_genetic_plot_data[grepl('d', level, fixed=T)], linewidth = 1, 
            aes(group = interaction(study, thresh), linetype = thresh)) +
  #geom_errorbar(aes(ymin=propvar - 1.96*se.scld, ymax=propvar+1.96*se.scld)) +
  
  labs(color = "") +
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), 'scl', 'chr'), 
                   labels = c(as.character(-16:0), 'scl', 'chrom')) + 
  labs(x = expression(Scale: log[2] (Morgan)), 
       #y = "Variance") + 
       y = "Proportion of variance") + 
  geom_segment(aes(x = 0, xend = 17.05, y = -Inf, yend = -Inf), color = 'black')  +
  #geom_line(data = wvtheory, aes(x = level, y = propvar.freq), group = 1, color = 'red') +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10))) 


# ===== Correlations =====

# ---- freq vs rec, physical map


p1 <- ggplot(wc_freq_rec[!grepl('s',level,fixed=T) & units == 'physical' & assembly == 'hg19' & study == 'stein'], 
             aes(x = level, y = cor_n)) + 
  geom_abline(slope = 0, intercept = 0, col = 'darkgrey') +
  #geom_point(size=2, color = "#329a9c") + 
  geom_point(size=2, aes(shape = thresh)) + 

  #scale_color_manual(values = c("#19CEBF", "#EFCA2F", "#E87DE0")) +  

  #geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax = cor_jack+1.96*cor_jack_se), width = 0, size = 1) + 
  #geom_point(size=2 ) + 
  
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(0:16, 'chrom')) +
  #scale_y_continuous(limits = c(-.45,0.2)) +
  labs(x = expression(Scale: log[2] ("kb")), 
       y = "Correlation",
       title = "") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf), color = 'black')  +
  theme(aspect.ratio = 1,
      text = element_text(size=15),
      axis.ticks.x = element_line(size=1),
      axis.line.x = element_blank(),
      axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
      axis.text.y = element_text(size=12),
      axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
      plot.title = element_text(hjust = -.1)) 


p1




# ---- freq vs log10rec, physical map

p1.1 <- ggplot(wc_freq_log10rec[!grepl('s',level,fixed=T) & units == 'physical' & assembly == 'hg19' & thresh == 'thresh'], aes(x = level, y = cor_n, color = study)) + 
  geom_abline(slope = 0, intercept = 0, col = 'darkgrey') +
  #geom_point(size=2, color = "#329a9c") + 
  geom_point(size=2, aes(shape = thresh)) + 
  scale_color_manual(values = c("#19CEBF", "#EFCA2F", "#E87DE0")) +  
  
  #geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax = cor_jack+1.96*cor_jack_se), width = 0, size = 1) + 
  geom_point(size=2, aes(shape = thresh)) + 
  
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(0:16, 'chrom')) +
  #scale_y_continuous(limits = c(-.45,0.2)) +
  labs(x = expression(Scale: log[2] ("kb")), 
       y = "Correlation",
       title = "") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf), color = 'black')  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1)) 


p1.1


p2 <- ggplot(wc_freq_cds[!grepl('s',level,fixed=T) & units == 'physical' & assembly == 'hg19' & thresh == 'thresh'], aes(x = level, y = cor_n, color = study)) + 
  geom_abline(slope = 0, intercept = 0, col = 'darkgrey') +
  #geom_point(size=2, color = "#329a9c") + 
  geom_point(size=2, aes(shape = thresh)) + 
  scale_color_manual(values = c("#19CEBF", "#EFCA2F", "#E87DE0")) +  
  
  #geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax = cor_jack+1.96*cor_jack_se), width = 0, size = 1) + 
  geom_point(size=2, aes(shape = thresh)) + 
  
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(0:16, 'chrom')) +
  #scale_y_continuous(limits = c(-.45,0.2)) +
  labs(x = expression(Scale: log[2] ("kb")), 
       y = "Correlation",
       title = "") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf), color = 'black')  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1)) 


p2

p3 <- ggplot(wc_rec_cds[!grepl('s',level,fixed=T) & units == 'physical' & assembly == 'hg19' & thresh == 'thresh' & study == 'skov'], aes(x = level, y = cor_n)) + 
  geom_abline(slope = 0, intercept = 0, col = 'darkgrey') +
  #geom_point(size=2, color = "#329a9c") + 
  geom_point(size=2) + 

  #geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax = cor_jack+1.96*cor_jack_se), width = 0, size = 1) + 
  geom_point(size=2) + 
  
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(0:16, 'chrom')) +
  #scale_y_continuous(limits = c(-.45,0.2)) +
  labs(x = expression(Scale: log[2] ("kb")), 
       y = "Correlation",
       title = "") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf), color = 'black')  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1)) 


p3


p4 <- ggplot(wc_freq_B[!grepl('s',level,fixed=T) & units == 'physical' & assembly == 'hg19' & thresh == 'thresh' ], aes(x = level, y = cor_n)) + 
  geom_abline(slope = 0, intercept = 0, col = 'darkgrey') +
  #geom_point(size=2, color = "#329a9c") + 
  geom_point(size=2, aes(color = study)) + 
  
  #geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax = cor_jack+1.96*cor_jack_se), width = 0, size = 1) + 

  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(0:16, 'chrom')) +
  #scale_y_continuous(limits = c(-.45,0.2)) +
  labs(x = expression(Scale: log[2] ("kb")), 
       y = "Correlation",
       title = "") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf), color = 'black')  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1)) 


p4

library(ggallin)

p5 <- ggplot(wc_freq_cdsM[!grepl('s',level,fixed=T) & units == 'physical' & assembly == 'hg19' & thresh == 'thresh'], aes(x = level, y = cor_n)) + 
  geom_abline(slope = 0, intercept = 0, col = 'darkgrey') +
  #geom_point(size=2, color = "#329a9c") + 
  geom_point(size=2, aes(color = study)) + 
  
  #geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax = cor_jack+1.96*cor_jack_se), width = 0, size = 1) + 
  
  theme_classic() + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(0:16, 'chrom')) +
  scale_y_continuous(trans= pseudolog10_trans) +
  #scale_y_continuous(limits = c(-.45,0.2)) +
  labs(x = expression(Scale: log[2] ("kb")), 
       y = "Correlation",
       title = "") + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf), color = 'black')  +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1)) 


p5





# correlation of recomb and gene density in same plot
wc_freq_rec_physical_hg19[, vars := 'frq_rec']
wc_rec_cds_physical_hg19[, vars := 'rec_cds']
wc_freq_cds_physical_hg19[, vars := 'frq_cds']

ggplot(rbindlist(list(wc_freq_rec_physical_hg19, 
             wc_rec_cds_physical_hg19,
             wc_freq_cds_physical_hg19))[!grepl('s',level)], 
       aes(x = level, y = cor_n, color = vars)) + geom_point() + facet_wrap(~study) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle=90)) + 
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(1:17, 'chrom')) #+
  #geom_errorbar(aes(ymin=cor_n -1.96*cor_jack_se, 
        #        ymax=cor_n + 1.96*cor_jack_se))
 

p3 <- ggplot(wc_rec_cds_physical_hg19[!grepl('s',level,fixed=T)], aes(x = level, y = cor_jack)) + 
  geom_abline(slope = 0, intercept = 0, col = 'darkgrey') +
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax = cor_jack+1.96*cor_jack_se), width = 0, size = 1) + 
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
  geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax = cor_jack+1.96*cor_jack_se), width = 0, size = 1) + 
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


p5 <- ggplot(wc_freq_B_physical_hg19[!grepl('s',level,fixed=T)], aes(x = level, y = cor_jack)) + 
  geom_abline(slope = 0, intercept = 0, col = 'darkgrey') + facet_wrap(~study) +
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax = cor_jack+1.96*cor_jack_se), width = 0, size = 1) + 
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


ggplot(lm_physical[model=='log10rec_cds_density' & !grepl('s',level,fixed=T)], aes(x = level, y = rsqrd_n)) + geom_point() +
  geom_errorbar(aes(ymin=rsqrd_n - 1.96*rsqrd_jack_se, ymax=rsqrd_n + 1.96*rsqrd_jack_se)) + 
  geom_segment(aes(x = 0, xend = 17, y = -Inf, yend = -Inf))  +
  theme_classic() +
  scale_x_discrete(breaks = c(paste0("d", 1:17), "chr"), 
                   labels = c(0:16, 'chrom')) +
  labs(x = expression(Scale: log[2] ("kb")), 
       y = "R squared", 
       title = "") + 
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        plot.title = element_text(hjust = -.1),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10))) 





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





