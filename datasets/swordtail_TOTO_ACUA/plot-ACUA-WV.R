library(data.table)
library(tidyverse)
library(gridExtra)
library(grid)
library(cubature)
source("~/workspace/gnomwav/R/theory.R")
library(showtext)
font_add_google("Nanum Gothic")
showtext_auto()

setwd("~/workspace/selection-against-introgression/datasets/swordtail_TOTO_ACUA/")

# 1. ===== Load Data =====
chrLen <- fread("xbir10x_chrlengths.txt")
#(chrLen[, floor(log2(V2/1000))])

loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]}


acua2006var <- data.table(loadFrom("ACUA_2006/wavelet_results.RData", "wavvar")); acua2006var[,year := "2006"]
acua2008var <- data.table(loadFrom("ACUA_2008/wavelet_results.RData", "wavvar")); acua2008var[,year := "2008"]
acua2013var <- data.table(loadFrom("ACUA_2013/wavelet_results.RData", "wavvar")); acua2013var[,year := "2013"]
acua2015var <- data.table(loadFrom("ACUA_2015/wavelet_results.RData", "wavvar")); acua2015var[,year := "2015"]
acua2018var <- data.table(loadFrom("ACUA_2018/wavelet_results.RData", "wavvar")); acua2018var[,year := "2018"]
wavvar <- rbindlist(list(acua2006var,acua2008var,acua2013var,acua2015var,acua2018var))

wavvarm <- melt(wavvar, id.vars = c("units", "year","level"))
wavvarm[, level := factor(level, levels = c(paste0("d", 1:11), "s8", "s9", "s11", "chr"))]


acua2006wc <- data.table(loadFrom("ACUA_2006/wavelet_results.RData", "wavcor")); acua2006wc[,year := "2006"]
acua2008wc <- data.table(loadFrom("ACUA_2008/wavelet_results.RData", "wavcor")); acua2008wc[,year := "2008"]
acua2013wc <- data.table(loadFrom("ACUA_2013/wavelet_results.RData", "wavcor")); acua2013wc[,year := "2013"]
acua2015wc <- data.table(loadFrom("ACUA_2015/wavelet_results.RData", "wavcor")); acua2015wc[,year := "2015"]
acua2018wc <- data.table(loadFrom("ACUA_2018/wavelet_results.RData", "wavcor")); acua2018wc[,year := "2018"]
wavcor <- rbindlist(list(acua2006wc,acua2008wc,acua2013wc,acua2015wc,acua2018wc))

acua2006rs <- data.table(loadFrom("ACUA_2006/wavelet_results.RData", "rsqrd")); acua2006rs[,year := "2006"]
acua2008rs <- data.table(loadFrom("ACUA_2008/wavelet_results.RData", "rsqrd")); acua2008rs[,year := "2008"]
acua2013rs <- data.table(loadFrom("ACUA_2013/wavelet_results.RData", "rsqrd")); acua2013rs[,year := "2013"]
acua2015rs <- data.table(loadFrom("ACUA_2015/wavelet_results.RData", "rsqrd")); acua2015rs[,year := "2015"]
acua2018rs <- data.table(loadFrom("ACUA_2018/wavelet_results.RData", "rsqrd")); acua2018rs[,year := "2018"]
rsqrd <- rbindlist(list(acua2006rs,acua2008rs,acua2013rs,acua2015rs,acua2018rs))


# ==== Power spectrum =====

# collapse scaling variances into one category
wavvarm_collapsed  <- rbind(wavvarm[grepl('s',level,fixed=T),
                                    .(level = 'scl', value = sum(value)),
                                    by = .(units,year,variable)],
                            wavvarm[!grepl('s',level,fixed=T), .(units,year,variable,level,value)]
)

# ----- plot on genetic map
wavvarm_collapsedG <- wavvarm_collapsed[units == 'genetic']
wavvarm_collapsedG[, level := droplevels(wavvarm_collapsedG$level)]
wavvarm_collapsedG[, level := factor(level, levels = c(paste0("d", 1:13), "scl", "chr"))]

lineDataG <- wavvarm_collapsedG[!level %in% c("scl", "chr")]


panel1 <- wavvarm_collapsedG[variable == "variance.meanFreq"] %>%
  ggplot(aes(x = level, y = value, group = year, color = year)) +
  geom_point(size=2) +
  geom_line(data = lineDataG[variable %in% c("variance.meanFreq")],
            aes(group = year), size=1, key_glyph = 'point') +
  labs(x = expression(Scale: log[2](Morgans)),
      title = 'A',
       y = "Variance",
       color = "Year", shape = "") +
  scale_x_discrete(breaks = c(paste0("d",1:11),"scl","chr"), labels = c(as.character(-12:-2),"scl", 'chrom')) +
  #scale_x_discrete(breaks = c(paste0("d",1:15),"s13", "s14", "s15", "chr"), labels = c(as.character(0:14),"13 (scaling var)", "14 (scaling var)", "15 (scaling var)", "chromosome")) +
  scale_shape_discrete(labels = c("Mean ancestry", "Individual ancestry"))+
  theme_classic() +
  scale_colour_viridis_d(option = 'G') +
  scale_y_continuous(label= function(x) {ifelse(x==0, "0", gsub("\\-0", "\\-", scientific_format()(x)))} ) +
  geom_segment(aes(x=.95,xend=11.05,y=-Inf,yend=-Inf),color="black")+
  theme(aspect.ratio = 1,
        legend.position = c(0.15,0.8),
        legend.key.size = unit(0.05, 'cm'),
        legend.text = element_text(size = 11),
        text = element_text(size=15, family = 'Nanum Gothic'),
        axis.ticks.x = element_line(size=0.5),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        axis.title.x = element_text(size = 13, vjust = 4),
        plot.title = element_text(hjust = -.1))

panel1


wavvar[units == 'genetic'] %>%
  ggplot(aes(x = level, y = variance.meanFreq, 
             group = year, color = year)) +
  geom_point(size=2) +
  geom_errorbar(data= wavvar[units == 'genetic' & !grepl('s', level)], 
                aes(ymin= variance.meanFreq - 1.96*variance.meanFreq.jack.se, 
                    ymax = variance.meanFreq + 1.96*variance.meanFreq.jack.se),
                linewidth=1, width = 0) +
  geom_line(data = wavvar[units == 'genetic' & grepl('d',level)],
            aes(group = year), size=1, key_glyph = 'point') +
  labs(x = expression(Scale: log[2](Morgans)),
       title = 'A',
       y = "Variance",
       color = "Year", shape = "") +
  scale_x_discrete(breaks = c(paste0("d",1:11),"s11","chr"), labels = c(as.character(-12:-2),"scl", 'chrom')) +
  scale_shape_discrete(labels = c("Mean ancestry", "Individual ancestry"))+
  theme_classic() +
  scale_colour_viridis_d(option = 'E') +
  scale_y_continuous(label= function(x) {ifelse(x==0, "0", gsub("\\-0", "\\-", scientific_format()(x)))} ) +
  geom_segment(aes(x=.95,xend=11.05,y=-Inf,yend=-Inf),color="black")+
  theme(aspect.ratio = 1,
        legend.position = c(0.15,0.8),
        legend.key.size = unit(0.05, 'cm'),
        legend.text = element_text(size = 11),
        text = element_text(size=15, family = 'Nanum Gothic'),
        axis.ticks.x = element_line(size=0.5),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        axis.title.x = element_text(size = 13, vjust = 4),
        plot.title = element_text(hjust = -.1))











# ----- compare to theory -----
# generation time? 2 gens per year vs 3 gens per year
g <- 3
a <- 114
wvtheory <- wavelet_variance_equilibrium(n.pop = 1000, n.sample = c(48,97)*2, unit.dist = 2^-12, gen = c(a-12*g, a), level = 1:11, alpha = 0.3)
wvtheory[, year := as.character((1/g)*gen + (2018 - a/g))]
wvtheory[, propvar := variance/sum(variance), by = .(year)]
wvtheory <- wvtheory[ (n.sample == 96 & year == 2006) | n.sample == 194 & year == 2018]

wvdetail <- wavvarm_collapsedG[variable == "variance.meanFreq" & grepl('d', level, fixed=T)]
wvdetail[, propvar := value/sum(value), by = year]

wvdetail[year %in% c(2006,2018)] %>%
  ggplot(aes(x = level, y = value, group = year, color = year)) +
  geom_point(size=3) +
  geom_line(aes(group = year), size=1.5) +
  geom_line(data = wvtheory, aes(x = level, y = variance), lty=11, size=1.5) +
  labs(x = expression(Scale: log[2](Morgans)),
       y = "Proportion of variance",
       color = "Year", shape = "") +
  scale_x_discrete(breaks = c(paste0("d",1:11)), labels = c(as.character(-12:-2))) +
  theme_classic() +
  #geom_segment(aes(x=.95,xend=13.05,y=-Inf,yend=-Inf),color="black")+
  theme(aspect.ratio = 1,
        text = element_text(size=15, family = "Nanum Gothic"),
        axis.ticks.x = element_line(size=1),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)))




# ---- physical
wavvarm_collapsedP <- wavvarm_collapsed[units == 'physical']
wavvarm_collapsedP[, level := droplevels(wavvarm_collapsedP$level)]
wavvarm_collapsedP[, level := factor(level, levels = c(paste0("d", 1:15), "scl", "chr"))]

lineDataP <- wavvarm_collapsedP[!level %in% c("scl", "chr")]

wavvarm_collapsedP[units == "physical" & variable == "variance.meanFreq"] %>%
  ggplot(aes(x = level, y = value, group = year, color = year)) +
  geom_point(size=3) +
  geom_line(data = lineDataP[variable %in% c("variance.meanFreq")],
            aes(group = year), size=02) +
  labs(#x = expression(Scale: log[2](Morgans)),
       x = expression(Scale: log[2](kb %*% 50)),
       y = "Variance",
       color = "Year", shape = "") +
  #scale_x_discrete(breaks = c(paste0("d",1:15),"s13","chr"), labels = c(as.character(-14:0),"scl", 'chrom')) +
  scale_x_discrete(breaks = c(paste0("d",1:15),"scl", "chr"), labels = c(as.character(1:15),"scl", "chrom")) +
  scale_shape_discrete(labels = c("Mean ancestry", "Individual ancestry"))+
  theme_classic() +
  scale_colour_viridis_d(option = 'E') +
  geom_segment(aes(x=.95,xend=9.05,y=-Inf,yend=-Inf),color="black")+
  theme(aspect.ratio = 1,
        text = element_text(size=15, family = 'Nanum Gothic'),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)))






# ===== plot wav cors ====

wavcor[, significant := ifelse(cor_jack - 1.96*cor_jack_se > 0, 1, 0)]


# ---- collapse  -----
wavcor_collapsed <- wavcor[!grepl('s',level,fixed=T)]

wavcor_collapsedG <- wavcor_collapsed[units == 'genetic']
wavcor_collapsedG[, level := droplevels(level)]
wavcor_collapsedG[, level := factor(level, levels = c(paste0("d", 1:13), 'scl', "chr"))]

wavcor_collapsedP <- wavcor_collapsed[units == 'physical']
wavcor_collapsedP[, level := droplevels(level)]
wavcor_collapsedP[, level := factor(level, levels = c(paste0("d", 1:15), 'scl', "chr"))]


# --- freq, recomb ----

# physical map
ggplot(wavcor_collapsedP[vars == 'meanFreq_r' & year == '2018'],
                  aes(x = level, y = cor_n, group = year))+ #, color = year)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=cor_n-1.96*cor_jack_se, ymax=cor_n + 1.96*cor_jack_se), width = 0.5, size=1)+
  scale_x_discrete(breaks = c(paste0("d",1:15),"chr"), labels = c(as.character(0:14),"chrom")) +
  labs(#x = expression(Scale: log[2]("Morgan")),
    x = expression(Scale: log[2](kb %*% 50)),
  y = "Correlation",
  title = "A") +
  geom_segment(aes(x=.95,xend=9.05,y=-Inf,yend=-Inf),color="black")+
  scale_color_viridis_d(option = 'E') +
  theme_classic() +
  theme(aspect.ratio = 1,
        text = element_text(size=15, family = "Nanum Gothic"),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1))

# genetic map
ggplot(wavcor_collapsedG[year == "2018" & 
                           vars == 'meanFreq_r'],
                   aes(x = level, y = cor_n, group = year))+ #, color = year)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=cor_n-1.96*cor_jack_se, ymax=cor_n + 1.96*cor_jack_se), width = 0, size=1)+
  scale_x_discrete(breaks = c(paste0("d",1:11),"chr"), labels = c(as.character(-12:-2),"chrom")) +
  labs(x = expression(Scale: log[2]("Morgan")),
    y = "Correlation",
    title = "A") +
  geom_segment(aes(x=.95,xend=11.05,y=-Inf,yend=-Inf),color="black")+

  theme_classic() +
  theme(aspect.ratio = 1,
        text = element_text(size=15, family = "Nanum Gothic"),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1))




# ----- mean freq, cds ----
ggplot(wavcor_collapsedP[year == "2018" & vars == 'meanFreq_cds_density'],
                   aes(x = level, y = cor_jack, group = year))+ #, color = year)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax=cor_jack + 1.96*cor_jack_se), width = 0, size=1)+
  scale_x_discrete(breaks = c(paste0("d",1:15),"chr"), labels = c(as.character(0:14), "chrom")) +
  #scale_x_discrete(breaks = c(paste0("d",1:15),"chr"), labels = c(as.character(0:14),  "chrom")) +
  labs(#x = expression(Scale: log[2]("Morgan")),
    x = expression(Scale: log[2](kb %*% 50)),
    # y = "Pearson cor (mean freq, CDS density)") +
    y = "Correlation",
    title = "A") +
  #y = "Pearson cor (log10 r, CDS density)" ) +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  theme_classic() +
  theme(aspect.ratio = 1,
        text = element_text(size=15, family = "Nanum Gothic"),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1))



ggplot(wavcor_collapsedG[year == "2018" & vars == 'meanFreq_cds_density'],
                   aes(x = level, y = cor_jack, group = year))+ #, color = year)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax=cor_jack + 1.96*cor_jack_se), width = 0, size=1)+
  scale_x_discrete(breaks = c(paste0("d",1:13),"chr"), labels = c(as.character(-14:-2),  "chrom")) +
  labs(x = expression(Scale: log[2]("Morgan")),
    y = "Correlation",
    title = "A") +
  geom_segment(aes(x=.95,xend=13.05,y=-Inf,yend=-Inf),color="black")+
  theme_classic() +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1))


# ----- rec, cds -----
ggplot(wavcor_collapsedP[year == "2018" & vars == 'r_cds_density'],
                   aes(x = level, y = cor_jack, group = year))+ #, color = year)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax=cor_jack + 1.96*cor_jack_se), width = 0, size=1)+
  scale_x_discrete(breaks = c(paste0("d",1:15),"chr"), labels = c(as.character(0:14), "chrom")) +
  #scale_x_discrete(breaks = c(paste0("d",1:15),"chr"), labels = c(as.character(0:14),  "chrom")) +
  labs(#x = expression(Scale: log[2]("Morgan")),
    x = expression(Scale: log[2]("1kb")),
    # y = "Pearson cor (mean freq, CDS density)") +
    y = "Correlation",
    title = "B") +
  #y = "Pearson cor (log10 r, CDS density)" ) +
  geom_segment(aes(x=.95,xend=9.05,y=-Inf,yend=-Inf),color="black")+
  theme_classic() +
  theme(aspect.ratio = 1,
        text = element_text(size=15, family = "Nanum Gothic"),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1))



ggplot(wavcor_collapsedG[year == "2018" & vars == 'r_cds_density'],
                   aes(x = level, y = cor_jack, group = year))+ #, color = year)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=cor_jack-1.96*cor_jack_se, ymax=cor_jack + 1.96*cor_jack_se), width = 0, size=1)+
  scale_x_discrete(breaks = c(paste0("d",1:13),"chr"), labels = c(as.character(-14:-2),  "chrom")) +
  labs(x = expression(Scale: log[2]("Morgan")),
    y = "Correlation",
    title = "B") +
  geom_segment(aes(x=.95,xend=13.05,y=-Inf,yend=-Inf),color="black")+
  theme_classic() +
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1))




# ===== plot r squared =====
rsqrdP <- rsqrd[units == 'physical' & !grepl("s",level,fixed=T)]
rsqrdP[, level := droplevels(level)]
rsqrdP[, level := factor(level, c(paste0("d", 1:15), "chr"))]

rsqrdG <- rsqrd[units == 'genetic' & !grepl("s",level,fixed=T)]
rsqrdG[, level := droplevels(level)]
rsqrdG[, level := factor(level, c(paste0("d", 1:13), "chr"))]


rsqrd_plot <- ggplot(rsqrdP[year == '2018' & model %in% 'r'],
       aes(x = level, y = rsqrd_jack, group = year)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin = rsqrd_jack - 1.96*rsqrd_jack_se, ymax= rsqrd_jack + 1.96*rsqrd_jack_se), width = 0, size=1) +

  #scale_x_discrete(breaks = c(paste0("d",1:13),"s13","chr"), labels = c(as.character(-14:-2), "-2s", "chromosome")) +
  scale_x_discrete(breaks = c(paste0("d",1:15),"chr"), labels = c(as.character(0:14), "chrom")) +
  labs(#x = expression(Scale: log[2]("Morgan")),
        x = expression(Scale: log[2]("1kb")),
       y = "R squared",
       title= "C") +
       #y = "Pearson cor (mean freq, log10 r)" ) +
  theme_classic() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.ticks.x = element_line(size=1),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1))

rsqrd_plot

rsqrdG[, significant := ifelse(rsqrd_jack-1.96*rsqrd_jack_se > 0, 1, 0)]

panel3 <- ggplot(rsqrdG[year == '2006' & model == "r" 
                              #& significant
                              ],
                     aes(x = level, y = rsqrd_jack)) + #, color = year)) + # y = cor_meanFreq_log10r)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin = rsqrd_jack - 1.96*rsqrd_jack_se, ymax= rsqrd_jack + 1.96*rsqrd_jack_se), width = 0, size=1) +
  scale_x_discrete(breaks = c(paste0("d",1:11),"chr"), labels = c(as.character(-12:-2), "chrom")) +
  labs(x = expression(Scale: log[2]("Morgan")),
    y = "R squared",
    title= "C") +
  theme_classic() +
  geom_segment(aes(x=.95,xend=11.05,y=-Inf,yend=-Inf),color="black")+
  theme(aspect.ratio = 1,
        text = element_text(size=15, family = "Nanum Gothic"),
        axis.ticks.x = element_line(size=0.5),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        axis.title.x = element_text(size = 13, vjust = 4),
        
        plot.title = element_text(hjust = -.1)
        )

panel3




# ===== Contribution to correlation =====

allwav <- merge(wavcor, wavvar, by  = c("units", "level", "year"))
allwav[vars== 'meanFreq_r', contribution := cor_n*sqrt(propvar.meanFreq*propvar.r)]

# collapse scaling
allwav_collapsed <- rbind(
  allwav[vars == 'meanFreq_r' & grepl('s',level,fixed=T), .(level = 'scl', contribution=sum(contribution)), by =.(units, year)],
      allwav[vars == 'meanFreq_r' & !grepl('s',level,fixed=T), .(units, level, year, contribution)]
  )

allwav_collapsed[, normcor := contribution/sum(contribution), by = .(units, year)]

allwav_collapsedP <- allwav_collapsed[units == 'physical']
allwav_collapsedP[, level := droplevels(level)]
allwav_collapsedP[, level := factor(level, levels = c(paste0("d", 1:15), "scl", "chr"))]

allwav_collapsedG <- allwav_collapsed[units == 'genetic']
allwav_collapsedG[, level := droplevels(level)]
allwav_collapsedG[, level := factor(level, levels = c(paste0("d", 1:13), "scl", "chr"))]


library(scales)

ggplot(allwav_collapsedP[year == "2018"]) +
  geom_bar(aes(fill = level, x = year, y = normcor), position = "stack", stat = "identity") +
  scale_fill_manual(values = c(viridis_pal()(9), 'lightgrey', 'darkgrey'), labels  = c(as.character(0:8), "scl", "chrom")) +
  labs(x = "",
       y = "Contribution to correlation",
       fill = expression(Scale: log[2](kb %*% 50)), title = "B") +
  theme_classic() +
  theme(aspect.ratio = 5,
        text = element_text(size=15),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        legend.key.size = unit(0.5, 'cm'),
        plot.title = element_text(hjust = -1))




panel2P <- ggplot(allwav_collapsedP, 
       aes(x = level, y = normcor, color = year, group = year)) + 
  geom_point() +  
  geom_line(data= allwav_collapsedP[grepl('d', level, fixed=T)]) + 
  scale_colour_viridis_d(option = 'E')  +
  scale_x_discrete(breaks = c(paste0("d",1:15),"scl","chr"), labels = c(as.character(0:14),"scl", "chrom")) +
  theme_classic() +
  geom_segment(aes(x=.95,xend=9.05,y=-Inf,yend=-Inf),color="black")+
  labs( x= expression(Scale: log[2](kb %*% 50)), y = 'Normalized contribution\n to overall correlation') +
  theme_classic() +
  theme(aspect.ratio = 1,
        text = element_text(size=15, family = "Nanum Gothic"),
        axis.ticks.x = element_line(size=0.5),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        plot.title = element_text(hjust = -.1))
  
panel2P

panel2G <- ggplot(allwav_collapsedG, 
       aes(x = level, y = normcor, color = year, group = year)) + 
  geom_point(size=2) +  
  geom_line(data= allwav_collapsedG[grepl('d', level, fixed=T)]) + 
  scale_colour_viridis_d(option = 'G')  +
  scale_x_discrete(breaks = c(paste0("d",1:11),"scl","chr"), labels = c(as.character(-12:-2),"scl", "chrom")) +
  theme_classic() +
  scale_y_continuous(limits = c(-.01, 0.3), labels = c(0,10,20,30)) +
  geom_segment(aes(x=.95,xend=11.05,y=-Inf,yend=-Inf),color="black")+
  labs( x= expression(Scale: log[2](Morgan)), 
        y = 'Percent contribution\nto overall correlation',
        title = 'B') +
  theme_classic() +
  theme(aspect.ratio = 1,
        legend.position = 'None',
        text = element_text(size=15, family = "Nanum Gothic"),
        axis.ticks.x = element_line(size=0.5),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(hjust=.4,margin=margin(r=10)),
        axis.title.x = element_text(size = 13, vjust = 4),
        plot.title = element_text(hjust = -.1))

panel2G

panel1 + panel2G
panel3


