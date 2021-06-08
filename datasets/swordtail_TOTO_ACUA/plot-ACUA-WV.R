library(data.table)
library(tidyverse)

# 1. ===== Load Data =====

loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 

# 1.1. ----- load wavelet variance files -----
acua2006wv <- loadFrom("ACUA_2006/anc_rec_cds_varDecomp.RData", "wvFinalAll"); acua2006wv[,year := "2006"]
acua2008wv <- loadFrom("ACUA_2008/anc_rec_cds_varDecomp.RData", "wvFinalAll"); acua2008wv[,year := "2008"]
acua2013wv <- loadFrom("ACUA_2013/anc_rec_cds_varDecomp.RData", "wvFinalAll"); acua2013wv[,year := "2013"]
acua2015wv <- loadFrom("ACUA_2015/anc_rec_cds_varDecomp.RData", "wvFinalAll"); acua2015wv[,year := "2015"]
acua2018wv <- loadFrom("ACUA_2018/anc_rec_cds_varDecomp.RData", "wvFinalAll"); acua2018wv[,year := "2018"]
allWV <- rbindlist(list(acua2006wv,acua2008wv,acua2013wv,acua2015wv,acua2018wv))

# 1.2. ----- load chromosome-level variance files -----
acua2006chrVar <- loadFrom("ACUA_2006/anc_rec_cds_varDecomp.RData", "chrVarAll"); acua2006chrVar[,year := "2006"]
acua2008chrVar <- loadFrom("ACUA_2008/anc_rec_cds_varDecomp.RData", "chrVarAll"); acua2008chrVar[,year := "2008"]
acua2013chrVar <- loadFrom("ACUA_2013/anc_rec_cds_varDecomp.RData", "chrVarAll"); acua2013chrVar[,year := "2013"]
acua2015chrVar <- loadFrom("ACUA_2015/anc_rec_cds_varDecomp.RData", "chrVarAll"); acua2015chrVar[,year := "2015"]
acua2018chrVar <- loadFrom("ACUA_2018/anc_rec_cds_varDecomp.RData", "chrVarAll"); acua2018chrVar[,year := "2018"]
allChrVar <- rbindlist(list(acua2006chrVar,acua2008chrVar,acua2013chrVar,acua2015chrVar,acua2018chrVar))

# combine wavelet and chromosome data - add chromosome as fake scale number for display
allChrVar[, scale := as.character(17)][, scale := as.numeric(scale)]
allWV[, scale := as.numeric(scale)]
allVar <- merge(allWV, allChrVar, all=TRUE)

# 1.3. ----- load cor decomp files -----
acua2006cor <- loadFrom("ACUA_2006/anc_rec_cds_corDecomp.RData", "allVarsCorDecomp"); acua2006cor[,year := "2006"]
acua2008cor <- loadFrom("ACUA_2008/anc_rec_cds_corDecomp.RData", "allVarsCorDecomp"); acua2008cor[,year := "2008"]
acua2013cor <- loadFrom("ACUA_2013/anc_rec_cds_corDecomp.RData", "allVarsCorDecomp"); acua2013cor[,year := "2013"]
acua2015cor <- loadFrom("ACUA_2015/anc_rec_cds_corDecomp.RData", "allVarsCorDecomp"); acua2015cor[,year := "2015"]
acua2018cor <- loadFrom("ACUA_2018/anc_rec_cds_corDecomp.RData", "allVarsCorDecomp"); acua2018cor[,year := "2018"]
allCor <- rbindlist(list(acua2006cor,acua2008cor,acua2013cor,acua2015cor,acua2018cor))
allCor[, scale := as.character(scale)]
allCor[scale == 'chr', scale := '17']
allCor[, scale := as.numeric(scale)]



# 2. ===== Ancestry Variance Decomposition, Genetic Scale =========
lineData <- allVar[scale < 17 & decomp == "pop_mean"]

allVar[distance == "genetic" & decomp == "pop_mean"] %>% ggplot(aes(x = scale, y = anc_variance, group = interaction(decomp, year), color = year)) +
  geom_point(aes(shape = decomp),size=2.2) +
  geom_line(data = lineData[distance == "genetic"], size=0.5, linetype=2) + 
  labs(x = expression(Scale: log[2](Morgans)), 
       y = "Variance",
       color = "Year", shape = "Signal") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(-15:-1),"chromosome\nlevel")) + 
  scale_shape_discrete(labels = c("Mean\nancestry"))+
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black")+
  theme(text = element_text(size=15),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=.4,margin=margin(t=-20)))

# 2.1. ----- Plot proportion ---------------
allVar[, propVar := anc_variance/sum(anc_variance), by = .(decomp,distance,year)]
lineDataProp <- allVar[scale < 17 & decomp == "pop_mean"]
  
allVar[decomp == "pop_mean" & distance == "genetic"] %>% ggplot(aes(x = scale, y = propVar, group = interaction(decomp, year), color = year)) +
  geom_point(aes(shape = decomp),size=2.2) +
  geom_line(data = lineDataProp[decomp == "pop_mean" & distance == "genetic"], size=0.5,linetype=2) + 
  labs(x = expression(Scale: log[2](Morgans)), 
       y = "Proportion of variance",
       color = "Year", shape = "Signal") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(-15:-1),"chromosome\nlevel")) + 
  scale_shape_discrete(labels = c("Individual","Pop mean"))+
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black")+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5),#,
        axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=.4,margin=margin(t=-20)))


# 2.2. ----- individual chromosomes -----

acua2006_wvChrG <- loadFrom("ACUA_2006/anc_rec_cds_varDecomp.RData", "wvChrsAll"); acua2006_wvChrG[,year := "2006"]

acua2006_wvChrG[distance == "genetic"] %>%
  ggplot(aes(x = scale, y = anc_variance, group = chr, color = chr)) + 
  geom_point() + geom_line() + 
  theme(legend.position = "none")


# 2.3. ----- Chromosomes with mapped incompatibilities-------

# which of these chromosomes is the outlier with huge large-scale variance?
acua2006_wvChrG[scale == 14][anc_variance==max(anc_variance, na.rm=T)]

selChrs <- c("ScyDAA6-1439-HRSCAF-1708",
             "ScyDAA6-1592-HRSCAF-1896",
             "ScyDAA6-1934-HRSCAF-2318",
             "ScyDAA6-2393-HRSCAF-2888",
             "ScyDAA6-5078-HRSCAF-5686",
             "ScyDAA6-695-HRSCAF-847",
             "ScyDAA6-7-HRSCAF-50")

acua2006_wvChrG[chr %in% selChrs, BDMI := "yes"]
acua2006_wvChrG[!chr %in% selChrs, BDMI := "no"]

acua2006_wvChrG[distance == "genetic"] %>%
  ggplot(aes(x = scale, y = anc_variance, group = chr, color = BDMI)) + 
  geom_point() + geom_line() #+ 
  #theme(legend.position = "none")


# 3. ===== Signal Variance Decompositions, Physical Scale =========

# 3.1. ----- Ancestry -----
allVar[distance == "physical"] %>% ggplot(aes(x = scale, y =anc_variance, group = interaction(decomp, year), color = year)) +
  geom_point(aes(shape = decomp),size=2.2) +
  geom_line(data = lineData[distance == "physical"], size=0.5,linetype=2) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Variance",
       color = "Year", shape = "Signal") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(0:14),"chromosome\nlevel")) + 
  scale_shape_discrete(labels = c("Individual","Pop mean"))+
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black")+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5),#,
        axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=.4,margin=margin(t=-20)))

# 3.1.1. ----- Ancestry: proportion of variance -----

allVar[decomp == "mean_individual" & distance == "physical"] %>% ggplot(aes(x = scale, y = propVar, group = interaction(decomp, year), color = year)) +
  geom_point(aes(shape = decomp),size=2.2) +
  geom_line(data = lineDataProp[decomp == "mean_individual" & distance == "physical"], size=0.5,linetype=2) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Proportion of variance",
       color = "Year", shape = "Signal") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(0:14),"chromosome\nlevel")) + 
  scale_shape_discrete(labels = c("Individual","Pop mean"))+
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black")+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5),#,
        axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=.4,margin=margin(t=-20)))

# 3.2 ----- Rec Var Decomp -----

allVar[distance == "physical" & year == "2018" & decomp == "pop_mean"] %>% ggplot(aes(x = scale, y =rec_variance)) +
  geom_point(size=2.2) +
  geom_line(data = lineData[distance == "physical" & decomp == "pop_mean"], size=0.5,linetype=2,group=year) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Variance",
       color = "Year") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(0:14),"chromosome\nlevel")) + 
  theme_classic() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black")+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5),#,
        axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=.4,margin=margin(t=-20)))

# 3.3 ----- CDS Var Decomp -----

allVar[distance == "physical" & year == "2018" & decomp == "pop_mean"] %>% ggplot(aes(x = scale, y =cds_variance)) +
  geom_point(size=2.2) +
  geom_line(data = lineData[distance == "physical" & decomp == "pop_mean"], size=0.5,linetype=2,group=year) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Variance",
       color = "Year") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(0:14),"chromosome\nlevel")) + 
  theme_classic() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black")+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5),#,
        axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=.4,margin=margin(t=-20)))


# 3.4 ----- CDS/CM Var Decomp -----

allVar[distance == "physical" & year == "2018" & decomp == "pop_mean"] %>% ggplot(aes(x = scale, y =cdsPerCm_variance)) +
  geom_point(size=2.2) +
  geom_line(data = lineData[distance == "physical" & decomp == "pop_mean"], size=0.5,linetype=2,group=year) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Variance",
       color = "Year") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(0:14),"chromosome\nlevel")) + 
  theme_classic() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black")+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5),#,
        axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=.4,margin=margin(t=-20)))


# 4. ===== Correlation Decompositions ===========
lineDataCor <- allCor[scale < 17]


# 4.1. ----- Rec x Ancestry -----
allCor %>%
  ggplot(aes(x = scale, y = rec_anc_cor, group = year, color = year)) + 
  geom_point() + 
  geom_line(data = lineDataCor, size=0.5,linetype=2) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Pearson correlation coefficient",
       color = "Year") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(0:14),"chromosome\nlevel")) + 
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 10)) +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black") +
  theme(axis.line.x = element_blank())

# 4.2. ----- CDS x Ancestry -----
allCor %>%
  ggplot(aes(x = scale, y = cds_anc_cor, group = year, color = year)) + 
  geom_point() + 
  geom_line(data = lineDataCor, size=0.5,linetype=2) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Pearson correlation coefficient",
       color = "Year") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(0:14),"chromosome\nlevel")) + 
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black") +
  theme(axis.line.x = element_blank())


# 4.3. ----- CDS x Rec -----
allCor %>%
  ggplot(aes(x = scale, y = cds_rec_cor)) + 
  geom_point() + 
  geom_line(data = lineDataCor, size=0.5,linetype=2) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Pearson correlation coefficient",
       color = "Year") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(0:14),"chromosome\nlevel")) + 
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black") +
  theme(axis.line.x = element_blank())


# 4.3. ----- cdsPerCm x Anc -----
allCor %>%
  ggplot(aes(x = scale, y = cdsPerCm_anc_cor)) + 
  geom_point() + 
  geom_line(data = lineDataCor, size=0.5,linetype=2) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Pearson correlation coefficient",
       color = "Year") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(0:14),"chromosome\nlevel")) + 
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black") +
  theme(axis.line.x = element_blank())


head(allCor)
head(allVar)

mrgd <- merge(allCor, allVar[distance == "physical" & decomp == "pop_mean"], by = c("year", "scale"))

allVar[year == 2018 & scale == 1 & distance == "physical"]
head(mrgd[year == 2018])
head(mrgd[year == 2015])


mrgd <- mrgd[scale < 17]
mrgd[, contribution := sqrt(anc_variance*rec_variance)*rec_anc_cor, by = .(year,scale)]
mrgd[, norm_cor := contribution/sum(contribution), by = .(year)]


mrgd[scale < 17] %>%
  ggplot() + 
  geom_bar(aes(fill = reorder(as.factor(scale), desc(as.factor(scale))), x = year, y = norm_cor),
               position = "stack", stat = "identity", color = "black") + 
  scale_fill_viridis_d(option = "plasma") + 
  labs(x = "Year", 
       y = "Normalized wavelet correlation",
       fill = expression(Scale: log[2](kb))) + 
  theme(text = element_text(size = 15))
  theme_classic()



# wav cor: take average over chromosomes
wavcorGnomWide <- allWC[, weighted.mean(variance, weight, na.rm=T), by = .(signal, scale, year)]
setnames(wavcorGnomWide, "V1", "correlation")

# wavelet correlations
# wav var: plot proportion
wavcorGnomWide[signal == "pop_mean"] %>% ggplot(aes(x = year, y = correlation, group = interaction(signal, scale), color = scale)) +
  geom_point(size=1.8) +
  geom_line(aes(linetype=signal), size=0.5) + 
  labs(x = "Year", 
       y = "Wavelet correlation",
       color = "scale", linetype = "Signal") +
  #scale_linetype_discrete(labels = c("Individual","Pop mean"))+
  theme_classic() +
  scale_colour_viridis_d()

wavcorGnomWide[, norm_cor := correlation/ sum(correlation, na.rm=T), by = .(signal,year)]
# stacked barplot
ggplot(wavcorGnomWide[signal == "pop_mean",]) +
  geom_bar(aes(fill = scale, x = year, y = norm_cor), position = "stack", stat = "identity", color = "black") +
  scale_fill_viridis_d(option = "plasma") + 
  labs(x = "Year", 
       y = "Normalized wavelet correlation",
       fill = expression(Scale: log[2](kb))) + 
  theme_classic()


