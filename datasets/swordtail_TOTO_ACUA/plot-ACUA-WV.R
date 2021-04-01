library(data.table)
library(tidyverse)

loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 

# load wavelet variance files
acua2006wv <- loadFrom("ACUA_2006/varDecompAll.RData", "wvFinalAll"); acua2006wv[,year := "2006"]
acua2008wv <- loadFrom("ACUA_2008/varDecompAll.RData", "wvFinalAll"); acua2008wv[,year := "2008"]
acua2013wv <- loadFrom("ACUA_2013/varDecompAll.RData", "wvFinalAll"); acua2013wv[,year := "2013"]
acua2015wv <- loadFrom("ACUA_2015/varDecompAll.RData", "wvFinalAll"); acua2015wv[,year := "2015"]
acua2018wv <- loadFrom("ACUA_2018/varDecompAll.RData", "wvFinalAll"); acua2018wv[,year := "2018"]
allWV <- rbindlist(list(acua2006wv,acua2008wv,acua2013wv,acua2015wv,acua2018wv))

# load chrromosome-level variance files
acua2006chrVar <- loadFrom("ACUA_2006/varDecompAll.RData", "chrVarAll"); acua2006chrVar[,year := "2006"]
acua2008chrVar <- loadFrom("ACUA_2008/varDecompAll.RData", "chrVarAll"); acua2008chrVar[,year := "2008"]
acua2013chrVar <- loadFrom("ACUA_2013/varDecompAll.RData", "chrVarAll"); acua2013chrVar[,year := "2013"]
acua2015chrVar <- loadFrom("ACUA_2015/varDecompAll.RData", "chrVarAll"); acua2015chrVar[,year := "2015"]
acua2018chrVar <- loadFrom("ACUA_2018/varDecompAll.RData", "chrVarAll"); acua2018chrVar[,year := "2018"]
allChrVar <- rbindlist(list(acua2006chrVar,acua2008chrVar,acua2013chrVar,acua2015chrVar,acua2018chrVar))


# load wav cor files
acua2006wc <- loadFrom("ACUA_2006/wavcorFinal.RData", "wavcorFinal"); acua2006wc[,year := "2006"]
acua2008wc <- loadFrom("ACUA_2008/wavcorFinal.RData", "wavcorFinal"); acua2008wc[,year := "2008"]
acua2013wc <- loadFrom("ACUA_2013/wavcorFinal.RData", "wavcorFinal"); acua2013wc[,year := "2013"]
acua2015wc <- loadFrom("ACUA_2015/wavcorFinal.RData", "wavcorFinal"); acua2015wc[,year := "2015"]
acua2018wc <- loadFrom("ACUA_2018/wavcorFinal.RData", "wavcorFinal"); acua2018wc[,year := "2018"]
allWC <- rbindlist(list(acua2006wc,acua2008wc,acua2013wc,acua2015wc,acua2018wc))

# combine wavelet and chromosome data - add chromosome as fake scale number for display
allChrVar[, scale := as.character(17)][, scale := as.numeric(scale)]
allWV[, scale := as.numeric(scale)]
allVar <- merge(allWV, allChrVar, all=TRUE)


# Plot variance decomposition genetic scale =========
lineData <- allVar[scale < 17]

allVar[distance == "genetic"] %>% ggplot(aes(x = scale, y = variance, group = interaction(decomp, year), color = year)) +
  geom_point(aes(shape = decomp),size=2.2) +
  geom_line(data = lineData[distance == "genetic"], size=0.5,linetype=2) + 
  labs(x = expression(Scale: log[2](Morgans)), 
       y = "Variance",
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

# Plot proportion ---------------
allVar[, propVar := variance/sum(variance), by = .(decomp,distance,year)]
lineDataProp <- allVar[scale < 17]
  
allVar[decomp == "mean_individual" & distance == "genetic"] %>% ggplot(aes(x = scale, y = propVar, group = interaction(decomp, year), color = year)) +
  geom_point(aes(shape = decomp),size=2.2) +
  geom_line(data = lineDataProp[decomp == "mean_individual" & distance == "genetic"], size=0.5,linetype=2) + 
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

# Plot variance decomp on physical scale =========

allVar[distance == "physical"] %>% ggplot(aes(x = scale, y = variance, group = interaction(decomp, year), color = year)) +
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

# Plot proportion


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


