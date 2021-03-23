library(data.table)
library(tidyverse)

loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 

# load wav variance files
acua2006wv <- loadFrom("ACUA_2006/wvFinal.RData", "wvFinal"); acua2006wv[,year := "2006"]
acua2008wv <- loadFrom("ACUA_2008/wvFinal.RData", "wvFinal"); acua2008wv[,year := "2008"]
acua2013wv <- loadFrom("ACUA_2013/wvFinal.RData", "wvFinal"); acua2013wv[,year := "2013"]
acua2015wv <- loadFrom("ACUA_2015/wvFinal.RData", "wvFinal"); acua2015wv[,year := "2015"]
acua2018wv <- loadFrom("ACUA_2018/wvFinal.RData", "wvFinal"); acua2018wv[,year := "2018"]
allWV <- rbindlist(list(acua2006wv,acua2008wv,acua2013wv,acua2015wv,acua2018wv))

# load chr mean files
acua2006chrVar <- loadFrom("ACUA_2006/wvFinal.RData", "chrVar"); acua2006chrVar[,year := "2006"]
acua2008chrVar <- loadFrom("ACUA_2008/wvFinal.RData", "chrVar"); acua2008chrVar[,year := "2008"]
acua2013chrVar <- loadFrom("ACUA_2013/wvFinal.RData", "chrVar"); acua2013chrVar[,year := "2013"]
acua2015chrVar <- loadFrom("ACUA_2015/wvFinal.RData", "chrVar"); acua2015chrVar[,year := "2015"]
acua2018chrVar <- loadFrom("ACUA_2018/wvFinal.RData", "chrVar"); acua2018chrVar[,year := "2018"]
allChrVar <- rbindlist(list(acua2006chrVar,acua2008chrVar,acua2013chrVar,acua2015chrVar,acua2018chrVar))

# for plotting chromosome-level  variance alongside wavelet variance
#allWV[, level := "within-chrom"]
allWV[, scale := as.numeric(scale)]
allWV[decomp == "mean_individual", decomp := "individual"]
#allChrVar[, level := "among-chrom"] 
allChrVar[, scale := 17]
setnames(allChrVar, "chrVar", "variance")

# load wav cor files
acua2006wc <- loadFrom("ACUA_2006/wavcorFinal.RData", "wavcorFinal"); acua2006wc[,year := "2006"]
acua2008wc <- loadFrom("ACUA_2008/wavcorFinal.RData", "wavcorFinal"); acua2008wc[,year := "2008"]
acua2013wc <- loadFrom("ACUA_2013/wavcorFinal.RData", "wavcorFinal"); acua2013wc[,year := "2013"]
acua2015wc <- loadFrom("ACUA_2015/wavcorFinal.RData", "wavcorFinal"); acua2015wc[,year := "2015"]
acua2018wc <- loadFrom("ACUA_2018/wavcorFinal.RData", "wavcorFinal"); acua2018wc[,year := "2018"]
allWC <- rbindlist(list(acua2006wc,acua2008wc,acua2013wc,acua2015wc,acua2018wc))

# wav var: take average over chromosomes
wvGnomWide <- allWV[, weighted.mean(variance, weight, na.rm=T), by = .(decomp, scale, year)]
setnames(wvGnomWide, "V1", "variance")

wvGnomWide <- merge(wvGnomWide, allChrVar, all=TRUE)
wvGnomWide[,propVar := variance/sum(variance),by=.(decomp,year)]
#wvGnomWide[, level := factor(level, levels = c("within-chrom", "among-chrom"))]


lineData <- wvGnomWide[scale < 17]
# wav var: plot

wvGnomWide[] %>% ggplot(aes(x = scale, y = variance, group = interaction(decomp, year), color = year)) +
  geom_point(aes(shape = decomp),size=2.2) +
  geom_line(data = lineData, size=0.5,linetype=2) + 
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

 e# 
  


# wav var: plot proportion
wvGnomWide[decomp=="pop_mean"] %>% ggplot(aes(x = scale, y = propVar, group = interaction(decomp, year), color = year)) +
  geom_point(size=1.3) +
  geom_line(aes(linetype=decomp), size=0.5, linetype=2) + 
  labs(x = expression(Scale: log[2](Morgans)), 
       y = "Proportion of variance",
       color = "Year", group = "Signal") +
  scale_x_discrete(labels = -14:0) + 
  scale_linetype_discrete(labels = c("Individual","Pop mean"))+
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=1,xend=17,y=-Inf,yend=-Inf))+
  theme(axis.line=element_blank())

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


