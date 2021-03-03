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

# load wav cor files
acua2006wc <- loadFrom("ACUA_2006/wavcorFinal.RData", "wavcorFinal"); acua2006wc[,year := "2006"]
acua2008wc <- loadFrom("ACUA_2008/wavcorFinal.RData", "wavcorFinal"); acua2008wc[,year := "2008"]
acua2013wc <- loadFrom("ACUA_2013/wavcorFinal.RData", "wavcorFinal"); acua2013wc[,year := "2013"]
acua2015wc <- loadFrom("ACUA_2015/wavcorFinal.RData", "wavcorFinal"); acua2015wc[,year := "2015"]
acua2018wc <- loadFrom("ACUA_2018/wavcorFinal.RData", "wavcorFinal"); acua2018wc[,year := "2018"]
allWC <- rbindlist(list(acua2006wc,acua2008wc,acua2013wc,acua2015wc,acua2018wc))
allWC <- melt(allWC, id.vars = c("signal", "year", "chr"), 
              variable.name = "scale", value.name = "correlation")

# wav var: take average over chromosomes
wvGnomWide <- allWV[, mean(variance, na.rm=T), by = .(decomp,year,scale)]
setnames(wvGnomWide, "V1", "variance")
wvGnomWide[,propVar := variance/sum(variance),by=.(decomp,year)]

# wav var: plot
wvGnomWide[decomp == "mean_individual"] %>% ggplot(aes(x = scale, y = variance, group = interaction(decomp, year), color = year)) +
  geom_point() +
  geom_line(aes(linetype=decomp), size=0.5) + 
  labs(x = expression(Scale: log[2](Morgans)), 
       y = "Wavelet variance",
       color = "Year", linetype = "Signal") +
  scale_x_discrete(labels = -14:0) + 
  scale_linetype_discrete(labels = c("Individual","Pop mean"))+
  theme_classic() +
  scale_colour_viridis_d()

# wav var: plot proportion
wvGnomWide %>% ggplot(aes(x = scale, y = propVar, group = interaction(decomp, year), color = year)) +
  geom_point(size=1.8) +
  geom_line(aes(linetype=decomp), size=0.5) + 
  labs(x = expression(Scale: log[2](Morgans)), 
       y = "Proportion of variance",
       color = "Year", linetype = "Signal") +
  scale_x_discrete(labels = -14:0) + 
  scale_linetype_discrete(labels = c("Individual","Pop mean"))+
  theme_classic() +
  scale_colour_viridis_d()

# wav cor: take average over chromosomes
wavcorGnomWide <- allWC[, mean(correlation, na.rm=T), by = .(signal,year,scale)]
setnames(wavcorGnomWide, "V1", "correlation")
wavcorGnomWide[, scale := factor(gsub("d", "",wavcorGnomWide$scale), levels = rev(1:15))]
wavcorGnomWide$scale
wavcorGnomWide[, norm_cor := correlation/ sum(correlation, na.rm=T), by = .(signal,year)]
# stacked barplot
ggplot(wavcorGnomWide[signal == "pop_mean",]) +
  geom_bar(aes(fill = scale, x = year, y = norm_cor), position = "stack", stat = "identity", color = "black") +
  scale_fill_viridis_d(option = "plasma") + 
  labs(x = "Year", 
       y = "Normalized wavelet correlation",
       fill = expression(Scale: log[2](kb))) + 
  theme_classic()


