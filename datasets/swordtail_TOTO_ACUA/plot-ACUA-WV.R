library(data.table)
library(tidyverse)

loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 

acua2006 <- loadFrom("ACUA_2006/wvFinal.RData", "wvFinal"); acua2006[,year := "2006"]
acua2008 <- loadFrom("ACUA_2008/wvFinal.RData", "wvFinal"); acua2008[,year := "2008"]
acua2013 <- loadFrom("ACUA_2013/wvFinal.RData", "wvFinal"); acua2013[,year := "2013"]
acua2015 <- loadFrom("ACUA_2015/wvFinal.RData", "wvFinal"); acua2015[,year := "2015"]
acua2018 <- loadFrom("ACUA_2018/wvFinal.RData", "wvFinal"); acua2018[,year := "2018"]

allWV <- rbindlist(list(acua2006,acua2008,acua2013,acua2015,acua2018))

# take average over chromosomes
wvGnomWide <- allWV[, mean(variance, na.rm=T), by = .(decomp,year,scale)]
setnames(wvGnomWide, "V1", "variance")
wvGnomWide[,propVar := variance/sum(variance),by=.(decomp,year)]

# plot this mf
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

# plot proportion
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


