library(data.table)
library(tidyverse)

# 1. ===== Load Data =====
chrLen <- fread("xbir10x_chrlengths.txt")
(chrLen[, floor(log2(V2/1000))])

loadFrom=function(file, name){e=new.env();load(file,env=e);e[[name]]} 


acua2006var <- data.table(loadFrom("ACUA_2006/ancestry_allVarDecomp.RData", "allVarDecomp")); acua2006var[,year := "2006"]
acua2008var <- data.table(loadFrom("ACUA_2008/ancestry_allVarDecomp.RData", "allVarDecomp")); acua2008var[,year := "2008"]
acua2013var <- data.table(loadFrom("ACUA_2013/ancestry_allVarDecomp.RData", "allVarDecomp")); acua2013var[,year := "2013"]
acua2015var <- data.table(loadFrom("ACUA_2015/ancestry_allVarDecomp.RData", "allVarDecomp")); acua2015var[,year := "2015"]
acua2018var <- data.table(loadFrom("ACUA_2018/ancestry_allVarDecomp.RData", "allVarDecomp")); acua2018var[,year := "2018"]
allVar <- rbindlist(list(acua2006var,acua2008var,acua2013var,acua2015var,acua2018var))


load("ACUA_2018/wavcor_results.RData")
ggplot(wavcor[!grepl("15", level)], aes(x = level, y = pearson.cor.meanFreq.cmTr)) + geom_point()
data.table(loadFrom("ACUA_2018/wavcor_results.RData")); acua2006var[,year := "2006"]


# =====plot simulation and real data together ===== 
simVarDecomp <- loadFrom("varDecomp_simulation_rhoCap0.005.RData", "allVarDecomp")
simVarDecomp <- simVarDecomp[N_diploids == 100]

allVar[, data_source := "real"]

# average across replicates
simVarDecomp <- simVarDecomp[, .(propVar = mean(propVar), variance = mean(variance)), by = .(gen, level, N_diploids,data_source, units)]


allVar[, gen := year]
sim_and_real <- rbind(simVarDecomp[, .(level,propVar,variance,N_diploids,gen,data_source,units)], allVar[signal == "mean", .(level,propVar,variance,N_diploids=1,gen, data_source,units)])

sim_and_real[, alpha:= ifelse(data_source== "real", 1, 0.75)]
sim_and_real[, level := factor(level, levels = c(paste0("d",1:15), paste0("s", 12:15), "chr"))]
sim_and_real[, levels(level)]

lineData <- sim_and_real[!level %in% c("s12","s13","s14","s15","chr")]
lineData[, levels(level)]


ggplot(sim_and_real[units == "genetic"], aes(x = level, 
                         y = variance, 
                         group = gen,
                         color = interaction(data_source, gen))) +
  geom_point(aes(shape = data_source)) + 
  geom_line(data = lineData[units == "genetic"], aes(group=gen, linetype = data_source), size=0.5) + 
  labs(x = expression(Scale: log[2](Morgans)), 
       # x = expression(Scale: log[2] (kbp)),
       y = "Ancestry variance",
       color = "") +
  #scale_x_discrete(breaks = c(paste0("d",1:14),"s12","s13","s14","chr"), labels = c(as.character(-14:-1), paste0(-3:-1," (scaling var)"), "chromosome")) + 
  scale_x_discrete(breaks = c(paste0("d",1:13),"s13","chr"), labels = c(as.character(-14:-2), " -2 (scaling var)", "chromosome")) + 
  
  #scale_x_discrete(breaks = c(paste0("d",1:15),"s13","s14","s15","chr"), labels = c(as.character(1:15), paste0(13:15," (scaling var)"), "chromosome")) + 
  theme_classic() +
  scale_colour_viridis_d() +
  #geom_segment(aes(x=.95,xend=13.05,y=-Inf,yend=-Inf),color="black") + 
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        #axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=0.4,margin=margin(t=-5)))


# ===== correlations =====

acua2018cor <- data.table(loadFrom("ACUA_2018/allCorDecomp_anc_rec_cds.RData", "allCorDecomp_anc_rec_cds"))

ggplot(acua2018cor[!grepl("s", level)], 
       aes(x=level, y=rec_anc_cor)) +
  geom_point() + 
  scale_x_discrete(breaks = c(paste0("d",1:15),"chr"), labels = c(as.character(1:15),  "chrom")) +
  theme_classic() + 
  labs(x = expression(Scale: log[2](kb)), y = "Correlation")
  





# ===== plot real data ====
levels(allVar$level)
lineData <- allVar[!level %in% c("s13", "s14", "chr")]

allVar[units == "genetic"] %>% ggplot(aes(x = level, y = variance, group = interaction(signal, year), color = year)) +
  geom_point(aes(shape = signal), size=2.2) +
  geom_line(data = lineData[units == "genetic"], size=0.5, linetype=2) + 
  labs(x = expression(Scale: log[2](Morgans)), 
       y = "Variance",
       color = "Year", shape = "Signal") +
  scale_x_discrete(breaks = c(paste0("d",1:13),"s13","chr"), labels = c(as.character(-14:-2),"-2 (scaling var)", "chromosome")) + 
  scale_shape_discrete(labels = c("Individual", "Mean"))+
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=13.05,y=-Inf,yend=-Inf),color="black")+
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=.4,margin=margin(t=-20)))


allVar[units == "genetic"] %>% ggplot(aes(x = level, y = propVar, group = interaction(signal, year), color = year)) +
  geom_point(aes(shape = signal), size=2.2) +
  geom_line(data = lineData[units == "genetic"], size=0.5, linetype=2) + 
  labs(x = expression(Scale: log[2](Morgans)), 
       y = "Proportion of genome-wide variance",
       color = "Year", shape = "Signal") +
  scale_x_discrete(breaks = c(paste0("d",1:13),"s13","chr"), labels = c(as.character(-14:-2),"-2 (scaling var)", "chromosome")) + 
  scale_shape_discrete(labels = c("Individual", "Mean"))+
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=13.05,y=-Inf,yend=-Inf),color="black")+
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=.4,margin=margin(t=-20)))

# physical units

allVar[, level := factor(level, levels = c(paste0("d",1:15), paste0("s",13:15), "chr"))]
lineData <- allVar[!level %in% c("s13", "s14", "s15", "chr")]
allVar[units == "physical"] %>% ggplot(aes(x = level, y = variance, group = interaction(signal, year), color = year)) +
  geom_point(aes(shape = signal), size=2.2) +
  geom_line(data = lineData[units == "physical"], size=0.5, linetype=2) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Variance",
       color = "Year", shape = "Signal") +
  scale_x_discrete(breaks = c(paste0("d",1:15),"s13","s14","s15","chr"), labels = c(as.character(0:14),"12(scaling)", "13(scaling)", "14(scaling)", "chromosome")) + 
  scale_shape_discrete(labels = c("Individual", "Mean"))+
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=.4,margin=margin(t=-20)))


allVar[units == "physical"] %>% ggplot(aes(x = level, y = propVar, group = interaction(signal, year), color = year)) +
  geom_point(aes(shape = signal), size=2.2) +
  geom_line(data = lineData[units == "physical"], size=0.5, linetype=2) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Proportion of total genomic variance",
       color = "Year", shape = "Signal") +
  scale_x_discrete(breaks = c(paste0("d",1:15),"s13","s14","s15","chr"), labels = c(as.character(0:14),"12(scaling)", "13(scaling)", "14(scaling)", "chromosome")) + 
  scale_shape_discrete(labels = c("Individual", "Mean"))+
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  theme(aspect.ratio = 1,
        text = element_text(size=15),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=.4,margin=margin(t=-20)))













































acua2006var <- data.table(var=loadFrom("ACUA_2006/totalVar.RData", "totalVar")); acua2006var[,year := "2006"]
acua2008var <- data.table(var=loadFrom("ACUA_2008/totalVar.RData", "totalVar")); acua2008var[,year := "2008"]
acua2013var <- data.table(var=loadFrom("ACUA_2013/totalVar.RData", "totalVar")); acua2013var[,year := "2013"]
acua2015var <- data.table(var=loadFrom("ACUA_2015/totalVar.RData", "totalVar")); acua2015var[,year := "2015"]
acua2018var <- data.table(var=loadFrom("ACUA_2018/totalVar.RData", "totalVar")); acua2018var[,year := "2018"]
allVar <- rbindlist(list(acua2006var,acua2008var,acua2013var,acua2015var,acua2018var))


ggplot(allVar, aes(x = year, y = var)) + 
  geom_point(size=2) + 
  theme_classic() + 
  labs(y="Total variance")


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
lineData <- allVar[scale < 17]

allVar[distance == "genetic"] %>% ggplot(aes(x = scale, y = anc_variance, group = interaction(decomp, year), color = year)) +
  geom_point(aes(shape = decomp), size=2.2) +
  geom_line(data = lineData[distance == "genetic"], size=0.5, linetype=2) + 
  labs(x = expression(Scale: log[2](Morgans)), 
       y = "Variance",
       color = "Year", shape = "Signal") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(-15:-1),"chromosome\nlevel")) + 
  scale_shape_discrete(labels = c("Individual", "Mean"))+
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black")+
  theme(aspect.ratio = 1,
    text = element_text(size=15),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5,size=12),
        axis.text.y = element_text(size=12),
        axis.ticks.x = element_line(size=c(rep(1,15),0)),
        axis.title.x = element_text(hjust=.4,margin=margin(t=-20)))

# 2.1. ----- Plot proportion ---------------
allVar[, propVar := anc_variance/sum(anc_variance), by = .(decomp,distance,year)]
lineDataProp <- allVar[scale < 17]
  
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


allVar[distance == "physical"] %>% ggplot(aes(x = scale, y = anc_variance, group = interaction(decomp, year), color = year)) +
  geom_point(aes(shape = rev(decomp)),size=2.2) +
  geom_line(data = lineDataProp[distance == "physical"], size=0.5,linetype=2) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Variance",
       color = "Year", shape = "Signal") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(0:14),"chrom")) + 
  scale_shape_discrete(labels = c("Sample mean", "Individual"))+
  theme_classic() +
  scale_colour_viridis_d() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black")+
  theme(aspect.ratio=1, text=element_text(size=15),
    axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5),#,
        axis.ticks.x = element_line(size=c(rep(1,15),0)))



chrLen[, scale:= floor(log2(V2/1000))]
ggplot(data=chrLen, aes(x = scale, y = 1)) + geom_boxplot()


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
  geom_point(size=2.2) + 
  geom_line(data = lineDataCor, size=0.5,linetype=2) + 
  labs(x = expression(Scale: log[2](kbp)), 
       y = "Wavelet correlation",
       color = "Year") +
  scale_x_continuous(breaks = c(1:15,17), labels = c(as.character(0:14),"chrom")) +  
  theme_classic() +
  geom_segment(aes(x=.95,xend=15.05,y=-Inf,yend=-Inf),color="black")+
  geom_segment(aes(x=16.5,xend=17.5,y=-Inf,yend=-Inf),color="black") +
  scale_color_viridis_d() +
  theme(aspect.ratio=1, text=element_text(size=15),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5),#,
        axis.ticks.x = element_line(size=c(rep(1,15),0)))


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
  geom_bar(aes(fill = reorder(as.factor(scale), desc(as.factor(scale))), x = year, y = contribution),
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


