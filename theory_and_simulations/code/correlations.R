source("datasets/wavelet_functions.R")
library(tidyverse)
library(data.table)
file1 <- "theory_and_simulations/results/equilibrium/ancestry_master.txt"
file2 <- "theory_and_simulations/results/add-sel-const-recomb/ancestry_master.txt"
file3 <- "theory_and_simulations/results/add-sel-periodic-recomb/ancestry_master.txt"

a1 <- read.table(file1, row.names = 1)
rownames(a1) <- paste0("neutral_", rownames(a1))

a2 <- read.table(file2, row.names = 1)
rownames(a2) <- paste0("sel-const-recomb_", rownames(a2))

a3 <- read.table(file3, row.names = 1)
rownames(a3) <- paste0("sel-periodic-recomb_", rownames(a3))


# combine data sets from simulations
all.sim <- rbind.data.frame(a1,a2,a3)

# reformat data for calculation and plotting
a <- as.data.frame(t(all.sim))
a$pos_absolute <- 1:1024

# create vector that describes recombination landscape (only applies to sims 3, 4)
x <- 1:1024
signal <- rep(0,1024)
for(i in 1:10){
  signal <- signal + sin(x*2*pi/2^i)
}
signal <- signal - min(signal)

const <- 1e-8/mean(exp(signal))
r <- const*exp(signal)

#plot(log(r), type = "n", ylab = "Crossover rate / bp", xlab = "Position", cex.lab = 1.5)
#lines(log(r))

a$recomb <- r

# add genetic distance corresponding to variable recombination rate (expected number of crossovers, binomial n*p)
a$pos_gen <- cumsum((1e8/1024)*r)

# tidy data
b  <-  a %>%
  gather(key = sim_rep_gen,
         value = freq, -c(pos_absolute,recomb,pos_gen)) %>% 
  separate(sim_rep_gen, c("sim_rep.id", "gen"), sep = "_gen") %>% 
  separate(sim_rep.id, c("sim", "rep.id"), sep = "_replicate")
setDT(b)

library(waveslim)

haar_modwt_coeffs <- function(x,variable,allcols){
  # can work for wavelet decomp on genetic and physical scale
  dt <- setDT(
    x[,brick.wall(modwt(get(variable), "haar", n.levels = floor(log2(length(get(variable))))), "haar")]
  )
  dt[, position := seq_len(.N)]
  dt <- melt(dt, id.vars = "position", variable.name = "level", value.name = "w")
  return(dt)
}

haar_modwt_var <- function(x,variable,allcols){
  # can work for wavelet decomp on genetic and physical scale
  dt <- setDT(
    x[,wave.variance(brick.wall(modwt(get(variable), "haar", n.levels = floor(log2(length(get(variable))))), "haar"))]
  )
  dt[, level := 1:floor(log(2))]
  dt[, position := seq_len(.N)]
  dt <- melt(dt, id.vars = "position", variable.name = "level", value.name = "w")
  return(dt)
}

# get wavelet coefficients

w1 <- b[,haar_dwt_coeffs(.SD, variable = "freq"), by = .(sim,rep.id,gen)]
w2 <- b[,haar_dwt_coeffs(.SD, variable = "recomb"), by = .(sim,rep.id,gen)]

W <- merge(w1,w2, by = c("sim", "rep.id", "gen", "level", "k"), suffixes = c("Freq", "Recomb"))

# calculate variances by scale

WV <- W[,  lapply(.SD, function(x){sum(x^2, na.rm=T)/length(x)}), .SDcols = c("wFreq","wRecomb"), by = .(sim,rep.id,gen,level)]
setnames(WV, c("wFreq", "wRecomb"), c("ancVar","recVar"))

# average proportion of variance
WV[, c("totVarFreq", "totVarRecomb")  := lapply(.SD, function(x)sum(x)), .SDcols = c("ancVar", "recVar"), by = .(sim,rep.id,gen)]
WV[, propVarFreq := ancVar/totVarFreq]
WV[, propVarRecomb := recVar/totVarRecomb]

# double check that variance calculations make sense
WV[, level := as.factor(level)]
WV[, level := factor(level, levels = 1:10)]
WV[sim == "neutral"] %>% ggplot(aes(x = level, y = ancVar, color = rep.id)) + facet_wrap(~gen) +
  geom_point() + geom_line(aes(group = rep.id)) 

# geometric mean of proportion of var
WV[, varWeight := sqrt(propVarFreq*propVarRecomb)]

# average over replicates
weights <- WV[, .( weight = mean(varWeight)), by = .(sim, gen, level)]

# calculate total correlations

cors <- b[sim == "sel-periodic-recomb", .(cor = cor(log10(recomb), freq)), by = .(rep.id, gen)]

cors[, z := 0.5*log((1+cor)/(1-cor))]
cors

avg.cors <- cors[, .(cor = (exp(2*mean(z))-1)/(exp(2*mean(z))+1)), by = gen]

plotA <- ggplot(cors, aes(x = log10(as.numeric(gen)), y = cor)) + 
  geom_line(aes(group = rep.id),  color = 'gray', alpha = 0.5) + geom_point(color = 'gray',alpha = 0.5) + 
  theme_classic() + 
  geom_line(data = avg.cors, size = 2) + geom_point(data = avg.cors, size = 2) + 
  labs(x = expression(log[10] (Generation)), 
       y = "Total correlation") + 
  theme(aspect.ratio=1,
        text=element_text(size=13))
 
# decompose correlation by level
rho <- W[sim == "sel-periodic-recomb", .(rho = sum(wFreq*wRecomb, na.rm=T)/sqrt(sum(wFreq^2, na.rm=T)*sum(wRecomb^2,na.rm=T))), by = .(gen, level, rep.id)]
rho2 <- W[sim == "sel-periodic-recomb", .(rho = cor(wRecomb, wFreq, use='pairwise.complete.obs')), by = .(gen, level, rep.id)]

rho[, z:= 0.5*log((1+rho)/(1-rho))]
rho2[, z:= 0.5*log((1+rho)/(1-rho))]

rho[rho > 0.999, z:= 5]
rho[rho < -0.999, z:=-5]
rho2[rho > 0.999, z:= 5]
rho2[rho < -0.999, z:=-5]

avg.rho <- rho[, .(rho = (exp(2*mean(z))-1)/(exp(2*mean(z))+1)), by = .(gen, level)]
avg.rho2 <- rho2[, .(rho = (exp(2*mean(z))-1)/(exp(2*mean(z))+1)), by = .(gen, level)]


avg.rho[, level := as.factor(level)]
avg.rho[, level := factor(level, levels = 1:10)]
avg.rho2[, level := as.factor(level)]
avg.rho2[, level := factor(level, levels = 1:10)]


library(grid)
txt <- textGrob("*", gp=gpar(fontsize=13, fontface="bold"))

plotB <- avg.rho %>%
  ggplot(aes(x = log10(as.numeric(gen)), y = rho, color = level))  +
  geom_point() + geom_line(aes(group = level)) + 
  scale_color_viridis_d(direction = -1, labels = -10:-1) + 
  theme_classic() +
  theme(aspect.ratio =1, 
        text = element_text(size = 13),
        plot.caption = element_text(size=20), 
        legend.position = "non") + 
  labs(color = "",
       x = expression(log[10] (Generation)), 
       y = expression(paste("\u03c1 ( ", w[x], " , ", w[y], ")"))) + annotation_logticks(sides = "b") + 
  annotation_custom(txt,xmin=log10(2),xmax=log10(2),ymin=-.9,ymax=-.9) + 
  coord_cartesian(clip = "off")
  #annotate(geom="text", x=log10(2), y=-.6, label="*", size = 15)
plotB

# stacked barplot
corPlotData <- merge(avg.rho, weights[sim == "sel-periodic-recomb"], by = c("gen","level"))
corPlotData[, contribution := rho*weight]

corPlotData[, sum.abs.contribution := sum(abs(contribution)), by = gen]
corPlotData[, abs.contribution := abs(contribution)/sum.abs.contribution]
corPlotData <- corPlotData[gen != "0001"]
corPlotData[, level := factor(level, levels = 10:1)]

plotC.leg <- corPlotData %>% 
  ggplot(aes(x = gen, fill = level, y = abs.contribution)) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() + 
  theme(aspect.ratio = 1,
        text = element_text(size = 13), 
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5),
        legend.title.align=0.5) + # ,
    #legend.direction = "horizontal") +
  scale_fill_viridis_d(labels = -1:-10) + 
  labs(fill = expression(Scale: log[2](Morgans)), 
       x = "Generation", y = "Contribution to total correlation") + 
  scale_x_discrete(labels = c(2,3,4,5,10,25,50,100,250,500,1000)) 
  #theme(legend.position="bottom",
  #    legend.spacing.x = unit(0, 'cm'))#+
plotC.leg

leg <- ggpubr::get_legend(plotC.leg)
#leg <- as_ggplot(leg)


plotC <- corPlotData %>% 
  ggplot(aes(x = gen, fill = level, y = abs.contribution)) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() + 
  theme(aspect.ratio = 1,
        text = element_text(size = 13), 
        legend.position = "",
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5)) + 
  scale_fill_viridis_d(labels = -1:-10) + 
  labs(fill = expression(Scale: log[2](Morgans)), 
       x = "Generation", y = "Contribution to total correlation") + 
  scale_x_discrete(labels = c(2,3,4,5,10,25,50,100,250,500,1000)) 


library(ggpubr)

toprow <- ggarrange(NULL, plotA, NULL, labels = c("","a",""),ncol = 3, widths = c(1,1.5,1))
toprow
bottomrow <- ggarrange(NULL, plotB, plotC, as_ggplot(leg), widths = c(0.5,1,1,0.5),labels = c("","b", "c",""), ncol = 4)
bottomrow
figure <- ggarrange(toprow, bottomrow, nrow = 2, labels = c("",""))
figure


