source("wavelet_functions.R")
source("~/workspace/gnomwav/R/multi_modwts.R")
source("~/workspace/gnomwav/R/variance_decomp.R")
source("~/workspace/gnomwav/R/correlation_decomp.R")

library(tidyverse)
library(data.table)


g1 <- "theory_and_simulations/results/add-sel-periodic-recomb_gen1/ancestry_master.txt"
g1_1000 <- "theory_and_simulations/results/add-sel-periodic-recomb_gen1-1000/ancestry_master.txt"
g100 <- "theory_and_simulations/results/add-sel-periodic-recomb_gen100/ancestry_master.txt"
g500 <- "theory_and_simulations/results/add-sel-periodic-recomb_gen500/ancestry_master.txt"

g1 <- read.table(g1, row.names = 1)
rownames(g1) <- paste0("g1_", rownames(g1))

g1_1000 <- read.table(g1_1000, row.names = 1)
rownames(g1_1000) <- paste0("g1_1000_", rownames(g1_1000))

g100 <- read.table(g100, row.names = 1)
rownames(g100) <- paste0("g100_", rownames(g100))

g500 <- read.table(g500, row.names = 1)
rownames(g500) <- paste0("g500_", rownames(g500))



# combine data sets from simulations
#all.sim <- rbind.data.frame(a1,a2,a3)

all.sim <- rbind.data.frame(g1,g100,g500,g1_1000)

# reformat data for calculation and plotting
a <- as.data.frame(t(all.sim))
a$pos_absolute <- 1:1000

# create vector that describes recombination landscape (only applies to sims 3, 4)
x <- 1:1000
signal <- rep(0,1000)
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
a$pos_gen <- cumsum((1e8/1000)*r)

# tidy data
b  <-  a %>%
  gather(key = sim_rep_gen,
         value = freq, -c(pos_absolute,recomb,pos_gen)) %>% 
  separate(sim_rep_gen, c("sim_rep.id", "gen"), sep = "_gen") %>% 
  separate(sim_rep.id, c("sim", "rep.id"), sep = "_replicate")
setDT(b)

#====== Calculate total correlations and covariances
# (and average over replicates)
totalcors <- b[, .(cor = cor(freq, log10(recomb))), by = .(sim, rep.id, gen)]
totalcovs <- b[, .(cov = cov(log10(recomb), freq)), by = .(rep.id, gen, sim)]

totalcors[, z := 0.5*log((1+cor)/(1-cor))]

avg.cors <- totalcors[, .(cor = (exp(2*mean(z))-1)/(exp(2*mean(z))+1)), by = .(gen, sim)]
avg.covs <- totalcovs[, .(cov = mean(cov)), by = .(gen, sim)]

# plot total correlations
ggplot(totalcors, aes(x = gen, y = cor)) + 
  geom_line(aes(group = rep.id),  color = 'gray', alpha = 0.5) + 
  geom_point(color = 'gray', alpha = 0.5) + 
  theme_classic() + 
  geom_line(data = avg.cors, aes(group = 1), size = 2) + 
  geom_point(data = avg.cors, size = 2) + 
  labs(x = expression(log[10] (Generation)), 
       y = "Total correlation") + 
  theme(aspect.ratio=1,
        text=element_text(size=13),
        axis.text.x = element_text(angle=90)) + facet_wrap(~sim) 


# ===== compute wavelet correlations

gcd <- b[, cor_tbl(data = .SD, 
                    chromosome = NA, 
                    signals = c("recomb", "freq"),
                    rm.boundary = TRUE), 
  by = .(rep.id, gen, sim)]

# plot all replicates
gcd[, gen := as.numeric(gen)]
ggplot(gcd, 
       aes(x = log10(gen), y = cor, group = interaction(level, rep.id), color = level)) +
  facet_wrap(~sim) + geom_point() + geom_line()
ggplot(gcd, 
       aes(x = log10(gen), y = cor, group = interaction(level, rep.id), color = level)) +
  facet_wrap(~sim) + geom_point() + geom_line()


# ---- average over replicates

# transform correlations
fisherz <- function(x){0.5*log((1+x)/(1-x))}
invfisherz <- function(x){ (exp(2*x)-1)/(exp(2*x)+1) }

gcd[, z := fisherz(cor)]
gcd[cor > 0.999, z:= 5]
gcd[cor < -0.999, z:=-5]

# average over transformed values, then reconvert
gcdm <- gcd[, .(cor = invfisherz(mean(z)) ), by = .(gen, sim, level)]

# compute confidence intervals on transformed scale then reconvert
se <- function(x){sd(x)/sqrt(length(x))}
gcdci <- gcd[, .(lower = mean(z) - 1.96*se(z), 
                 upper = mean(z) + 1.96*se(z)),
             by = .(gen,sim,level) ]
gcdci <- gcdci[, lapply(.SD, invfisherz),
               .SDcols = c("lower", "upper"), 
      by = .(gen, sim, level)]

# combine means and cis
gcdm <- merge(gcdm, gcdci, by = c("gen", "sim", "level"))

# plot
gcdm[, gen := as.numeric(gen)]

ggplot(gcdm, 
       aes(x = log10(gen), y = cor, group = level, color = level)) +
  geom_errorbar(aes(ymin=lower, ymax=upper))+
  facet_wrap(~sim) + geom_point() + geom_line()


# ===== Wavelet covariances ====

w <- b[, multi_modwts(.SD, chromosome = NA, signals = c("freq", "recomb"), rm.boundary = T),
  by = .(sim, rep.id, gen)]

meanFreq <- b[, .(p = mean(freq)), by = .(sim,gen)]
meanFreq[, gen := as.numeric(gen)]

covtbl <- w[, .(cov = cov(coefficient.freq, coefficient.recomb)), by = .(sim,rep.id,gen,level)]
covtbl[, gen := as.numeric(gen)]

covtblm <- covtbl[, .(cov = mean(cov)), by = .(sim,gen,level)]
covtblm <- merge(covtblm, meanFreq)

ggplot(covtblm, 
       aes(x = log10(gen), y = cov, group = level, color = level)) +
  #geom_errorbar(aes(ymin=lower, ymax=upper))+
  facet_wrap(~sim, scales = "free_y") + geom_point() + geom_line()

ggplot(covtblm, 
       aes(x = log10(gen), y = cov/sqrt(p*(1-p)), group = level, color = level)) +
  #geom_errorbar(aes(ymin=lower, ymax=upper))+
  facet_wrap(~sim) + geom_point() + geom_line()



# ===== Stacked Barplot of contribution to correlation =====

wv <- b[, gnom_var_decomp(.SD, chromosome = NA, signals = c("recomb", "freq")),
                      by = .(sim, rep.id, gen)]

# average variance across replicates
wv <- wv[, lapply(.SD, mean), .SDcols = c("variance.recomb", "variance.freq"),
   by = .(sim, gen, level)]
wv[, gen := as.numeric(gen)]
gcdm <- merge(gcdm, wv)

gcdm[, c("totalvar.recomb", "totalvar.freq") := lapply(.SD, sum), 
     .SDcols = c("variance.recomb", "variance.freq"), by = .(gen, sim)]

gcdm[, propvar.recomb := variance.recomb/totalvar.recomb]
gcdm[, propvar.freq := variance.freq/totalvar.freq]
gcdm[, contribution := cor*sqrt(propvar.freq*propvar.recomb)]


ggplot(gcdm) +
  geom_bar(aes(fill = level, x = as.factor(gen),
               y = contribution), position = "stack",
           stat = "identity", color = "black") +
  scale_fill_viridis_d(option = "plasma", direction = -1, 
                       labels = c(as.character(1:9), "9 (scaling)")) + labs(x = "Gen", 
       y = "Contribution",
       fill = "Scale: log[2](100 kb)") +
       #fill = expression(Scale: log[2](Morgan))) + 
  theme_classic() + facet_wrap(~sim, scales = "free_y") + 
  theme(aspect.ratio = 1, 
        axis.text.x = element_text(angle = 90))



#file1 <- "theory_and_simulations/results/equilibrium/ancestry_master.txt"
#file2 <- "theory_and_simulations/results/add-sel-const-recomb/ancestry_master.txt"
#file3 <- "theory_and_simulations/results/add-sel-periodic-recomb/ancestry_master.txt"

# a1 <- read.table(file1, row.names = 1)
# rownames(a1) <- paste0("neutral_", rownames(a1))
# 
# a2 <- read.table(file2, row.names = 1)
# rownames(a2) <- paste0("sel-const-recomb_", rownames(a2))
# 
# a3 <- read.table(file3, row.names = 1)
# rownames(a3) <- paste0("sel-periodic-recomb_", rownames(a3))











freq.modwt <- b[, brick.wall(wf= "haar", 
               modwt(freq, wf="haar", n.levels = 10)), 
  by = .(sim, rep.id, gen)] %>%
  melt(value.name = "freq.coefficient", 
       variable.name = "level",
       id.vars = c("sim", "rep.id", "gen"))
freq.modwt[, pos := seq_len(.N), by = .(sim, rep.id, gen)]

rec.modwt <- b[, brick.wall(wf= "haar", 
                             modwt(recomb, wf="haar", n.levels = 10)), 
                by = .(sim, rep.id, gen)] %>%
  melt(value.name = "rec.coefficient", 
       variable.name = "level",
       id.vars = c("sim", "rep.id", "gen"))
rec.modwt[, pos := seq_len(.N), by = .(sim, rep.id, gen)]

allmodwt <- merge(freq.modwt, rec.modwt, by = c("sim", "rep.id", "gen", "level", "pos"))

wv_freq <- b[, wave.variance(brick.wall(wf= "haar", 
               modwt(freq, wf="haar", n.levels = 10))), 
  by = .(sim, rep.id, gen)] 
wv_freq[, level := seq_len(.N), by = .(sim,rep.id,gen)]

wv_rec <- b[, wave.variance(brick.wall(wf= "haar", 
                                       modwt(recomb, wf="haar", n.levels = 10))), 
            by = .(sim, rep.id, gen)] 
wv_rec[, level := seq_len(.N), by = .(sim,rep.id,gen)]


allmodwt <- na.omit(allmodwt)


covtbl <- allmodwt[, .(cov = cov(freq.coefficient, rec.coefficient)), by = .(sim,rep.id,gen,level)]

covtbl[, level := NULL]
covtbl[, level := seq_len(.N), by = .(sim,rep.id,gen)]

allwt <- merge(merge(wv_rec,covtbl, by = c("sim", "rep.id", "gen", "level")), wv_freq)

allcor <- allwt[, .(cor = cov/sqrt(wavevar.x*wavevar.y)), by = .(sim,rep.id,gen, level)]
allcor

allcor[, z:= 0.5*log((1+cor)/(1-cor))]

allcor[cor > 0.999, z:= 5]
allcor[cor < -0.999, z:=-5]

avg.cors <- allcor[level %in% 1:9, .(avg.cor = (exp(2*mean(z))-1)/(exp(2*mean(z))+1)), by = .(sim, gen, level)]
avg.cors[, level := as.factor(level)]

ggplot(avg.cors, aes(x= gen, y = avg.cor, group = level, color = level)) + 
  scale_color_viridis_d(direction=-1) +
  geom_point() + geom_line() + facet_wrap(~sim)

















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

w1 <- b[, wave.variance(brick.wall(modwt(freq, n.levels=10, wf="haar"),"haar")), by = .(sim,rep.id,gen)]
w2 <- b[, wave.variance(brick.wall(modwt(recomb, n.levels=10, wf="haar"),"haar")), by = .(sim,rep.id,gen)]
w1[, level := 1:11, by = .(sim,rep.id,gen)]
w2[, level := 1:11, by = .(sim,rep.id,gen)]

W <- merge(w1,w2, by = c("sim", "rep.id", "gen", "level", "k"), suffixes = c("Freq", "Recomb"))

W <- merge(w1,w2, by = c("sim", "rep.id", "gen", "level"), suffixes = c("Freq", "Recomb"))


# calculate variances by scale

WV <- W[,  lapply(.SD, function(x){sum(x^2, na.rm=T)/length(x)}), .SDcols = c("wFreq","wRecomb"), by = .(sim,rep.id,gen,level)]
WV <- W[,  lapply(.SD, function(x){mean(x^2, na.rm=T)/length(x)}), .SDcols = c("wavevarFreq","wavevarRecomb"), by = .(sim,rep.id,gen,level)]

setnames(WV, c("wavevarFreq", "wavevarRecomb"), c("ancVar","recVar"))

# average proportion of variance
WV[, c("totVarFreq", "totVarRecomb")  := lapply(.SD, function(x)sum(x)), .SDcols = c("ancVar", "recVar"), by = .(sim,rep.id,gen)]
WV[, propVarFreq := ancVar/totVarFreq]
WV[, propVarRecomb := recVar/totVarRecomb]

# double check that variance calculations make sense
WV[, level := as.factor(level)]
WV[, level := factor(level, levels = 1:10)]
WV[sim == "g1" & level != 11] %>% ggplot(aes(x = level, y = ancVar, color = rep.id)) + 
  facet_wrap(~gen, scales = "free_y") +
  geom_point() + geom_line(aes(group = rep.id)) 

# geometric mean of proportion of var
WV[, varWeight := sqrt(propVarFreq*propVarRecomb)]

# average over replicates
weights <- WV[, .( weight = mean(varWeight)), by = .(sim, gen, level)]



# calculate total correlations

cors <- b[, .(cor = cor(log10(recomb), freq)), by = .(rep.id, gen, sim)]
vars <- b[, lapply(.SD, var), .SDcols = c("freq", "recomb"), by = .(rep.id, gen, sim)]

avg.vars <- vars[, lapply(.SD, mean), .SDcols = c("freq"), by = .(gen,sim)]

cors[, z := 0.5*log((1+cor)/(1-cor))]
cors

covs <- b[, .(cov = cov(log10(recomb), freq)), by = .(rep.id, gen, sim)]

avg.cors <- cors[, .(cor = (exp(2*mean(z))-1)/(exp(2*mean(z))+1)), by = .(gen,sim)]
avg.covs <- covs[, .(cov = mean(cov)), by = .(gen,sim)]



cors[, sim := factor(sim, levels = c("g1", "g2", "g1.5","g5.10", "g10.50", "g50.100"))]
covs[, sim := factor(sim, levels = c("g1", "g2", "g1.5","g5.10", "g10.50", "g50.100"))]

corcomponents <- merge(avg.covs, avg.vars)
corcomponents <- melt(corcomponents, measure.vars = c("cov", "freq"), id.vars = c("gen", "sim"))

nms <- list(
  "g1"="1",
  "g2"="2",
  "g1.5"="1-5",
  "g5.10"="5-10",
  "g10.50"="10-50",
  "g50.100"="50-100")
lblr <- function(variable, value){
  return(nms[value])
}

# plot cor components
ggplot(corcomponents, aes(x = log10(as.numeric(gen)), 
                          y = value, group = variable, color = variable)) + 
         geom_point() + geom_line() + facet_wrap(~sim) +
  theme_classic()+
  scale_colour_brewer(type = "qual", palette = "Set2",
                      labels = c("Covariance", "Ancestry variance")) +
  labs(x = expression(log[10] (Generation)), 
       y = "Value") +
  theme(aspect.ratio=1,
        text=element_text(size=13)) + facet_wrap(~sim, labeller=lblr) 


plotA <- ggplot(cors, aes(x = log10(as.numeric(gen)), y = cor)) + 
  geom_line(aes(group = rep.id),  color = 'gray', alpha = 0.5) + geom_point(color = 'gray',alpha = 0.5) + 
  theme_classic() + 
  geom_line(data = avg.cors, size = 2) + geom_point(data = avg.cors, size = 2) + 
  labs(x = expression(log[10] (Generation)), 
       y = "Total correlation") + 
  theme(aspect.ratio=1,
        text=element_text(size=13)) + facet_wrap(~sim, labeller=lblr) 
plotA
 
# decompose correlation by level
rho <- W[, .(rho = sum(wFreq*wRecomb, na.rm=T)/sqrt(sum(wFreq^2, na.rm=T)*sum(wRecomb^2,na.rm=T))), by = .(gen, level, sim, rep.id)]
rho <- allmodwt[, .(rho = cov(freq.coefficient,rec.coefficient)/sqrt(var(freq.coefficient)*var(rec.coefficient))), by = .(gen, level, sim, rep.id)]

#rho <- allmodwt[, .(rho = cov(freq.coefficient,rec.coefficient)), by = .(gen, level, sim, rep.id)]


rho2 <- W[, .(rho = cor(wRecomb, wFreq, use='pairwise.complete.obs')), by = .(gen, level, rep.id,sim)]
rho2 <- allmodwt[, .(rho = cor(freq.coefficient, rec.coefficient, use='pairwise.complete.obs')), by = .(gen, level, rep.id,sim)]



rho[, z:= 0.5*log((1+rho)/(1-rho))]
rho2[, z:= 0.5*log((1+rho)/(1-rho))]

rho[rho > 0.999, z:= 5]
rho[rho < -0.999, z:=-5]
rho2[rho > 0.999, z:= 5]
rho2[rho < -0.999, z:=-5]

avg.rho <- rho[, .(rho = (exp(2*mean(z))-1)/(exp(2*mean(z))+1)), by = .(gen, sim,level)]
avg.rho2 <- rho2[, .(rho = (exp(2*mean(z))-1)/(exp(2*mean(z))+1)), by = .(gen, sim, level)]

# for just cov
#avg.rho <- rho[, .(rho = mean(rho)), by = .(gen, sim,level)]

avg.rho[, level := as.factor(level)]
avg.rho[, level := factor(level, levels = 1:10)]
avg.rho2[, level := as.factor(level)]
avg.rho2[, level := factor(level, levels = 1:10)]


library(grid)
txt <- textGrob("*", gp=gpar(fontsize=13, fontface="bold"))

avg.rho
plotB <- avg.rho[level %in% paste0("d", 1:9)] %>%
  ggplot(aes(x = log10(as.numeric(gen)), y = rho, color = level))  +
  geom_point() + geom_line(aes(group = level)) + 
  scale_color_viridis_d(direction = -1, labels = -10:-1) + 
  theme_classic() +
  theme(aspect.ratio =1, 
        text = element_text(size = 13),
        plot.caption = element_text(size=20))+ #, 
        #legend.position = "non") + 
  labs(color = expression(Scale: log[2](Morgans)),
       x = expression(log[10] (Generation)), 
       y = expression(paste("\u03c1 ( ", w[x], " , ", w[y], ")"))) + annotation_logticks(sides = "b") + 
  annotation_custom(txt,xmin=log10(2),xmax=log10(2),ymin=-.9,ymax=-.9) + 
  coord_cartesian(clip = "off") + facet_wrap(~sim, labeller = lblr) + #, scales = "free_y") +
  scale_fill_viridis_d(labels = -1:-10) 
  #annotate(geom="text", x=log10(2), y=-.6, label="*", size = 15)
plotB

# stacked barplot
avg.rho[, level := NULL]
avg.rho[, level := 1:11, by = .(gen,sim)][, level := as.factor(level)][]
corPlotData <- merge(avg.rho[level %in% 1:9], weights, by = c("sim", "gen","level"))
corPlotData[, contribution := rho*weight]

corPlotData[, sum.abs.contribution := sum(abs(contribution)), by = .(gen,sim)]
corPlotData[, abs.contribution := abs(contribution)/sum.abs.contribution]

corPlotData <- corPlotData[gen != "0001"]
corPlotData[, level := factor(level, levels = 10:1)]

plotC.leg <- corPlotData %>% 
  ggplot(aes(x = gen, fill = level, y = abs.contribution)) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() + 
  facet_wrap(~sim) +
  theme(aspect.ratio = 1,
        text = element_text(size = 13), 
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5),
        legend.title.align=0.5) + # ,
    #legend.direction = "horizontal") +
  scale_fill_viridis_d(labels = -1:-10) + 
  labs(fill = expression(Scale: log[2](Morgans)), 
       x = "Generation", y = "Contribution to total correlation")  
  #scale_x_discrete(labels = c(2,3,4,5,10,25,50,100,250,500,1000)) 
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


