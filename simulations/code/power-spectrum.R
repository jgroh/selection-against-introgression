library(plyr)
library(dplyr)
library(tidyverse)
library(wavethresh)
library(stringr)
library(cubature)


###############################################################
## Infinite population approximation for expectation
###############################################################
M <- 20000
f <- function(x) {
  (1/M)*(
    a*exp(-r*t*(1e9/1024)*abs(x[2]-x[1])) + 
      (a^2)*(1-exp(-r*t*(1e9/1024)*abs(x[2]-x[1])))
  ) + 
    ((M-1)/M)*a^2
}

###############################################################
## General expectation formula
###############################################################
g <- function(x) {
  M <- 20000
  r <- 1e-8
  u <- M*r*(1e9/1024)*abs(x[2]-x[1])
  e <- exp(-(t/M)*(1+u))
  
  (1/M)*(
   a*(1+u*e)/(1+u) + a^2*(u*(1-e))/(1+u)
    )
   + 
    ((M-1)/M)*(
      a*(1-e)/(1+u) + a^2*(u+e)/(1+u)
    )
}


###############################################################
## Power spectrum of introgressed allele frequency through time
###############################################################
ps <- function(x){
  w <- wd(x$avg.frq, family = "DaubExPhase", filter.number = 1)
  x.var <- mean((x$avg.frq-mean(x$avg.frq))^2)
  temp <- vector(length = w$nlevels);
  for(i in 1:w$nlevels){
    temp[i] <- (sum((accessD(w,level=(w$nlevels - i)))^2)/length(x$avg.frq))/x.var
  }
  return(temp)
}

####################################################################
# compute correlation of wavelet coefficients at a single time point
####################################################################

wavelet_correlation <- function(data, yvar, xvar, fam = "DaubExPhase", filt = 1, detail = TRUE){
  
  x <- data[[xvar]]
  y <- data[[yvar]]
  total.cor <- cor(x, y, method = "pearson")
  
  # calculate empirical variance of each signal (n instead of n-1 in denominator)
  x.variance <- sum((x-mean(x))^2)/length(x) 
  y.variance <- sum((y-mean(y))^2)/length(y)
  
  # wavelet decomposition of each signal
  x.w <- wd(x, family = fam, filter.number = filt)
  y.w <- wd(y, family = fam, filter.number = filt)
  
  # fill lists while looping over scales of the transformation
  cor.list.detail <- list()
  cor.list.smooth <- list()
  alpha.x <- list()
  alpha.y <- list()
  weight <- list()
  p <- list()
  
  for(i in 1:(x.w$nlevels)){
    
    # extract wavelet coefficients at specified level
    
    # if(i < (x.w$nlevels+1)){
    # detail coefficients
    x.i.d <- accessD(x.w, level = (i-1))
    y.i.d <- accessD(y.w, level = (i-1))
    cor.list.detail[[i]] <- cor(x.i.d, y.i.d, method = "pearson")
    
    # smooth coefficients
    x.i.s <- accessC(x.w, level = (i-1))
    y.i.s <- accessC(y.w, level = (i-1))
    cor.list.smooth[[i]]  <- cor(x.i.s, y.i.s, method = "pearson")
      
    # calculate *proportion* of variance explained by each scale from detail coefficients
    alpha.x[[i]] <- (sum(x.i.d^2)/length(x)) /x.variance
    alpha.y[[i]] <- (sum(y.i.d^2)/length(y)) /y.variance
    
    weight[[i]] <- sqrt(alpha.x[[i]]*alpha.y[[i]])
    
    # calculate weight
    # fit linear model with intercept forced through origin, obtain r squared
    #s.lm <- summary(lm(y.i.d ~ x.i.d -1))
    #prp.expl[[i]] <- s.lm$adj.r.squared
    #p[[i]] <- cor(x.i.d, y.i.d)
    
    p[[i]] <- (x.i.d %*% y.i.d)/sqrt(sum(x.i.d^2)*sum(y.i.d^2))
    
    # } else{
    #   cor.list.detail[[i]] <- NA
    #   cor.list.smooth[[i]] <- NA
    #   alpha.x[[i]] <- x.variance
    #   alpha.y[[i]] <- y.variance
    #   s.lm <- summary(lm(y ~ x - 1))
    #   weight[[i]] <- s.lm$r.squared
    #   
    # }
    
    #alpha.x.s[[i]] <- (sum(x.i.s^2)/length(x))/x.variance
    #alpha.y.s[[i]] <- (sum(y.i.s^2)/length(y))/y.variance
    
  }
  
  
  names(cor.list.detail) <- paste0("detail-cor:scale",1:(x.w$nlevels))
  names(cor.list.smooth) <- paste0("smooth-cor:scale",1:(x.w$nlevels))
  
  names(alpha.x) <- paste0("alpha-", xvar, ":scale", 1:(x.w$nlevels))
  names(alpha.y) <- paste0("alpha-", yvar, ":scale", 1:(y.w$nlevels))
  
  names(weight) <- paste0("weight:scale", 1:(x.w$nlevels))
  names(p) <- paste0("p:scale", 1:(x.w$nlevels))
  return(c(unlist(cor.list.detail), unlist(cor.list.smooth), unlist(alpha.x), unlist(alpha.y), total.cor = total.cor, unlist(weight), unlist(p)))
}



####################################################################
# read and format input data
####################################################################

file1 <- "simulations/results/neutral-const-recomb/ancestry_master.txt"
file2 <- "simulations/results/neutral-periodic-recomb/ancestry_master.txt"
file3 <- "simulations/results/sel-const-recomb/ancestry_master.txt"
file4 <- "simulations/results/sel-periodic-recomb/ancestry_master.txt"


a1 <- read.table(file1, row.names = 1)
rownames(a1) <- paste0("neutral-const-recomb_", rownames(a1))

a2 <- read.table(file2, row.names = 1)
rownames(a2) <- paste0("neutral-periodic-recomb_", rownames(a2))

a3 <- read.table(file3, row.names = 1)
rownames(a3) <- paste0("sel-const-recomb_", rownames(a3))

a4 <- read.table(file4, row.names = 1)
rownames(a4) <- paste0("sel-periodic-recomb_", rownames(a4))


# combine data sets from simulations
all.sim <- rbind.data.frame(a1,a2,a3,a4)

# reformat data for calculation and plotting
a <- as.data.frame(t(all.sim))
a$pos_absolute <- 1:1024
# a$pos_on_chrom <- rep(1:64, 16)
# a$chrom <- rep(1:16, each=64)

# create vector that describes recombination landscape (only applies to sims 3, 4)
x <- 1:1024
signal <- rep(0,1024)
for(i in 1:10){
  signal <- signal + sin(x*2*pi/2^i)
}
signal <- signal - min(signal)

const <- 1e-8/mean(exp(signal))
r <- const*exp(signal)

plot(r, type = "n", ylab = "Crossover rate / bp", xlab = "Position", cex.lab = 1.5)
lines(r)

a$recomb <- r

b  <-  a %>%
  gather(key = sim_rep_gen,
         value = freq, -c(pos_absolute,recomb)) %>% 
  separate(sim_rep_gen, c("sim_rep.id", "gen"), sep = "_gen") %>% 
  separate(sim_rep.id, c("sim", "rep.id"), sep = "_replicate")


# calcualte average introgressed frequency per simulation per generation
d  <- b %>% 
   group_by(sim,pos_absolute,gen,recomb) %>% 
   dplyr::summarise(avg.frq = mean(freq))


########################################
# visualize allele frequency data
########################################


# look at frequency in a few snapshots to get a sense for the data
# b %>%
#   filter(sim == "neutral-no-recomb-var", gen == "0002") %>%
#   ggplot(aes(x = pos_absolute, y = freq)) +
#   geom_point() + facet_wrap(~rep.id)
#
# # There is an issue with how ancestry is calcualated - something about how fixed mutations are counted?
# doesn't make sense why frequency should be >.5 in all replicates in generation 2, and also to decay over time.


# plot average frequency over time, excluding neutral simulations for now
d %>% 
  #filter(sim == "neutral-const-recomb") %>% 
  filter(!gen %in% c( "0004", "50", "0250", "0500", "0750")) %>%
  filter(str_detect(sim, "sel-periodic")) %>% 
  ggplot(aes(x=pos_absolute, y=avg.frq)) +
  geom_point() + 
  facet_wrap(~gen, scales = "free_y") + xlab("Position") + ylab("Avg. introgressed ancestry proportion")

# examine effect of distance to center
d %>% 
  filter(sim == "sel-const-recomb") %>% 
  ggplot(aes(x=pos_absolute, y=avg.frq, group = sim, col = sim)) +
  geom_point() + 
  geom_smooth(method = "loess", color = "black") + 
  facet_wrap(~gen, scales = "free_y")





cor_text = d %>% 
  filter(str_detect(sim, "sel-periodic")) %>% 
  group_by(sim,gen) %>% 
  dplyr::summarise(correlation = round(cor.test(avg.frq,recomb)$estimate, 3) ) %>% #,
  #p.value = round(cor.test(freq,recomb)$p.value, 6)) %>% 
  mutate(text = paste0("r = ", correlation)) #, "\n", "P = ", p.value))


p <- d %>%
  filter(str_detect(sim, "sel-periodic")) %>% 
  filter(gen == "1000") %>%
  ggplot(aes(recomb, avg.frq)) +
  geom_point() +
  geom_smooth(method = "lm") +
  # scale_y_continuous(trans = "sqrt" ) +
  #scale_x_continuous(trans = "log10" ) +
  #facet_wrap(~gen, scales = "free_y") +
  theme(aspect.ratio=1, axis.text.x = element_text(angle = 90), axis.title = element_text(size=15)) + 
  ylab("Avg. introgressed ancestry proportion") + xlab("Crossover rate / bp")

p 
#p + geom_text(cor_text,
#              mapping = aes(x = Inf, y = Inf, label = text),
#              hjust=1,vjust=1)


########################################
# calculate and visualize power spectrum
########################################

# calculate power spectrum per simulation, generation
frq_ps <- ddply(d, .(sim,gen), ps)
colnames(frq_ps)[-c(1:2)] <- rev(round((1e9/2^(1:10))/1000000))  

# reformat 
frq_ps <- frq_ps %>% 
  gather(key = scale, value = ps, -c(sim,gen))

# order scale as factor for plotting
frq_ps$scale <- as.factor(frq_ps$scale)
frq_ps$scale <- factor(frq_ps$scale, levels = as.character(rev(round((1e9/2^(1:10))/1000000))))

# plot power spectrum
frq_ps %>% 
  filter(!gen %in% c("0001", "0002")) %>% 
  ggplot(aes(x = as.factor(scale), y = ps, group = sim, color = sim, pch = sim)) + 
  geom_line(position=position_dodge(w=0.1)) +
  geom_point() +
  facet_wrap(~gen) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Scale (Mb)")

#################################################
## plot expectation vs simulated data
#################################################
genlist <- list()
gen <- c("0003","0010","0050","0100","0500","1000")
for(i in 1:6){
  
  t <- as.numeric(gen[i])
  genvec_approx <- vector()
  genvec_exact <- vector()
  
  for(j in 1:10){
    
    part1 <- adaptIntegrate(f, lowerLimit = c(0,0), 
                            upperLimit = c(2^(j-1),2^(j-1)))
    part2 <- adaptIntegrate(f, lowerLimit = c(0,2^(j-1)),
                            upperLimit = c(2^(j-1),(2^j)))
    genvec_approx[j] <- ((part1$integral - part2$integral)/(2^(2*j-1)))
    
    part3 <- adaptIntegrate(g, lowerLimit = c(0,0), 
                            upperLimit = c(2^(j-1),2^(j-1)))
    part4 <- adaptIntegrate(g, lowerLimit = c(0,2^(j-1)),
                            upperLimit = c(2^(j-1),(2^j)))
    genvec_exact[j] <- ((part3$integral - part4$integral)/(2^(2*j-1)))
    
  }
  genvec_approx <- genvec_approx/sum(genvec_approx)
  genvec_exact <- genvec_exact/sum(genvec_exact)
  
  genlist[[i]] <- cbind(genvec_approx, genvec_exact)
}
names(genlist) <- gen

df <- do.call(cbind.data.frame, genlist)
df$scale <- as.character(rev( round(( (1e9/1024)*(1024/(2^(1:10))) ) /1e6)))

df2 <- df %>% gather(key = gen_formula,
                     value = variance, -c(scale)) %>% 
  separate(col = gen_formula, into = c("gen", "formula"), sep = ".genvec_")
df2

# order scale and generation as factors
df2$scale <- as.factor(df2$scale)
df2$scale <- factor(df2$scale, levels = as.character(rev( round(( (1e9/1024)*(1024/(2^(1:10))) ) /1e6))))
df2$gen <- as.factor(df2$gen)
df2$gen <- factor(df2$gen, levels = gen)
df2
pl <- ggplot(df2, aes(x = scale, y =variance, group = formula, color = formula)) + 
  geom_point() + geom_line() + facet_wrap(~gen) + ylab("Proportion of variance") + xlab("Scale (Mb)")
pl


sim.dat <- frq_ps %>% 
  filter(gen %in% c("0003","0010","0050","0100","0500","1000")) %>% 
  filter(sim == "neutral-const-recomb") %>% select(-sim) 
names(sim.dat)[names(sim.dat) == "ps"] <- "variance"
sim.dat$formula <- "simulation"

combined <- rbind(sim.dat, df2)

gen.labs <- c("Gen 3", "Gen 10", "Gen 50", "Gen 100", "Gen 500", "Gen 1000")
names(gen.labs) <- gen
combined %>%
  ggplot(aes(x = as.factor(scale), y = variance, group = formula, color = formula)) + 
  geom_line() + geom_point() + 
  facet_wrap(~gen, labeller = labeller(gen = gen.labs)) + 
  xlab("Scale (Mb)") + ylab("Proportion of variance") + 
  scale_color_manual(values = c("red", "blue", "black"), labels = c("\nApprox. for large N\n", "Exact", "\nSimulation\n(2N = 20000)")) + 
  theme(aspect.ratio = 1, axis.title = element_text(size = 14), 
        axis.text.x = element_text(angle = 90), legend.text = element_text(size= 14), 
        legend.key=element_blank(), legend.background = element_blank(),
        legend.title = element_blank()) 


val <- vector()
y <- (1e-8)*3*(1e9/1024)
for(j in 1:10){
  val[j] <- -4*0.5^2*(2^j/(2^(2*j)*y) + 
                          4*exp(-1/2*2^j*y)/(2^(2*j)*y^2) - 
                          exp(-2^j*y)/(2^(2*j)*y^2) - 3/(2^(2*j)*y^2)) + 
              4*0.5*(2^j/(2^(2*j)*y) + 
               4*exp(-1/2*2^j*y)/(2^(2*j)*y^2) - 
               exp(-2^j*y)/(2^(2*j)*y^2) - 3/(2^(2*j)*y^2))
}
val
val <- val/sum(val)
plot(val)


################################################################################################
# Wavelet correlation analysis
#################################################################################################

# compute correlations using average allele frequency
wav_cor <- d %>% 
  filter(str_detect(sim, "periodic-recomb"), !gen %in% c("0001","0002")) %>% 
  ddply( .(sim,gen), wavelet_correlation, xvar = "recomb", yvar = "avg.frq") %>% 
  select_if(~sum(!is.na(.)) > 0) %>% 
  gather(key = temp, value = temp_val, -c(gen,sim,total.cor)) %>% 
  separate(col = temp, into = c("thing", "scale"), sep = ":scale") %>% 
  spread(thing, temp_val) %>%
  mutate(weighted.cor = `weight`*`p`)

wav_cor <- wav_cor %>% 
  group_by(sim,gen) %>% 
  mutate(sum.weighted.cor = sum(weighted.cor), sum.abs.weighted.cor = sum(abs(weighted.cor)))


  
wav_cor$scale <- factor(wav_cor$scale, levels = as.character(1:10))


# plot correlation of detail coefficients 
wav_cor %>% 
  filter(sim == "sel-periodic-recomb") %>% #, !scale %in% c(1,2)) %>% 
  ggplot(aes(x=gen, y=`detail-cor`,group=scale,color=scale)) + 
  geom_point() + 
  geom_line() + 
  theme(aspect.ratio = 1) +
  facet_wrap(~sim) + ylab("Generation") 

# use the other correlation statistic from the decomposition formula
wav_cor %>% 
  filter(sim == "sel-periodic-recomb") %>% #, !scale %in% c(1,2)) %>% 
  ggplot(aes(x=gen, y=`p`,group=scale,color=scale)) + 
  geom_point() + 
  geom_line() + 
  theme(aspect.ratio = 1) +
  facet_wrap(~sim) + ylab("Detail coefficient 'correlation'") + xlab("Generation")



# plot correlation of smooth coefficients
wav_cor %>% 
  filter(sim == "sel-periodic-recomb", !scale %in% c(1,2)) %>% 
  ggplot(aes(x=gen, y=`smooth-cor`,group=scale,color=scale)) + 
  geom_point() + 
  geom_line() + 
  theme(aspect.ratio = 1) +
  facet_wrap(~sim)

# stacked barplot
wav_cor <- filter(wav_cor, sim == "sel-periodic-recomb")
ggplot(wav_cor) +
  geom_bar(aes(fill = scale, x = gen, y = weighted.cor), position = "stack", stat = "identity", color = "black") +
  scale_fill_viridis_d(option = "plasma") +
  geom_point(aes(x = gen, y = total.cor), size=4) +
  facet_wrap(~sim) + xlab("Generation") + ylab("Contribution to correlation")

  #geom_line(aes(x = gen, y = total.cor, group = 1)) +
  #geom_line(aes(x = gen,y = sum.weighted.cor, group = 1)) # this line confirms that the sum of weighted correlations equals the total correlation
  

# normalize by total correlation so that each generation sums to 1

wav_cor <- wav_cor %>% 
  group_by(gen) %>% 
  mutate(normalized.weighted.cor = abs(weighted.cor)/sum.abs.weighted.cor)


wav_cor %>% 
  ggplot() +
  geom_bar(aes(fill = scale, x = gen, y = normalized.weighted.cor), 
           position = "stack", stat = "identity", color = "black") +
  scale_fill_viridis_d(option = "plasma") +
  facet_wrap(~sim)
  
  







ggplot(wav_cor, aes(x = gen, y = total.cor, group = 1)) + geom_line() + facet_wrap(~sim)



# plot wavelet correlations 

detail_cor %>% 
  ggplot(aes(x=gen, y=correlation,group=scale,color=scale)) + 
  geom_point() + 
  geom_line() + 
  theme(aspect.ratio = 1)

smooth_cor %>% 
  ggplot(aes(x=gen, y=correlation,group=scale,color=scale)) + 
  geom_point() + 
  geom_line() + 
  theme(aspect.ratio = 1)


# compute correlation averaging over wavelet coefficients? ...


