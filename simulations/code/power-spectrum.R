library(plyr)
library(dplyr)
library(tidyverse)
library(wavethresh)
library(stringr)
library(cubature)


#trace(integral3, edit = T)

# Expected wavelet variance: approximate integrand ------------------------------------------
# assumes infinite population
wav_var_approx <- function(x, expected.crossovers.per.unit.dist, n.sample, alpha) {
  u <- expected.crossovers.per.unit.dist
  (1/n.sample)*(
    alpha*exp(-t*u*abs(x[2]-x[1])) + 
      (alpha^2)*(1-exp(-t*u*abs(x[2]-x[1])))
  ) + ((n.sample-1)/n.sample)*alpha^2
}

# Expected wavelet variance: exact  integrand ------------------------------------------

wav_var_exact <- function(x, expected.crossovers.per.unit.dist, n.pop, n.sample, alpha) {
  u <- expected.crossovers.per.unit.dist

  v <- n.pop*u*abs(x[2]-x[1])
  w <- exp(-(t/n.pop)*(1+v))
  (
  (1/n.sample)*(
   alpha*(1+v*w)/(1+v) + alpha^2*(v*(1-w))/(1+v)
    )
   + 
    ((n.sample-1)/n.sample)*(
      alpha*(1-w)/(1+v) + alpha^2*(v+w)/(1+v)) 
  )
}

# Continuous Haar wavelet function -----------------------------------------------
haarCts <- function(x, j){
  (x <= 0)*0 + (x > 0 & x <= 2^(j-1))*2^(-j/2) + (x > 2^(j-1) & x <= 2^j)*(-2^(-j/2))
}

# Expected wavelet variance: single sweep ----------------------------------------------
# assumes infinite population. x[1] and x[2] x[3] are positions of l1,l2,ls

wav_var_sweep <- function(x, j, r, n.sample, alpha, s, t) {
  
  #frequency of resident allele through time
  p <- (1-alpha)*exp(s*t/2) / ( alpha + (1-alpha)*exp(s*t/2) )
  q <- 1-p
  
  # recomb. probabilities
  f_prime <- function(a, b){
    exp(-r*abs(a-b) * 2*log(alpha + (1-alpha)*exp(s*t/2)) / s )
  }
  g_prime <- function(a,b){ 1 - f_prime(a,b) }
  f <- function(a,b){ exp(-t*r*abs(a-b)) }
  g <- function(a,b){1-f(a,b)}
  
  # integrand - nonzero only for l2 > l1, will multiply by two later
  if(x[2] > x[1]){
    if(x[3] <= x[1]){
      cov_ii <-  q*(f_prime(x[2],x[3]) + alpha*(g_prime(x[1],x[3])*f(x[1],x[2]) + 
                                                  f_prime(x[1],x[3])*g(x[1],x[2]) +
                                                  alpha*g_prime(x[1],x[3])*g(x[1],x[2]))) +
        p*(alpha*g_prime(x[1],x[3])*(f(x[1],x[2]) + alpha*g(x[1],x[2])))
    } else if(x[3] > x[1] && x[3] <= x[2]){
      cov_ii <-  q*(f_prime(x[1],x[3]) + alpha*g_prime(x[1],x[3]))*(f_prime(x[2],x[3]) + alpha*g_prime(x[2],x[3]))
      + p*alpha^2*g_prime(x[1],x[3])*g_prime(x[2],x[3])
    } else if(x[3] > x[2]){
      cov_ii <-  q*(f_prime(x[1],x[3]) + alpha*(g_prime(x[2],x[3])*f(x[1],x[2]) + 
                                                  f_prime(x[2],x[3])*g(x[1],x[2]) +
                                                  alpha*g_prime(x[2],x[3])*g(x[1],x[2]))) +
        p*(alpha*g_prime(x[2],x[3])*(f(x[1],x[2]) + alpha*g(x[1],x[2])))
    } else{cov_ii <- 0}
  } else{cov_ii <- 0}
  haarCts(x[1], j=j)*haarCts(x[2],j=j)*(1/n.sample)*cov_ii
}


# Power spectrum function ----------------------------------------------------------------
ps <- function(x, xvar){
  w <- wd(x[[xvar]], family = "DaubExPhase", filter.number = 1)
  temp <- vector(length = w$nlevels);
  x.var <- mean((x[[xvar]] - mean(x[[xvar]]))^2)
  for(i in 1:w$nlevels){
    temp[i] <- (sum((accessD(w,level=(w$nlevels - i)))^2)/length(x[[xvar]]))/x.var
  }
  return(temp)
}

# Wavelet correlation function ----------------------------------------------------------------

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
    p[[i]] <- (x.i.d %*% y.i.d)/sqrt(sum(x.i.d^2)*sum(y.i.d^2))
    
  }
  
  
  names(cor.list.detail) <- paste0("detail-cor:scale",1:(x.w$nlevels))
  names(cor.list.smooth) <- paste0("smooth-cor:scale",1:(x.w$nlevels))
  
  names(alpha.x) <- paste0("alpha-", xvar, ":scale", 1:(x.w$nlevels))
  names(alpha.y) <- paste0("alpha-", yvar, ":scale", 1:(y.w$nlevels))
  
  names(weight) <- paste0("weight:scale", 1:(x.w$nlevels))
  names(p) <- paste0("p:scale", 1:(x.w$nlevels))
  return(c(unlist(cor.list.detail), unlist(cor.list.smooth), unlist(alpha.x), unlist(alpha.y), total.cor = total.cor, unlist(weight), unlist(p)))
}

# Interpolate ancestry and recombination rate -------------------------------------------------------------------

getAncestry <- function(data, newPos){
  anc_pred <- vector() # will hold new values for ancestry
  rec <- vector() # will hold values for recombination 
  
  for(i in 1:length(newPos)){
  
    # next smallest value of genetic distance
    if (min(data$pos_gen) < newPos[i]){
      maxless <- max(data$pos_gen[data$pos_gen < newPos[i]])
      maxless.index <- which(data$pos_gen == maxless)[1]
    } else{
      maxless <- 0
    } 
    
    # next largest value of genetic distance
    if (max(data$pos_gen) == newPos[i]){
      mingreater <- max(data$pos_gen)
    } else{
      mingreater <- min(data$pos_gen[data$pos_gen > newPos[i]])
    } 
    mingreater.index <- which(data$pos_gen == mingreater)[1]
    
    # recombination measured by genetic distance between surrounding basepairs
    rec[i] <- abs(mingreater - maxless)
    
    # get two surrounding ancestry values for nearest bp's with known genetic distance
    # and then take weighted average weighting by proximity
    
    if ( mingreater == max(data$pos_gen) | maxless == 0 ){
      # if the new genetic position is below the minimum present in the data or is the maximum
      
      if( mingreater == max(data$pos_gen)){
        # if the new genetic position is the maximum, just use same ancestry value
        anc_pred[i] <- data$avg.frq[mingreater.index]
      }
      
      if( maxless == 0){
        # if the new genetic position is below the minimum, use the ancestry value
        # at the minimum genetic distance
        anc_pred[i] <- data$avg.frq[which(data$pos_gen == min(data$pos_gen))]
      }
      
    } else{
      # otherwise the new genetic position has a genetic position on either side with known ancestry value
      # compute weighted average weighted by distance
      
      d1 <- abs(maxless - newPos[i])
      d2 <- abs(mingreater - newPos[i])
      
      weight1 <- d1/(d1+d2)
      weight2 <- d2/(d1+d2)
      
      anc_pred[i] <- data$avg.frq[maxless.index]*weight1 + data$avg.frq[mingreater.index]*weight2
    }
  }
    ret <- paste(anc_pred, rec, sep = ":")
    names(ret) <- newPos
    return(ret)
}


# Read and format input data ----------------------------------------------------------------------------------------


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

# add genetic distance corresponding to variable recombination rate (expected number of crossovers, binomial n*p)
a$pos_gen <- cumsum((1e9/1024)*r)

# tidy data
b  <-  a %>%
  gather(key = sim_rep_gen,
         value = freq, -c(pos_absolute,recomb,pos_gen)) %>% 
  separate(sim_rep_gen, c("sim_rep.id", "gen"), sep = "_gen") %>% 
  separate(sim_rep.id, c("sim", "rep.id"), sep = "_replicate")

# average introgressed frequency per simulation,generation over replicates
d  <- b %>% 
   group_by(sim,pos_absolute,gen,recomb,pos_gen) %>% 
   dplyr::summarise(avg.frq = mean(freq))



# visualize allele frequency data: physical scale ----------------------------------------------------------------------------------------

# plot average frequency over time
d %>% 
  #filter(!gen %in% c( "0004", "50", "0250", "0500", "0750")) %>%
  #filter(str_detect(sim, "sel-periodic")) %>% 
  ggplot(aes(x=pos_absolute, y=avg.frq, group=sim,color=sim)) +
  geom_point() + 
  facet_wrap(~gen, scales = "free_y") + xlab("Position") + ylab("Avg. introgressed ancestry proportion")

# examine effect of distance to center
d %>% 
  filter(sim == "sel-const-recomb") %>% 
  ggplot(aes(x = pos_absolute, y = avg.frq)) +
  geom_point() + 
  geom_smooth(method = "loess") + 
  facet_wrap(~gen, scales = "free_y")

# calc correlation in different generations
cor_text = d %>% 
  filter(str_detect(sim, "sel-periodic")) %>% 
  group_by(sim,gen) %>% 
  dplyr::summarise(correlation = round(cor.test(avg.frq,recomb)$estimate, 3) ) %>% #,
  #p.value = round(cor.test(freq,recomb)$p.value, 6)) %>% 
  mutate(text = paste0("r = ", correlation)) #, "\n", "P = ", p.value))

# snapshot of correlation generation 1000
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


# visualize allele frequency data: genetic scale ----------------------------------------------------------------------------------------


# visualize relationship between physical position and genetic position
plot(a$pos_gen ~ a$pos_absolute, xlab = "physical position (~Mb)", ylab = "genetic position (M)")

plot(a$pos_absolute ~ a$pos_gen)
plot(rep(1, length(a$pos_gen)) ~ a$pos_gen)

# Highest density of data points for genetic distance where recombination is low
dens <- density(a$pos_gen, bw = "SJ")
plot(dens$y ~ dens$x, pch = 16)


# compare ancestry proportion along chrom using physical distance vs genetic distance
d %>% 
  filter(sim == "sel-periodic-recomb", gen == "0100") %>% 
  ggplot(aes(x=pos_absolute, y=avg.frq)) +
  geom_point() + 
  xlab("physical position (~1 Mb)") + 
  ylab("Average ancestry proportion")


d %>% 
  filter(sim == "sel-periodic-recomb", gen == "0100") %>% 
  ggplot(aes(x=pos_gen, y=avg.frq)) +
  geom_point() +
  xlab("genetic position (M)") + 
  ylab("Average ancestry proportion")


# compute ancestry and recombination values for evenly space intervals in genetic distance
# how to decide minimum distance between points?
newPos <- seq(min(d$pos_gen), max(d$pos_gen), length.out = length(a$pos_gen))

# what is the unit of distance between sites in Morgans?
unit <- newPos[2] - newPos[1]

gn <- d %>% 
  filter(str_detect(sim, "periodic-recomb")) %>% 
  ddply(.(sim,gen), getAncestry, newPos = newPos) %>% 
  gather(key = gen_pos_even, value = temp_val, -c(sim,gen), convert = TRUE) %>% 
  separate(col = temp_val, into = c("anc", "rec"), sep = ":", convert = TRUE) 

#examine interpolation: neutral case 
gn %>% 
  filter(sim == "neutral-periodic-recomb", gen == "1000") %>% 
  ggplot(aes(gen_pos_even, anc)) + 
  geom_point() + 
  labs(x = "Genetic position (M)",
       y = "Ancestry proportion (interpolated)",
       title = "Generation 1000")

# selection case
gn %>% 
  filter(sim == "sel-periodic-recomb", gen == "0100") %>% 
  ggplot(aes(gen_pos_even, anc)) + 
  geom_point() + 
  xlab("Genetic position (M)") + 
  ylab("Ancestry proportion (interpolated)")

# plot recombination rate on genetic scale
gn %>% 
  filter(gen == "0100") %>% 
  ggplot(aes(gen_pos_even,rec)) +
  geom_line() + 
  xlab("genetic position (M)") + 
  ylab("recombination rate (~ M/Mb)")


# calculate correlation between recombination and ancestry
gcor_text = gn %>%
  group_by(sim,gen) %>% 
  dplyr::summarise(correlation = round(cor.test(anc,rec)$estimate, 3) ) %>% 
  mutate(text = paste0("r = ", correlation))


h <- gn %>%
  filter(sim == "sel-periodic-recomb") %>% 
  ggplot(aes(rec,anc)) +
  geom_point() +
  facet_wrap(~gen, scales = "free_y") +
  geom_smooth(method = "lm") +
  # scale_y_continuous(trans = "sqrt" ) +
  #scale_x_continuous(trans = "log10" ) +
  theme(aspect.ratio=1) 

h + geom_text(filter(gcor_text, sim == "sel-periodic-recomb"),
              mapping = aes(x = Inf, y = Inf, label = text),
              hjust=1,vjust=1, col = "red")


# Power spectrum: physical scale  ----------------------------------------------------------------------------------------

# calculate power spectrum per simulation, generation
frq_ps <- ddply(d, .(sim,gen), ps, xvar = "avg.frq")
colnames(frq_ps)[-c(1:2)] <- rev(round((1e9/2^(0:9))/1000000))  

# reformat 
frq_ps <- frq_ps %>% 
  gather(key = scale, value = ps, -c(sim,gen))

# order scale as factor for plotting
frq_ps$scale <- as.factor(frq_ps$scale)
frq_ps$scale <- factor(frq_ps$scale, levels = as.character(rev(round((1e9/2^(1:10))/1000000))))

# plot power spectrum
frq_ps %>% 
  filter(! gen %in% c("0001", "0002")) %>% 
#  filter(str_detect(sim, "sel")) %>% 
  ggplot(aes(x = as.factor(scale), y = ps, group = sim, color = sim, pch = sim)) + 
  geom_line(position=position_dodge(w=0.1)) +
  geom_point() +
  facet_wrap(~gen) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Scale (Mb)")


# Power spectrum: genetic scale  ----------------------------------------------------------------------------------------

frq_ps_g <- ddply(gn, .(sim,gen), ps, xvar = "anc")
colnames(frq_ps_g)[-c(1:2)] <- rev(round( (10/2^(1:10))*100 )) # measuring in centimorgans
# reformat 
frq_ps_g <- frq_ps_g %>% 
  gather(key = scale, value = ps, -c(sim,gen))
# order scale as factor for plotting
frq_ps_g$scale <- factor(frq_ps_g$scale, levels = as.character(rev(round( (10/2^(1:10))*100 ))))

# plot power spectrum of interpolated ancestry 
frq_ps_g %>% 
  filter(! gen %in% c("0001", "0002")) %>% 
  #  filter(str_detect(sim, "sel")) %>% 
  ggplot(aes(x = as.factor(scale), y = ps, group = sim, color = sim, pch = sim)) + 
  geom_line(position=position_dodge(w=0.1)) +
  geom_point() +
  facet_wrap(~gen) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Scale (cM)")



# compute expected wavelet variance over a range of parameters: physical scale ---------------------------------------------------------------------------------

genlist <- list()
gen <- c("0003","0010","0050","0250","0500","1000")

# loop over generations 
for(i in 1:6){
  t <- as.numeric(gen[i])

  # loop over different population size, sample size, scale
  
  n.sample <- 2*c(0.5, 10, 100, 1000, 10000, 100000)
  n.pop <- 2*c(100, 1000, 10000, 100000, Inf)
  scale <- 1:10
  
  # make grid of parameters over which we evaluate the function
  grd <- expand.grid(n.sample=n.sample, n.pop=n.pop, scale=scale, stringsAsFactors = F)
  grd <- grd[grd$n.pop >= grd$n.sample,] # we only want evaluation where the sample is less than or equal to the population size
  
  grd$gen <- rep(t, nrow(grd)) # rep since we are inside the loop for a specific generation
  
  grd$variance <- vector(length = nrow(grd)) # this is the vector we fill in the calculation

  for(q in 1:nrow(grd)){
    j <- grd[q,]$scale
    n.sample <- grd[q,]$n.sample
    n.pop <- grd[q,]$n.pop

    if(n.pop == Inf){ # use infinite population approximation
      part1 <- adaptIntegrate(wav_var_approx, n.sample = n.sample, expected.crossovers.per.unit.dist=(1e9/1024)*1e-8, alpha=0.5, lowerLimit = c(0,0), 
                              upperLimit = c(2^(j-1),2^(j-1)))
      part2 <- adaptIntegrate(wav_var_approx, n.sample = n.sample, expected.crossovers.per.unit.dist=(1e9/1024)*1e-8, alpha=0.5, lowerLimit = c(0,2^(j-1)),
                              upperLimit = c(2^(j-1),(2^j)))
    } else { # use exact formula
      part1 <- adaptIntegrate(wav_var_exact, n.sample = n.sample, n.pop = n.pop, expected.crossovers.per.unit.dist=(1e9/1024)*1e-8, alpha=0.5, lowerLimit = c(0,0), 
                              upperLimit = c(2^(j-1),2^(j-1)))
      part2 <- adaptIntegrate(wav_var_exact, n.sample = n.sample, n.pop = n.pop, expected.crossovers.per.unit.dist=(1e9/1024)*1e-8, alpha=0.5, lowerLimit = c(0,2^(j-1)),
                              upperLimit = c(2^(j-1),(2^j)))
    }
    grd$variance[q] <- ((part1$integral - part2$integral)/(2^(2*j-1)))
  }
    genlist[[i]] <- grd
}

df <- do.call(rbind.data.frame, genlist)
#df$scale <- as.character(rev( round(( (1e9/1024)*(1024/(2^(1:10))) ) /1e6)))

#test_data <- filter(df, n.sample == 2000, n.pop == 2000, gen == 1000)
#sum(test_data$variance)

# compute proportion of variance by scale 
df <- df %>% group_by(n.sample,n.pop, gen) %>% 
  mutate(varsum=sum(variance)) %>% 
  mutate(prop.var = variance/varsum)


# Visualize effects of population and sample on expected wavelet variance -----------------------------------

# for plotting
df$scale <- factor(df$scale)
scale.labs <- as.character(rev( round(( (1e9/1024)*(1024/(2^(1:10))) ) /1e6)))

# Examine effect of sampling on variance. 
df %>% 
  filter(n.pop == Inf) %>% 
  ggplot(aes(x=scale,y=variance,group=as.factor(n.sample), color=as.factor(n.sample))) + 
  geom_point() + geom_line() + facet_wrap(~gen) + 
  scale_x_discrete(labels = scale.labs) + 
  labs(color = "Sample size", x = "Scale (Mb)", title = "Infinite population")

# Examine effect of sampling on *proportion* of variance. 
df %>% 
  filter(n.pop == Inf) %>% 
  ggplot(aes(x=scale,y=prop.var, group=n.sample, color=as.factor(n.sample))) + 
  scale_x_discrete(labels = scale.labs) + 
  geom_point() + geom_line() + facet_wrap(~gen) + 
  labs(color = "Sample size", x = "Scale (Mb)",
       y = "Proportion of variance", title = "Infinite Population")


# Examine effect of population size on variance. 
df %>% 
  filter(n.sample == 1) %>% 
  ggplot(aes(x=scale,y=variance, group=n.pop, color=as.factor(n.pop))) + 
  scale_x_discrete(labels = scale.labs) + 
  geom_point() + geom_line() + facet_wrap(~gen) + 
  labs(color = "Population", x = "Scale (Mb)",
       y = "Variance", title = "Sample size = 1")

# Examine effect of population size on *proportion* of variance. 
df %>% 
  filter(n.sample == 1) %>% 
  ggplot(aes(x=scale,y=prop.var, group=n.pop, color=as.factor(n.pop))) + 
  scale_x_discrete(labels = scale.labs) + 
  geom_point() + geom_line() + facet_wrap(~gen) + 
  labs(color = "Population", x = "Scale (Mb)",
       y = "Proportion of variance", title = "Sample size = 1")


# Effect of both population size and sample size on variance. 
 df %>% 
   filter(n.pop %in% c(200, Inf)) %>% 
   ggplot(aes(x=scale,y=variance, color = as.factor(n.sample), shape = as.factor(n.pop) , group = interaction(as.factor(n.pop), as.factor(n.sample)))) +   
   geom_point() + 
   geom_line(aes(linetype = as.factor(n.pop))) +
   facet_wrap(~gen) + labs(color = "Sample",
                           x = "Scale (Mb)",
                           shape = "Population", lty = "Population")
 
 # Effect of both population size and sample size on *proportion* of variance. 
 df %>% 
   filter(n.pop %in% c(200, Inf)) %>% 
   ggplot(aes(x=scale,y=prop.var, color = as.factor(n.sample), shape = as.factor(n.pop) , group = interaction(as.factor(n.pop), as.factor(n.sample)))) +   
   geom_point() + 
   geom_line(aes(linetype = as.factor(n.pop))) +
   facet_wrap(~gen) + labs(color = "Sample",
                           x = "Scale (Mb)",
                           y = "Proportion of variance",
                           shape = "Population", lty = "Population")


# Compute expected wavelet variance on genetic scale ---------------------------------------------
 
genlist <- list()
gen <- c("0003","0010","0050","1000")
 
# loop over generations 
for(i in 1:4){
   t <- as.numeric(gen[i])
   
   # loop over different population size, sample size, scale
   
   n.sample <- 2*c(0.5, 10, 100, 1000, 10000, 100000)
   n.pop <- 2*c(100, 1000, 10000, 100000, Inf)
   scale <- 1:10
   
   # make grid of parameters over which we evaluate the function
   grd <- expand.grid(n.sample=n.sample, n.pop=n.pop, scale=scale, stringsAsFactors = F)
   grd <- grd[grd$n.pop >= grd$n.sample,] # we only want evaluation where the sample is less than or equal to the population size
   
   grd$gen <- rep(t, nrow(grd)) # rep since we are inside the loop for a specific generation
   
   grd$variance <- vector(length = nrow(grd)) # this is the vector we fill in the calculation
   
   for(q in 1:nrow(grd)){
     j <- grd[q,]$scale
     n.sample <- grd[q,]$n.sample
     n.pop <- grd[q,]$n.pop
     
     if(n.pop == Inf){ # use infinite population approximation
       part1 <- adaptIntegrate(wav_var_approx, n.sample = n.sample, expected.crossovers.per.unit.dist=0.01, alpha=0.5, lowerLimit = c(0,0), 
                               upperLimit = c(2^(j-1),2^(j-1)))
       part2 <- adaptIntegrate(wav_var_approx, n.sample = n.sample, expected.crossovers.per.unit.dist=0.01, alpha=0.5, lowerLimit = c(0,2^(j-1)),
                               upperLimit = c(2^(j-1),(2^j)))
     } else { # use exact formula
       part1 <- adaptIntegrate(wav_var_exact, n.sample = n.sample, n.pop = n.pop, expected.crossovers.per.unit.dist=0.01, alpha=0.5, lowerLimit = c(0,0), 
                               upperLimit = c(2^(j-1),2^(j-1)))
       part2 <- adaptIntegrate(wav_var_exact, n.sample = n.sample, n.pop = n.pop, expected.crossovers.per.unit.dist=0.01, alpha=0.5, lowerLimit = c(0,2^(j-1)),
                               upperLimit = c(2^(j-1),(2^j)))
     }
     grd$variance[q] <- ((part1$integral - part2$integral)/(2^(2*j-1)))
   }
   genlist[[i]] <- grd
}
 
dfg <- do.call(rbind.data.frame, genlist)
 #df$scale <- as.character(rev( round(( (1e9/1024)*(1024/(2^(1:10))) ) /1e6)))
 
 #test_data <- filter(df, n.sample == 2000, n.pop == 2000, gen == 1000)
 #sum(test_data$variance)
 
# compute proportion of variance by scale 
dfg <- dfg %>% group_by(n.sample,n.pop, gen) %>% 
   mutate(varsum=sum(variance)) %>% 
   mutate(prop.var = variance/varsum)
 
 
# Visualize expected wavelet variance vs simulated data: genetic scale -----------------------------------

# for plotting
dfg$scale <- as.factor(dfg$scale)
levels(dfg$scale) <-  as.character(rev(round( (10/2^(1:10))*100 )))

# Examine effect of sampling on variance. 
dfg %>% 
  filter(n.pop == Inf) %>% 
  ggplot(aes(x=scale,y=variance,group=as.factor(n.sample), color=as.factor(n.sample))) + 
  geom_point() + geom_line() + facet_wrap(~gen) + 
  labs(color = "Sample size", x = "Scale (cM)", title = "Infinite population")
 

# reformat simulated data to combine with expectation
frq_ps_g$group <- "sim"
frq_ps_g$prop.var <- frq_ps_g$ps
dfg$group <- "expectation"

nm <- intersect(colnames(frq_ps_g), colnames(dfg))
sim.dat <- filter(frq_ps_g, sim == "neutral-periodic-recomb", gen %in% c("0003","0010","0050","1000"))[,nm]
exp.dat <-filter(dfg, n.sample == 20000, n.pop == 20000)[,nm] 
combined.g <- rbind.data.frame(sim.dat, exp.dat)
combined.g$gen <- as.numeric(combined.g$gen)

combined.g %>% 
  #filter(group == "expectation") %>% 
  ggplot(aes(x = scale, y = prop.var, group = group, color = group)) + 
  geom_point() +
  facet_wrap(~gen) + geom_line() + 
  labs(x = "Scale (cM)",
       y = "Proportion of variance",
       title = "Population = 20000, sample = 20000") +
  scale_color_manual(values = c("red", "black"), 
                     labels = c("Expectation",  "\nSimulated data \nw/ interpolation")) + 
  theme(aspect.ratio = 1, axis.title = element_text(size = 14), 
        axis.text.x = element_text(angle = 90), legend.text = element_text(size= 14), 
        legend.key=element_blank(), legend.background = element_blank(),
        legend.title = element_blank()) 




# Calculate expected wavelet variance with single sweep -------------------------

 
genlist <- list()
gen <- c("0005", "0010", "0050", "0100", "0500", "1000")
L <- 1024
# loop over generations 
for(i in 1:length(gen)){
  t <- as.numeric(gen[i])
  
  # loop over different population size, sample size, scale
  
  n.sample <- 2*c(0.5)
  scale <- 1:10
  
  # make grid of parameters over which we evaluate the function
  grd <- expand.grid(n.sample=n.sample, scale=scale, stringsAsFactors = F)

  grd$gen <- rep(t, nrow(grd)) # rep since we are inside the loop for a specific generation
  
  grd$variance <- vector(length = nrow(grd)) # this is the vector we fill in the calculation
  
  for(q in 1:nrow(grd)){
    j <- grd[q,]$scale
    n.sample <- grd[q,]$n.sample
    
    h <- hcubature(wav_var_sweep, c(0,0,0), c(3,3,3), j=j, n.sample=n.sample, alpha=0.5, s=0.01, t=t, r=0.01)
    
    grd$variance[q] <- h$integral/(1/2^(j-1))
    
  }
  genlist[[i]] <- grd
}


dfs <- do.call(rbind.data.frame, genlist)

# compute proportion of variance by scale 
dfs <- dfs %>% group_by(n.sample,gen) %>% 
  mutate(varsum=sum(variance)) %>% 
  mutate(prop.var = variance/varsum)

 
# for plotting
dfs$scale <- as.factor(dfs$scale)
levels(dfs$scale) <-  as.character(rev(round( (10/2^(1:10))*100 )))

# Examine effect of sampling on variance. 
dfs %>% 
  ggplot(aes(x=scale,y=variance,group=as.factor(n.sample), color=as.factor(n.sample))) + 
  geom_point() + geom_line() + facet_wrap(~gen) + 
  labs(color = "Sample size", x = "Scale (cM)", title = "Infinite population")


dfs %>% 
  #filter(group == "expectation") %>% 
  ggplot(aes(x = scale, y = variance, group = gen)) + 
  geom_point() +
  facet_wrap(~gen) + geom_line() + 
  labs(x = "Scale (cM)",
       y = "Proportion of variance",
       title = "Population = 20000, sample = 20000") +
  #scale_color_manual(values = c("red", "black"), 
                     #labels = c("Expectation",  "\nSimulated data \nw/ interpolation")) + 
  theme(aspect.ratio = 1, axis.title = element_text(size = 14), 
        axis.text.x = element_text(angle = 90), legend.text = element_text(size= 14), 
        legend.key=element_blank(), legend.background = element_blank(),
        legend.title = element_blank()) 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
# compare simulated data to expectation: pysical scale (need to fix) ----------------------------

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


# Power spectrum of ancestry on genetic scale

# calculate per simulation, generation
frq_ps_g <- ddply(g, .(sim,gen), ps, xvar = "anc")
colnames(frq_ps_g)[-c(1:2)] <- rev(round(max(a$pos_gen)/2^(0:9)*100))

# reformat 
frq_ps_g <- frq_ps_g %>% 
  gather(key = scale, value = ps, -c(sim,gen))

# order scale as factor for plotting
frq_ps_g$scale <- as.factor(frq_ps_g$scale)
frq_ps_g$scale <- factor(frq_ps_g$scale, levels = as.character(rev(round(max(a$pos_gen)/2^(0:9)*100))))

# plot power spectrum
frq_ps_g %>% 
  filter(! gen %in% c("0001", "0002")) %>% 
  #  filter(str_detect(sim, "sel")) %>% 
  ggplot(aes(x = as.factor(scale), y = ps, group = sim, color = sim, pch = sim)) + 
  geom_line(position=position_dodge(w=0.1)) +
  geom_point() +
  facet_wrap(~gen) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Scale (cM)")

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
  

##########################################
# Wavelet correlation on genetic scale
##########################################
  
# compute correlations using average allele frequency
wav_cor_g <- g %>% 
  filter(sim == "sel-periodic-recomb") %>% 
  ddply("gen", wavelet_correlation, xvar = "rec", yvar = "anc") %>% 
  select_if(~sum(!is.na(.)) > 0) %>% 
  gather(key = temp, value = temp_val, -c(gen,total.cor)) %>% 
  separate(col = temp, into = c("thing", "scale"), sep = ":scale") %>% 
  spread(thing, temp_val) %>%
  mutate(weighted.cor = `weight`*`p`)

wav_cor_g <- wav_cor_g %>% 
  group_by(gen) %>% 
  mutate(sum.weighted.cor = sum(weighted.cor), sum.abs.weighted.cor = sum(abs(weighted.cor)))



wav_cor_g$scale <- factor(wav_cor_g$scale, levels = as.character(1:10))


# plot correlation of detail coefficients 
wav_cor_g %>% 
  ggplot(aes(x=gen, y=`detail-cor`,group=scale,color=scale)) + 
  geom_point() + 
  geom_line() + 
  theme(aspect.ratio = 1) 

# plot correlation of smooth coefficients
wav_cor_g %>% 
  ggplot(aes(x=gen, y=`smooth-cor`,group=scale,color=scale)) + 
  geom_point() + 
  geom_line() + 
  theme(aspect.ratio = 1) 

ggplot(wav_cor_g) +
  geom_bar(aes(fill = scale, x = gen, y = weighted.cor), position = "stack", stat = "identity", color = "black") +
  scale_fill_viridis_d(option = "plasma") +
  geom_point(aes(x = gen, y = total.cor), size=4) 
#geom_line(aes(x = gen, y = total.cor, group = 1)) +
#geom_line(aes(x = gen,y = sum.weighted.cor, group = 1)) # this line confirms that the sum of weighted correlations equals the total correlation


# normalize by total correlation so that each generation sums to 1

wav_cor_g <- wav_cor_g %>% 
  group_by(gen) %>% 
  mutate(normalized.weighted.cor = abs(weighted.cor)/sum.abs.weighted.cor)

wav_cor_g %>% 
  ggplot() +
  geom_bar(aes(fill = scale, x = gen, y = normalized.weighted.cor), 
           position = "stack", stat = "identity", color = "black") +
  scale_fill_viridis_d(option = "plasma") 









