library(plyr)
library(dplyr)
library(tidyverse)
library(wavethresh)
library(stringr)


###############################################################
## Power spectrum of introgressed allele frequency through time
###############################################################
ps <- function(x, xvar){
  w <- wd(x[[xvar]], family = "DaubExPhase", filter.number = 1)
  temp <- vector(length = w$nlevels);
  for(i in 1:w$nlevels){
    temp[i] <- sum((accessD(w,level=(w$nlevels - i)))^2);
  }
  temp <- temp/sum(temp);
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
  #alpha.x.s <- list()
  #alpha.y.s <- list()
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
    s.lm <- summary(lm(y.i.d ~ x.i.d -1))
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
# function to get value of ancestry at evenly spaced intervals in genetic distance
# and calculate recombination rate
####################################################################
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
    
    # recombination given by genetic distance between surrounding basepairs
    rec[i] <- abs(mingreater - maxless)
    
    # get two surrounding ancestry values for nearest bp's with known genetic distance
    # and then take weighted average weighting by proximity
    
    if ( mingreater == max(data$pos_gen) | maxless == 0 ){
      # if the new genetic position is below the minimum or is the maximum
      
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

plot(r, type = "n", ylab = "Crossover rate / bp", xlab = "Position")
lines(r)

a$recomb <- r

# add genetic distance
a$pos_gen <- cumsum((1e9/1024)*r)


b  <-  a %>%
  gather(key = sim_rep_gen,
         value = freq, -c(pos_absolute,recomb,pos_gen)) %>% 
  separate(sim_rep_gen, c("sim_rep.id", "gen"), sep = "_gen") %>% 
  separate(sim_rep.id, c("sim", "rep.id"), sep = "_replicate")



# calcualte average introgressed frequency per simulation per generation
d  <- b %>% 
   group_by(sim,pos_absolute,gen,recomb,pos_gen) %>% 
   dplyr::summarise(avg.frq = mean(freq))


########################################
# visualize allele frequency data
########################################


# plot average frequency over time
d %>% 
  #filter(str_detect(sim, "sel")) %>% 
  ggplot(aes(x=pos_absolute, y=avg.frq, group = sim, col = sim)) +
  geom_point() + 
  facet_wrap(~gen, scales = "free_y")


# calculate correlation between recombination and ancestry
cor_text = d %>% 
  filter(str_detect(sim, "sel-periodic")) %>% 
  group_by(sim,gen) %>% 
  summarise(correlation = round(cor.test(avg.frq,recomb)$estimate, 3) ) %>% #,
  #p.value = round(cor.test(freq,recomb)$p.value, 6)) %>% 
  mutate(text = paste0("r = ", correlation)) #, "\n", "P = ", p.value))


p <- d %>%
  filter(str_detect(sim, "sel-periodic")) %>% 
  ggplot(aes(recomb, avg.frq)) +
  geom_point() +
  geom_smooth(method = "lm") +
  # scale_y_continuous(trans = "sqrt" ) +
  #scale_x_continuous(trans = "log10" ) +
  facet_wrap(~gen, scales = "free_y") +
  theme(aspect.ratio=1) 

p + geom_text(cor_text,
              mapping = aes(x = Inf, y = Inf, label = text),
              hjust=1,vjust=1)


########################################
# plot on genetic scale
########################################

# visualize relationship between physical position and genetic position
plot(a$pos_gen ~ a$pos_absolute, xlab = "physical position (~Mb)", ylab = "genetic position (M)")

plot(a$pos_absolute ~ a$pos_gen)
plot(rep(1, length(a$pos_gen)) ~ a$pos_gen)

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
newPos <- seq(min(d$pos_gen), max(d$pos_gen), length.out = length(a$pos_gen))
g <- d %>% 
  filter(str_detect(sim, "periodic-recomb")) %>% 
  ddply(.(sim,gen), getAncestry, newPos = newPos) %>% 
  gather(key = gen_pos_even, value = temp_val, -c(sim,gen), convert = TRUE) %>% 
  separate(col = temp_val, into = c("anc", "rec"), sep = ":", convert = TRUE) 

g %>% 
  filter(sim == "sel-periodic-recomb", gen == "0100") %>% 
  ggplot(aes(gen_pos_even, anc)) + 
  geom_point() + 
  xlab("Genetic position (M)") + 
  ylab("Ancestry proportion (interpolated)")

# plot recombination rate on genetic scale
g %>% 
  filter(gen == "0100") %>% 
  ggplot(aes(gen_pos_even,rec)) +
  geom_line() + 
  xlab("genetic position (M)") + 
  ylab("recombination rate (~ M/Mb)")


# calculate correlation between recombination and ancestry
gcor_text = g %>%
  group_by(sim,gen) %>% 
  dplyr::summarise(correlation = round(cor.test(anc,rec)$estimate, 3) ) %>% 
  mutate(text = paste0("r = ", correlation))


h <- g %>%
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


########################################
# calculate and visualize power spectrum
########################################

# calculate power spectrum per simulation, generation
frq_ps <- ddply(d, .(sim,gen), ps, xvar = "avg.frq")
colnames(frq_ps)[-c(1:2)] <- rev(round((1e9/2^(0:9))/1000000))  

# reformat 
frq_ps <- frq_ps %>% 
  gather(key = scale, value = ps, -c(sim,gen))

# order scale as factor for plotting
frq_ps$scale <- as.factor(frq_ps$scale)
frq_ps$scale <- factor(frq_ps$scale, levels = as.character(rev(round((1e9/2^(0:9))/1000000))))

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


#################################################################################################
# Wavelet correlation analysis
#################################################################################################

# compute correlations using average allele frequency
wav_cor <- d %>% 
  filter(sim == "sel-periodic-recomb") %>% 
  ddply("gen", wavelet_correlation, xvar = "recomb", yvar = "avg.frq") %>% 
  select_if(~sum(!is.na(.)) > 0) %>% 
  gather(key = temp, value = temp_val, -c(gen,total.cor)) %>% 
  separate(col = temp, into = c("thing", "scale"), sep = ":scale") %>% 
  spread(thing, temp_val) %>%
  mutate(weighted.cor = `weight`*`p`)

wav_cor <- wav_cor %>% 
  group_by(gen) %>% 
  mutate(sum.weighted.cor = sum(weighted.cor), sum.abs.weighted.cor = sum(abs(weighted.cor)))


  
wav_cor$scale <- factor(wav_cor$scale, levels = as.character(1:10))


# plot correlation of detail coefficients 
wav_cor %>% 
  ggplot(aes(x=gen, y=`detail-cor`,group=scale,color=scale)) + 
  geom_point() + 
  geom_line() + 
  theme(aspect.ratio = 1) 

# plot correlation of smooth coefficients
wav_cor %>% 
  ggplot(aes(x=gen, y=`smooth-cor`,group=scale,color=scale)) + 
  geom_point() + 
  geom_line() + 
  theme(aspect.ratio = 1) 

ggplot(wav_cor) +
  geom_bar(aes(fill = scale, x = gen, y = weighted.cor), position = "stack", stat = "identity", color = "black") +
  scale_fill_viridis_d(option = "plasma") +
  geom_point(aes(x = gen, y = total.cor), size=4) 
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
  scale_fill_viridis_d(option = "plasma") 
  

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









