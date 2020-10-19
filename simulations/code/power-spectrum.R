library(plyr)
library(dplyr)
library(tidyverse)
library(wavethresh)
library(stringr)


###############################################################
## Power spectrum of introgressed allele frequency through time
###############################################################
ps <- function(x){
  w <- wd(x$avg.frq, family = "DaubExPhase", filter.number = 1)
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

wavelet_correlation <- function(data, var1, var2, fam = "DaubExPhase", filt = 1, detail = TRUE){
  
  total.cor <- cor(data[[var1]], data[[var2]], method = "pearson")
  x.variance <- var(data[[var1]])
  y.variance <- var(data[[var2]])
  
  # designed to work for a given time point 
  x.w <- wd(data[[var1]]-mean(data[[var1]]), family = fam, filter.number = filt)
  y.w <- wd(data[[var2]]-mean(data[[var2]]), family = fam, filter.number = filt)
  
  cor.list.detail <- list()
  cor.list.smooth <- list()
  alpha.x <- list()
  alpha.y <- list()
  
  for(i in 1:x.w$nlevels){
    
    # extract wavelet coefficients at specified level
    
    # detail coefficients
    x.i.d <- accessD(x.w, level = (i-1))
    y.i.d <- accessD(y.w, level = (i-1))
    cor.list.detail[[i]] <- cor(x.i.d, y.i.d, method = "pearson")
    
    # smooth coefficients
    x.i.s <- accessC(x.w, level = (i-1))
    y.i.s <- accessC(y.w, level = (i-1))
    cor.list.smooth[[i]]  <- cor(x.i.s, y.i.s, method = "pearson")
      
    # calculate proportion of variance explained by each scale from detail coefficients
    alpha.x[[i]] <- sum(x.i.d^2)
    alpha.y[[i]] <- sum(y.i.d^2)
    
  }
  alpha.x <- unlist(alpha.x)/sum(unlist(alpha.x))
  alpha.y <- unlist(alpha.y)/sum(unlist(alpha.y))
  
  names(cor.list.detail) <- paste0("detail-cor:scale",1:x.w$nlevels)
  names(cor.list.smooth) <- paste0("smooth-cor:scale",1:x.w$nlevels)
  
  names(alpha.x) <- paste0("alpha-", var1, ":scale", 1:x.w$nlevels)
  names(alpha.y) <- paste0("alpha-", var2, ":scale", 1:y.w$nlevels)
  return(c(unlist(cor.list.detail), unlist(cor.list.smooth), alpha.x, alpha.y,total.cor = total.cor))
}



####################################################################
# read and format input data
####################################################################

file1 <- "simulations/results/neutral-no-recomb-var/ancestry_master.txt"
file2 <- "simulations/results/sel-no-recomb-var/ancestry_master.txt"
file3 <- "simulations/results/neutral-recomb-var/ancestry_master.txt"
file4 <- "simulations/results/sel-recomb-var/ancestry_master.txt"


a1 <- read.table(file1, row.names = 1)
rownames(a1) <- paste0("neutral-no-recomb-var_", rownames(a1))

a2 <- read.table(file2, row.names = 1)
rownames(a2) <- paste0("sel-no-recomb-var_", rownames(a2))

a3 <- read.table(file3, row.names = 1)
rownames(a3) <- paste0("neutral-recomb-var_", rownames(a3))

a4 <- read.table(file4, row.names = 1)
rownames(a4) <- paste0("sel-recomb-var_", rownames(a4))


# combine data sets from simulations
all.sim <- rbind.data.frame(a2,a4)

# reformat data for calculation and plotting
a <- as.data.frame(t(all.sim))
a$pos_absolute <- 1:1024
a$pos_on_chrom <- rep(1:64, 16)
a$chrom <- rep(1:16, each=64)

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

b  <-  a %>%
  gather(key = sim_rep_gen,
         value = freq, -c(pos_absolute,pos_on_chrom,chrom,recomb)) %>% 
  separate(sim_rep_gen, c("sim_rep.id", "gen"), sep = "_gen") %>% 
  separate(sim_rep.id, c("sim", "rep.id"), sep = "_replicate")


# calcualte average introgressed frequency per simulation per generation
d  <- b %>% 
   group_by(sim,pos_absolute,chrom,gen,pos_on_chrom,recomb) %>% 
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
  filter(sim == "sel-recomb-var", gen != "0001") %>% 
  ggplot(aes(x=pos_absolute, y=avg.frq, group = sim, col = sim)) +
  geom_point() + 
  facet_wrap(~gen, scales = "free_y")



cor_text = d %>% 
  filter(sim == "sel-recomb-var", gen != "0001") %>% 
  group_by(sim,gen) %>% 
  summarise(correlation = round(cor.test(avg.frq,recomb)$estimate, 3) ) %>% #,
  #p.value = round(cor.test(freq,recomb)$p.value, 6)) %>% 
  mutate(text = paste0("r = ", correlation)) #, "\n", "P = ", p.value))


p <- d %>%
  filter(sim == "sel-recomb-var", gen != "0001") %>% 
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
# calculate and visualize power spectrum
########################################

# calculate power spectrum per simulation, generation
frq_ps <- ddply(d, .(sim,gen), ps)
colnames(frq_ps)[-c(1:2)] <- rev(round((1e9/2^(0:9))/1000000))  

# reformat 
frq_ps <- frq_ps %>% 
  gather(key = scale, value = ps, -c(sim,gen))

# order scale as factor for plotting
frq_ps$scale <- as.factor(frq_ps$scale)
frq_ps$scale <- factor(frq_ps$scale, levels = as.character(rev(round((1e9/2^(0:9))/1000000))))

# plot power spectrum
frq_ps %>% 
  filter(str_detect(sim, "sel")) %>% 
  ggplot(aes(x = as.factor(scale), y = ps, group = sim, color = sim, pch = sim)) + 
  geom_line(position=position_dodge(w=0.1)) +
  geom_point() +
  facet_wrap(~gen) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Scale (Mb)")



#################################################################################################
# Wavelet correlation analysis
#################################################################################################

# compute correlations using average allele frequency
wav_cor <- d %>% 
  filter(sim == "sel-recomb-var", gen != "0001") %>% 
  ddply("gen", wavelet_correlation, var1 = "recomb", var2 = "avg.frq") %>% 
  select_if(~sum(!is.na(.)) > 0) %>% 
  gather(key = temp, value = temp_val, -c(gen,total.cor)) %>% 
  separate(col = temp, into = c("thing", "scale"), sep = ":scale") %>% 
  spread(thing, temp_val) %>%
  mutate(weighted.cor = ((`alpha-avg.frq`*`alpha-recomb`)^0.5)*`detail-cor`)
  
wav_cor$scale <- factor(wav_cor$scale, levels = as.character(1:10))

bars <- wav_cor %>%
  ggplot(aes(fill = scale, x = gen, y = weighted.cor)) + 
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_brewer(palette="Dark2")


total <- wav_cor %>% 
  group_by(gen) %>% 
  summarise(tot = sum(weighted.cor, na.rm = T)) %>% 
  ggplot(aes(x = gen, y = tot, group = 1)) + 
  geom_line()

total

ggplot(wav_cor, aes(x = gen, y = total.cor, group = 1)) + geom_line()



smooth_cor <- d %>% 
  filter(sim == "sel-recomb-var") %>% 
  ddply("gen", wavelet_correlation, var1 = "recomb", var2 = "avg.frq", detail = F) %>% 
  filter(gen != "0001") %>% 
  select_if(~sum(!is.na(.)) > 0) %>% 
  gather(key = scale, value = correlation, contains("scale"))

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


# decompose correlation

head(detail_cor)


wr  <- wd(x$recomb, family = "DaubExPhase", filter.number = 1)
wf <- wd(x$avg.frq, family = "DaubExPhase", filter.number = 1)

n <- 1024
s <- vector()
for(i in 1:wr$nlevels){
  wx <- accessD(wr, level = i-1)
  wy <- accessD(wf, level = i-1)
  s[i] <- wx %*% wy
}

sum(s)








d <- data.frame(temp = c("smooth:scale1", "detail:scale2", "alpha-avg.var1:scale3","alpha-var2:scale4"))
d
separate(d, col = temp, into = c("thing", "scale"), sep = ":")
