#!/usr/bin/Rscript --vanilla

library(tidyverse)
library(tidyselect)

args <- commandArgs(trailingOnly = TRUE)
scheme <- args[1]
filename <- paste0("results/single-chrom/", scheme, "/ancestry-results.txt")

a <- read.table(filename, row.names = 1)
#a <- read.table("results/single-chrom/scheme01/ancestry-results.txt", row.names = 1)
a <- as.data.frame(t(a))
a$pos <- seq(1,1000)


# plot frequency trajectory directly

## reorganize data
b <- a %>%
  gather(key = rep_gen,
         value = freq, -pos) %>%
  separate(rep_gen, c("rep.id", "gen"), sep = "-gen")

## make plot 
b %>%
  ggplot(aes(pos, freq, colour=rep.id)) +
  #geom_point() +
  geom_smooth(method = "loess", span = .1, se = FALSE) +
  facet_grid(gen~.)

ggsave(paste0("results/plots/", scheme, "-freq-trajectory.pdf"))



# plot frequency differentials between time steps

## calculate differentials 

gens <- c("gen005","gen010","gen050","gen100","gen500")
for(i in 1:10) {
  # loop over replicates
  
  for(j in 0:(length(gens)-1)) {
    # loop over time steps
   
    # select appropriate columns to take difference
    freq2 <- a[paste0("replicate-", i, "-", gens[j+1])]
    
    if (j == 0) {
      freq1 <- rep(0.2, 1000)
      n <- paste0("r", i, "_", gens[j+1], "-", "gen000")
    
    } else {
      freq1 <- a[paste0("replicate-", i, "-", gens[j])]
      n <- paste0("r", i, "_", gens[j+1], "-", gens[j])
    }
    
    # calculate frequency difference for time step and add new column
    a[, ncol(a) + 1] <- freq2-freq1
    names(a)[ncol(a)] <- n
  }
}


## reorganize data
d <- a %>%
  select(-contains("replicate")) %>% # removes cols with raw freq values
  gather(key = time.step.id,
         value = delta.freq, -pos) %>%
  separate(time.step.id, c("rep.id", "step.num"), sep = "_")

## make plot 
d %>%
  ggplot(aes(pos, delta.freq, colour=rep.id)) +
  #geom_point() +
  geom_smooth(method = "loess", span = .1, se = FALSE) +
  facet_grid(step.num~.)

ggsave(paste0("results/plots/", scheme, "delta-freq-trajectory.pdf"))

  
