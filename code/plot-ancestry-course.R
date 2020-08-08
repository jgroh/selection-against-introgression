library(tidyverse)
library(tidyselect)

a <- read.table("results/single-chrom/scheme01/ancestry-results.txt", row.names = 1)
a <- as.data.frame(t(a))
a$pos <- seq(1,1000)

gens <- c("gen005","gen010","gen050","gen100","gen500")


  
for(i in 1:10) {
  # loop over replicates
  
  for(j in 0:(length(gens)-1)) {
    # loop over time steps
   
    # select appropriate columns to take difference
    freq2 <- a[paste0("replicate-", i, "-", gens[j+1])]
    
    if (j != 0) {
      freq1 <- a[paste0("replicate-", i, "-", gens[j])]
      n <- paste0("r", i, "_", gens[j+1], "-", gens[j])
    
    } else {
      freq1 <- rep(0.2, 1000)
      n <- paste0("r", i, "_", gens[j+1], "-", "gen000")
    }
    
    # calculate frequency difference for time step and add new column
    a[, ncol(a) + 1] <- freq2-freq1
    names(a)[ncol(a)] <- n
  }
}

# reorganize data
b <- a %>%
  select(-contains("replicate")) %>% # removes cols with raw freq values
  gather(key = time.step.id,
         value = delta.freq, -pos) %>%
  separate(time.step.id, c("rep.id", "step.num"), sep = "_")

# make plot 
b %>%
  ggplot(aes(pos, delta.freq, colour=rep.id)) +
  #geom_point() +
  geom_smooth(method = "loess", span = .1, se = FALSE) +
  facet_grid(step.num~.)
  
