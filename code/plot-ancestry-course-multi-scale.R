#!/usr/bin/Rscript --vanilla

library(tidyverse)
library(tidyselect)

#args <- commandArgs(trailingOnly = TRUE)
#scheme <- args[1]
#filename <- paste0("results/single-chrom/", scheme, "/ancestry-results.txt")
filename <- "results/multi-scale/ancestry-results.txt"
a <- read.table(filename, row.names = 1)

#a <- read.table("workspace/selection-against-introgression/results/single-chrom/scheme01/ancestry-results.txt", row.names = 1)
#scheme <- "scheme01"

a <- as.data.frame(t(a))
a$pos <- rep(1:128, 8)
a$chrom <- rep(1:8, each=128)
colnames(a)


# recomb rate vector
baseRate <- 1e-8
r1 <- c(c(baseRate*c(0.1,0.1,1,1)*c(10,1,10,5)),c(baseRate*c(0.01,0.01,0.1,0.1)*c(10,1,10,5)));
recomb.rates <- rep(rep(r1, each = 32),4)
a$recomb <- recomb.rates

# show recombination scheme
plot(a$recomb, type = "n", xlab = "position", ylab = "rec. rate")
lines(a$recomb)
for(i in 1:7){
  lines(x = c(128*i,128*i), y =c(0,1e-7), lwd = 3, lty = 2)
}


# plot frequency trajectory directly

## reorganize data
b  <- a %>%
  gather(key = rep_gen,
         value = freq, -c(pos,chrom,recomb)) %>%
  separate(rep_gen, c("rep.id", "gen"), sep = "-gen")


## make plot of correlation on a per base scale across generations
b %>%
  ggplot(aes(jitter(recomb), freq)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_continuous(trans = "sqrt" ) +
  scale_x_continuous(trans = "log10" ) +
  facet_grid(gen~., scales = "free_y")

# wavelet transformation
library(wavethresh)

# loop over generations, average frequency over replicates

r3 <- vector()
r4 <- vector()
p3 <- vector()
p4 <- vector()
gens <- c("002", "005", "010", "050", "100", "500")
for(i in 1:length(gens)){
  d <- subset(b, gen == gens[i])

  d$pos2 <- 1:1024
  avg.freq <- tapply(d$freq, d$pos2, mean)
  recomb <- subset(d, rep.id =="replicate-1")$recomb

  w1 <- wd(recomb, family = "DaubExPhase", filter.number = 1)
  w2 <- wd(avg.freq, family  =  "DaubExPhase", filter.number = 1)
  plot(w1)
  plot(w2)
  
  x3 <- accessD(w1, level = 3)
  y3 <- accessD(w2, level = 3)
  r3[i] <- cor.test(x3,y3)$estimate
  p3[i] <- cor.test(x3,y3)$p.value
  
  x4 <- accessD(w1, level = 4)
  y4 <- accessD(w2, level = 4)
  r4[i] <- cor.test(x4,y4)$estimate
  p4[i] <- cor.test(x4,y4)$p.value

}

plot(r3, ylim = c(-1,1), col = "red", lwd = 2, xlab = "generation", ylab = expression(paste("Kendall's  ", tau)),  xaxt = "n")
lines(r3, lwd = 2, col = "red")
axis(1, at = 1:6, labels = gens)

points(r4, col = "blue", lwd =2, add = T)
lines(r4, lwd = 2, col = "blue")

legend("bottomright", col = c("red", "blue"), legend= c("1.25 Mb", "0.625 Mb"), pch = 1, box.lwd = 0, lwd = 2)

text(x = 1:6, y = r3+.1, labels = round(p3,4))
text(x = 1:6, y = r4+.1, labels = round(p4,4))

# try to plot continuous thing


d <- subset(b, gen == "500")
d$pos2 <- 1:1024
avg.freq <- tapply(d$freq, d$pos2, mean)
plot(avg.freq, type = "n", ylab = "Avg. introgressed ancestry", xlab="position")
lines(avg.freq)

library(Rwave)
Mod(DOG(avg.freq,9,nvoice=10,twoD=TRUE,plot=T,moments=2)) 
Mod(DOG(recomb,9,nvoice=10,twoD=TRUE,plot=T,moment=2)) 
?Mod





## make plot of allele frequency over time
b %>%
  ggplot(aes(seq, freq, colour=rep.id)) +
  #geom_point() +
  geom_smooth(aes(group = rep.id), method = "loess", span = .1, se = FALSE, size = .5) +
  facet_grid(gen~., scales = "free_y")






































# # plot frequency differentials between time steps
# 
# ## calculate differentials 
# 
# gens <- c("gen002", "gen005","gen010","gen050","gen100","gen500")
# for(i in 1:10) {
#   # loop over replicates
#   
#   for(j in 0:(length(gens)-1)) {
#     # loop over time steps
#    
#     # select appropriate columns to take difference
#     freq2 <- a[paste0("replicate-", i, "-", gens[j+1])]
#     
#     if (j == 0) {
#       freq1 <- rep(0.2, 100)
#       n <- paste0("r", i, "_", gens[j+1], "-", "gen001")
#     
#     } else {
#       freq1 <- a[paste0("replicate-", i, "-", gens[j])]
#       n <- paste0("r", i, "_", gens[j+1], "-", gens[j])
#     }
#     
#     # calculate frequency difference for time step and add new column
#     a[, ncol(a) + 1] <- freq2-freq1
#     names(a)[ncol(a)] <- n
#   }
# }
# 
# 
# ## reorganize data
# d <- a %>%
#   select(-contains("replicate")) %>% # removes cols with raw freq values
#   gather(key = time.step.id,
#          value = delta.freq, -pos) %>%
#   separate(time.step.id, c("rep.id", "step.num"), sep = "_")
# 
# ## make plot 
# d %>%
#   ggplot(aes(pos, delta.freq, colour=rep.id)) +
#   #geom_point() +
#   geom_smooth(method = "loess", span = .1, se = FALSE, size = .5) +
#   facet_grid(step.num~., scales = "free_y")
# 
# ggsave(paste0("/Users/jeff/workspace/selection-against-introgression/results/plots/", scheme, "-delta-freq-trajectory.pdf"))
# 
#   
