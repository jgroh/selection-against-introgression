library(tidyverse)

anc_005 <- read.csv("results/ancestry-005.csv", header = F, col.names = "gen005")
anc_010 <- read.csv("results/ancestry-010.csv", header = F, col.names = "gen010")
anc_050 <- read.csv("results/ancestry-050.csv", header = F, col.names = "gen050")
anc_100 <- read.csv("results/ancestry-100.csv", header = F, col.names = "gen100")
anc_500 <- read.csv("results/ancestry-500.csv", header = F, col.names = "gen500")


a <- cbind(ancestry_005, ancestry_010, ancestry_050, ancestry_100, ancestry_500)
tail(a)
a$pos <- seq(1,1000)
colnames(a)
head(a)

# course of allele frequency
a %>%
  gather(key = gen, 
         value = ancestry, gen005, gen010, gen050, gen100, gen500) %>%
  ggplot(aes(x = pos, y=ancestry)) + facet_grid(gen~.) + 
  geom_point() + geom_smooth(method = "loess", span = .1, se = FALSE)

# course of allele frequency change
a$d1 <- with(ancestry_course, gen005 - rep(0.2, 1000))
a$d2 <- with(ancestry_course, gen010 - gen005)
a$d3 <- with(ancestry_course, gen050 - gen010)
a$d4 <- with(ancestry_course, gen100 - gen050)
a$d5 <- with(ancestry_course, gen500 - gen100)

a %>%
  gather(key = delta_freq, 
         value = ancestry, d1, d2, d3, d4, d5) %>%
  ggplot(aes(x = pos, y=ancestry)) + facet_grid(delta_freq ~.) + 
  geom_point() + geom_smooth(method = "loess", span = .1, se = FALSE)

