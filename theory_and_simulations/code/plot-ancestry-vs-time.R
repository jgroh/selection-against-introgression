a <- read.table("results/single-chrom/scheme01/avgs-over-time.txt", stringsAsFactors = F)

colnames(a) <- c("replicate", 1:500)

b <- a %>%
  gather(key = gen, value = ancestry, -replicate, convert = T) %>%
  separate(ancestry, c("q1", "q2", "q3", "q4"), sep = ",", convert = T) %>%
  gather(key = portion, value = ancestry, q1, q2, q3, q4, factor_key = TRUE)

ggplot(b, aes(x=gen, y=ancestry, colour=portion, group=interaction(portion,replicate))) + geom_line()

       