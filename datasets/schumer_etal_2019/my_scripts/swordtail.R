library(data.table)
library(tidyverse)
library(wavethresh)

# read genotype files
load("datasets/schumer_etal_2019/swordtail_chr1.Rdata")
#a1 <- fread("datasets/schumer_etal_2019/Population1_TOTO/ancestry-probs-par1_TOTO_allgroups.tsv")
#a2 <- fread("datasets/schumer_etal_2019/Population1_TOTO/ancestry-probs-par2_TOTO_allgroups.tsv")

# subset to chr X
#cols1 <- names(a1)[grep("group1:", names(a1))]
#a1 <- a1[, .SD, .SDcols = c("V1",cols1)]
#a2 <- a2[, .SD, .SDcols = c("V1", cols1)]

# read recomb map for chr X
# probably want to edit header of this file before reading it in
r1 <- fread("datasets/schumer_etal_2019/lg-maps/group1.txt")
#head(r1)
#max(r1$right_snp)

# merge two files by ID and reformat ----------------------------------------
a <- merge(a1, a2, by = "V1", suffixes = c(":gen11", ":gen22"))
setnames(a, "V1", "ID")
a <-melt(a, id.vars = "ID", value.name = "post.prob")
a[,c("chr", "phys_pos", "genotype") := tstrsplit(variable, ":", fixed = TRUE)
  ][,variable:=NULL][, phys_pos := as.integer(phys_pos)]
a <- dcast(a, ... ~ genotype, value.var = "post.prob")

# add column for frequency of allele 1 in individual weighted by post. probs. 
a[,frq := (gen11 + 0.5*(1 - sum(gen11, gen22))), by = .(ID,chr,phys_pos)]

# examine ancestry pattern for a couple of individuals
#a[V1 %in% unique(a$V1)[1:2]] %>% ggplot(aes(x = phys_pos, y = frq, group = V1, color = V1)) + geom_point()

# Recombination map -----------------------
# generate vector of genetic distance using recombination LD map
head(r1)
r <- r1[, .(rate = rep(mean, times = (right_snp - (left_snp) ))), by = right_snp][, -c("right_snp")]
if(min(r1$left_snp) > 1){
  front <- data.table(rate = rep(r1[left_snp == min(left_snp), mean], times=min(r1$left_snp)))
} else{front <- data.table()}
if(max(r1$right_snp) < max(a$phys_pos)){
  back <- data.table(
    rate = rep(r1[right_snp == max(right_snp), mean],
               times=(max(a$phys_pos)-max(r1$right_snp))))
} else{back <- data.table()}

r <- rbindlist(list(front, r, back))

# assign genetic position of SNPS present in the data and per bp rate
a[, c("gen_pos", "rec_rate") := 
    list(cumsum(r$rate)[phys_pos], r$rate[phys_pos]),
  by = ID]
# sanity checks
# dim(r)[1] == max(a$phys_pos)
# plot gen distancew against physical
# ind <- a[V1 == V1[1]]$phys_pos
# plot(y = cumsum(r$rate)[ind], x = a[V1 == V1[1]]$phys_pos)

# Interpolate ancestry -------------------------

# interpolate ancestry at each SNP position
a[, c("frq_interp") := 
    approx(x = phys_pos, y = frq, xout = phys_pos, rule = 2)$y, 
  by = ID]

# average over individuals
a[, "mean_frq_interp" := mean(frq_interp), by = .(chr, phys_pos, rec_rate)]

# plot population mean of minor parent ancestry along chromosome
a[ID == ID[1]] %>% ggplot(aes(x = phys_pos, y = 1-mean_frq_interp)) + geom_point()

# plot overall correlation between rec. and population mean ancestry on per bp level
a[ID == ID[1], .(mean_frq_interp, rec_rate)] %>% 
  ggplot(aes(x = log(rec_rate), y = 1-mean_frq_interp)) + 
  geom_point() + geom_smooth(method = "lm") + 
  labs(x = "log(rho / bp)", y = "Minor parent ancestry")

a[, cor(rec_rate, 1-mean_frq_interp), by = ID]

# average rec. and ancestry in 50 kb windows
num.win <- ceiling(max(a$phys_pos)/50000)
a.50kb <- a[ID==ID[1], lapply(.SD, mean), 
  by = cut(phys_pos, breaks = num.win, include.lowest = TRUE, labels = 1:num.win), 
  .SDcols = c("rec_rate","mean_frq_interp")]
setnames(a.50kb, "cut", "window")
a.50kb

# plot window averages binned by recombination quantile
a.50kb[, .SD, by = cut(rec_rate, breaks=c(quantile(rec_rate, probs = seq(0, 1, by = 0.20))), 
              labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"), include.lowest = T)] %>%
  ggplot(aes(x = cut, y = 1 - mean_frq_interp)) + geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Chromosome 1", x= "Recomb. rate quantile", y="Minor parent freq") +
  stat_summary(fun.y=mean, geom="point", shape=18,
              size=3, color="red") 

#save(a1,a2,a,a.50kb, file = "datasets/schumer_etal_2019/swordtail_chr1.Rdata")

# plot minor parent ancestry vs physical distance
a[ID == ID[1]] %>% ggplot(aes(y = 1-mean_frq_interp, x = phys_pos))+geom_point()

# vs. genetic distance
a[ID == ID[1]] %>% ggplot(aes(y = 1-mean_frq_interp, x = gen_pos)) + 
  geom_point() + labs(x = "Genetic distance (4N Morgans)", 
                      y = "Avg. minor parent ancestry")


# interpolate of population mean ancestry on genetic scale
# level grid detail is arbitrary - but can pick largest dyadic sequence lower than the number of SNPs
# so that, in some sense, the average location is averaging over >= 1 SNP

xout = seq(0, max(a$gen_pos), length.out = 2^floor(log2(length(unique(a$phys_pos)))))

a.interp <- a[ID == ID[1], approx(x = gen_pos, y = mean_frq_interp, xout = xout, rule = 2)]
setnames(a.interp, c("x", "y"), c("gen_pos","mean_frq_interp"))
head(a.interp)
a.interp %>% ggplot(aes(x = gen_pos, y = 1-mean_frq_interp)) + geom_point()
#save(a.interp, file = "datasets/schumer_etal_2019/swordtail_chr1.Rdata")

w <- wd(a.interp$mean_frq_interp, filter.number = 1, family = "DaubExPhase")
plot(w)








# how to deal with missing data? 
# just averaging with na.rm=T ignores information about ancestry at linked sites in individuals
# could be biased if reads map don't map equally well to both parental genomes at a location
f <- a[, .(avg_frq = mean(frq, na.rm = T)), by = .(chr,phys_pos)]

f[chr == "group17"][3800:3900] %>% 
  ggplot(aes(x = 3800:3900, y = avg_frq)) + geom_point()
# frequency 
names(a1) == names(a2)
smpl <- sample(names(a1)[-1], size = 100000, replace = F)
a1_smpl <- a1[, ..smpl]
a2_smpl <- a2[, ..smpl]
x <- as.numeric(a1_smpl[1,])
y <- as.numeric(a2_smpl[1,])
par(mar = c(5,4,4,1))
plot(y = y, x = x)
hist(x+y, xlab = "Sum of ancestry posterior probabilities",
     freq = F, main = "", col = "red")

# look at tract lengths on chr 1
cols <- names(a1)[grep("group1:", names(a1))]
a1.c1 <- as.data.frame(a1[, ..cols])
a1.c1.pos <-  as.numeric(sub(".*:", "", x= names(a1.c1)))
a1.c1  <- as.numeric(a1.c1[1,])
plot(y = a1.c1,x = a1.c1.pos, xlab = "Position on chr1", ylab = "Ancestry posterior probability")

mean(a1.c1, na.rm = T)
