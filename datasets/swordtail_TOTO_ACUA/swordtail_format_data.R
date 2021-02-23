library(data.table)
library(tidyverse)
library(wavethresh)
library(waveslim)

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]
chrom <- args[2]

# edit for snakemake to apply to multiple pops
par1 <- fread("datasets/schumer_etal_2019/Population1_TOTO/ancestry-probs-par1_TOTO_allgroups.tsv")
par2 <- fread("datasets/schumer_etal_2019/Population1_TOTO/ancestry-probs-par2_TOTO_allgroups.tsv")

# Merge two files by ID and reformat ----------------------------------------
gnom <- merge(par1, par2, by = "V1", suffixes = c(":gen11", ":gen22"))
setnames(gnom, "V1", "ID")

# subset to chr 1 here
chrom1.cols <- c("ID", names(gnom)[grep("group1:", names(gnom))])
chrom1 <- gnom[, ..chrom1.cols]

# move all post. probability entries to a single column
# (throws warning about missing data, but OK)
chrom1 <- melt(chrom1, id.vars = "ID", value.name = "post_prob")

# separate by chromosome, position, genotype
chrom1[, c("chr", "phys_pos", "genotype") := transpose(stri_split_fixed(variable, ":"))]


chrom1[,variable:=NULL]
chrom1[, phys_pos := as.integer(phys_pos)]

# move genotypes to separate columns
chrom1 <- dcast(chrom1, ... ~ genotype, value.var = "post_prob")

# add column for frequency of allele 1 in individual weighted by post. probs. 
chrom1[, frq := (gen11 + 0.5*(1 - sum(gen11, gen22))), by = .(ID,chr,phys_pos)]

# Interpolate ancestry in individuals -------------------------

# interpolate ancestry at each SNP position
chrom1[, c("ind_frq_interp") := 
         approx(x = phys_pos, y = frq, xout = phys_pos, rule = 2)$y, 
       by = .(ID,chr)]

# average over individuals to get population mean
chrom1[, "pop_mean" := mean(ind_frq_interp), by = .(chr, phys_pos)]


# examine ancestry pattern for a couple of individuals
#a[V1 %in% unique(a$V1)[1:2]] %>% ggplot(aes(x = phys_pos, y = frq, group = V1, color = V1)) + geom_point()

# Recombination map -----------------------

# read recomb map for chr X (will update for parallelization)
# edited headers of these files before reading in
r.chrom1 <- fread("datasets/schumer_etal_2019/lg-maps/group1.txt")


# generate vector of genetic distance using recombination LD map

# propogate mean values for intervals to all intervening bps 
r.chrom1.vec <- r.chrom1[, .(rate = rep(mean, times = (right_snp - (left_snp) ))), by = right_snp][, -c("right_snp")]

# extend per bp values to ends of chromosome
if(min(r.chrom1$left_snp) > 1){
  front <- data.table(rate = rep(r.chrom1[left_snp == min(left_snp), mean], times=min(r.chrom1$left_snp)))
} else{front <- data.table()}
if(max(r.chrom1$right_snp) < max(chrom1$phys_pos)){
  back <- data.table(
    rate = rep(r.chrom1[right_snp == max(right_snp), mean],
               times=(max(chrom1$phys_pos)-max(r.chrom1$right_snp))))
} else{back <- data.table()}

# stick together for full rec. map per bp 
r.chrom1.vec <- rbindlist(list(front, r.chrom1.vec, back))

# assign genetic position of SNPS present in the data and per bp rate
chrom1[, c("gen_pos", "rec_rate") := 
         list(cumsum(r.chrom1.vec$rate)[phys_pos], r.chrom1.vec$rate[phys_pos]),
       by = .(chr)]
# sanity checks
# dim(r.chrom1.vec)[1] == max(chrom1$phys_pos)
# plot gen distancew against physical
# ind <- chrom1[ID == ID[1]]$phys_pos
# plot(y = cumsum(r.chrom1.vec$rate)[ind], x = chrom1[ID == ID[1]]$phys_pos)

# Interpolate ancestry to evenly spaced genetic coordinates -------------------------

# plot population mean of minor parent ancestry along chromosome
chrom1[ID == ID[1]] %>% ggplot(aes(x = phys_pos, y = 1-pop_mean)) + geom_point()

# plot overall correlation between rec. and population mean ancestry on per bp level
chrom1[ID == ID[1], .(pop_mean, rec_rate)] %>% 
  ggplot(aes(x = log(rec_rate), y = 1-pop_mean)) + 
  geom_point() + geom_smooth(method = "lm") + 
  labs(x = "log(rho / bp)", y = "Minor parent ancestry")

chrom1[ID == ID[1], cor(rec_rate, 1-pop_mean)]

# Average rec. and ancestry in 50 kb windows -----------------

# calculate number of windows by chromosome
chrom1[, num.win := ceiling(max(phys_pos)/50000), by = chr]

# add variable for window number
chrom1[, window_50kb := 
         cut(phys_pos, breaks = num.win, 
             include.lowest = TRUE, labels = 1:num.win),
       by = chr]
# get average by window    
chrom1.50kb <- chrom1[, lapply(.SD, mean), 
                      by = window_50kb, .SDcols = c("rec_rate", "pop_mean")]
head(chrom1.50kb)

# plot window averages binned by recombination quantile
chrom1.50kb[, .SD, by = cut(rec_rate, breaks=c(quantile(rec_rate, probs = seq(0, 1, by = 0.20))), 
                            labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"), include.lowest = T)] %>%
  ggplot(aes(x = cut, y = 1 - pop_mean)) + geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Chromosome 1", x= "Recomb. rate quantile", y="Minor parent freq") +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=3, color="red") 


# plot minor parent ancestry vs physical distance
chrom1[ID == ID[1]] %>% ggplot(aes(y = 1-pop_mean, x = phys_pos))+geom_point()

# vs. genetic distance (non-interpolated genetic distances)
chrom1[ID == ID[1]] %>% ggplot(aes(y = 1-pop_mean, x = gen_pos)) + 
  geom_point() + labs(x = "Genetic distance (4N Morgans)", 
                      y = "Avg. minor parent ancestry")


# interpolate of population mean ancestry at evenly spaced genetic distances

# what unit to interpolate at?
# check # of SNPs / total genetic distance (in units 4*N*Morgans):
# length(unique(chrom1$phys_pos))/max(chrom1$gen_pos) #~= 0.25 
# so we can interpolate to unit distance of N*Morgans so that we on average have 1 snp per unit distance

xout.l <- ceiling(max(chrom1$gen_pos))/4 # # of interpolation points to give unit distance of N*M

xout = seq(0, max(chrom1$gen_pos), length.out = xout.l)

# interpolate individual ancestry at genetic coordinates
chrom1.interp <- chrom1[, approx(x = gen_pos, 
                                 y = ind_frq_interp, 
                                 xout = xout, rule = 2), 
                        by = .(ID,chr)]
setnames(chrom1.interp, c("x", "y"), c("gen_pos","ind_frq_interp"))

# compute population mean
chrom1.interp[, pop_mean_interp := mean(ind_frq_interp), 
              by = .(chr, gen_pos)]

chrom1.interp[ID == ID[1]] %>% ggplot(aes(x = gen_pos, y = 1-pop_mean_interp)) + geom_point()

save(chrom1, chrom1.50kb, chrom1.interp, file = "datasets/schumer_etal_2019/TOTO_chrom1.Rdata")
load(file="datasets/schumer_etal_2019/TOTO_chrom1.Rdata")

## End Formatting

