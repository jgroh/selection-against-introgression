library(data.table)
library(tidyverse)
library(wavethresh)
library(waveslim)
library(stringi)

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

# MODWT ----------------
n.levels <- floor(log2(length(unique(chrom1.interp$gen_pos))))

# run wavelet decomposition on individuals
ind.modwt <- chrom1.interp[, waveslim::modwt(ind_frq_interp, "haar", 
                                                  n.levels = n.levels), by = ID]
# compute wavelet variance
wav_var <- function(x){sum(x^2)/(length(x))}
detail.cols <- names(ind.modwt)[grep("s", names(ind.modwt),invert=T)]
ind.wav.var <- ind.modwt[,..detail.cols][,lapply(.SD,wav_var), by = ID]
setnames(ind.wav.var, old = paste0("d",1:n.levels), new = as.character(1:n.levels))
ind.wav.var <- melt(ind.wav.var, id.vars = c("ID"), 
                    variable.name = "scale",
                    value.name = "variance")

# compute mean across individuals
ind_mean_wav_var <- ind.wav.var[, mean(variance), by = scale]
setnames(ind_mean_wav_var, "V1", "variance")


# run wavelet decomp on population mean
pop.mean.modwt <- chrom1.interp[ID==ID[1],
                               waveslim::modwt(x=pop_mean_interp, wf="haar", 
                               n.levels = n.levels)]
pop.mean.wav.var <- rbindlist(list(pop.mean.modwt))
smooth.col <- grep("s",names(pop.mean.wav.var))
pop.mean.wav.var <- pop.mean.wav.var[,-smooth.col,with=F][,lapply(.SD,wav_var)]
setnames(pop.mean.wav.var, paste0("d",1:n.levels), as.character(1:n.levels))


#pop.mean.wv <- setDT(wave.variance(pop.mean.modwt, type = "nongaussian"))
#setnames(pop.mean.wv, "wavevar", "variance")
#pop.mean.wv <- pop.mean.wv[-nrow(pop.mean.wv)]
#pop.mean.wv[,"scale" := as.numeric(1:15)]

pop.mean.wav.var <- melt(pop.mean.wav.var, measure.vars = as.character(1:n.levels),
     variable.name = "scale", value.name = "variance")
pop.mean.wav.var[,scale:=as.numeric(scale)]

ggplot(ind.wav.var, aes(x=as.numeric(scale), y=variance)) + 
  theme(legend.position = "none") + 
  geom_line(aes(group = ID, color = ID, alpha = 0.8)) +
  geom_line(data = pop.mean.wav.var, size = 1) + 
  geom_line(data = ind_mean_wav_var, color = "red", size = 1) + 
  #geom_ribbon(data = pop.mean.wv, aes(ymin=lower,ymax=upper)) +
  labs(x = "Scale log2 (N x Morgans)", y = "Wavelet variance")#   +
  geom_line(data = pop.mean.roll.wav.var, color = "blue", size = 1) 

# Circular permutation of sequence within individuals ---------
roll <- function(x){
  n <- sample(1:length(x),1)
  if(n == 0)
    return(x)
  c(tail(x,n), head(x,-n))
}

# create column that will get rolled 
chrom1.interp[, pos_roll := gen_pos, by = .(chr,ID)]

# data table that will house permuted pop mean wav var
pop.mean.roll.wav.var <- data.table(scale=1:n.levels)

# do this stuff many times:
for(i in paste0("rep",1:1000)){
  # circularly shift individual frequencies
  chrom1.interp[, pos_roll := roll(pos_roll), by = .(chr,ID)]
  # recalculate population mean
  chrom1.interp[, pop_mean_roll := mean(ind_frq_interp), by =.(chr, pos_roll)]
  
  # wavelet transform of rolled pop mean
  pop.mean.roll.modwt <- chrom1.interp[ID==ID[1], modwt(pop_mean_roll, "haar", n.levels=n.levels)]
  tmp <- rbindlist(list(pop.mean.roll.modwt))
  smooth.col <- grep(paste0("s",n.levels),names(tmp))
  tmp <- tmp[,-smooth.col,with=F][,lapply(.SD,wav_var)]
  pop.mean.roll.wav.var[, (i) := t(tmp)]
}

pop.mean.roll.wav.var <- melt(pop.mean.roll.wav.var, id.vars = "scale", 
                              variable.name = "replicate", value.name="variance")
pop.mean.null.dist <- pop.mean.roll.wav.var[, .(min(variance),max(variance)), by=scale]
setnames(pop.mean.null.dist, c("V1","V2"), c("min", "max"))

pop.mean.plot.data <- merge(pop.mean.null.dist,pop.mean.wav.var,by="scale")

ggplot(pop.mean.plot.data,aes(x=scale,y=variance)) +
  geom_line() + 
  geom_ribbon(aes(ymin=min,ymax=max), alpha = 0.5) + 
  labs(x = "Scale log2 (N x Morgans)", y = "Wavelet variance")#   +


save(pop.mean.plot.data, file = "datasets/schumer_etal_2019/TOTO_chrom1_pop_mean_WV.R")











# test on single individual: first plot wavelet variance for one ind.
ggplot(ind.wav.var[ID %in% ID[1:2]], aes(x=as.numeric(scale), y=variance)) + 
  geom_line(aes(group = ID, color = "red", alpha = 0.8,size=2)) +
  geom_line(data=test.ind.wav.var, aes(group = ID), color = "black")# + 
 # geom_line(data=test.pop.wav.var) + 
  #geom_line(data=test.ind.re_wav.var, aes(group = ID), color = "green")
  
# permute

permute_wav_coeff <- function(x){
  x[1]
}
str(ind.modwt[1,V1][[1]])


test.ind.permuted <- ind.modwt[ID %in% unique(ID)[1:2], lapply(.SD, sample, replace = F), by = ID]
test.ind.wav.var <- test.ind.permuted[,-c("s15")][,lapply(.SD,wav_var), by = ID]
setnames(test.ind.wav.var, old = paste0("d",1:15), new = as.character(1:15))
test.ind.wav.var <- melt(test.ind.wav.var, id.vars = c("ID"), 
                    variable.name = "scale",
                    value.name = "variance")
test.ind.wav.var
test.pop.wav.var <- test.ind.wav.var[, mean(variance), by = scale]
setnames(test.pop.wav.var, "V1", "variance")

# have determined that the permutation is not the issue
# it preserves the variance as desired. 
l <- as.list(test.ind.permuted[ID==ID[1]])
atts <- list("names" = names(l),
             "class" = "modwt", 
             "wavelet" = "haar", 
             "boundary"="periodic")
attributes(l) <- atts
l$ID <- NULL
l.x <- imodwt(l)
re_mod <- modwt(l.x, "haar", n.levels=n.levels)

# key problem: these two things should be equal!!!
test.ind.permuted$d1[1:10]
re_mod$d1[1:10]



# function to return inverse modwt for each individual 
inv_modwt <- function(x){
  l <- as.list(x)
  atts <- list("names" = names(l),
               "class" = "modwt", 
               "wavelet" = "haar", 
               "boundary"="periodic")
  attributes(l) <- atts
  l$ID <- NULL
  return(imodwt(l))
}

test.ind.imodwt <- test.ind.permuted[,inv_modwt(.SD), by = ID]
setnames(test.ind.imodwt, "V1", "inverted_frq")

# redo modwt to make sure this worked correctly
test.ind.re_modwt <- test.ind.imodwt[, modwt(inverted_frq, "haar", 
                                                 n.levels = n.levels), by=ID]
test.ind.re_wav.var <- test.ind.re_modwt[,-c("s15")][,lapply(.SD,wav_var), by = ID]
setnames(test.ind.re_wav.var, old = paste0("d",1:15), new = as.character(1:15))
test.ind.re_wav.var <- melt(test.ind.re_wav.var, id.vars = c("ID"), 
                         variable.name = "scale",
                         value.name = "variance")








inv.ind.modwt.permuted <- ind.modwt.permuted[, inv_modwt(.SD), by = ID]
setnames(inv.ind.modwt.permuted, "V1", "inverted_frq")

# add column to index position
inv.ind.modwt.permuted[, position := seq_len(nrow(.SD)), by = ID]

# recompute population mean
pop.mean.permuted <- inv.ind.modwt.permuted[, mean(inverted_frq), by = position]
setnames(pop.mean.permuted, "V1", "pop_mean")

# rerun decomp on individuals (safety check)

# run wavelet decomposition on individuals
ind.modwt.p <- inv.ind.modwt.permuted[, modwt(x=inverted_frq, "haar", n.levels = n.levels),
                           by = ID]
ind.wav.var.p <- ind.modwt.p[,-c("s15")][,lapply(.SD,wav_var), by = ID]
setnames(ind.wav.var.p, old = paste0("d",1:15), new = as.character(1:15))
ind.wav.var.p <- melt(ind.wav.var.p, id.vars = c("ID"), 
                    variable.name = "scale",
                    value.name = "variance")

# compute mean across individuals
ind_mean_wav_var.p <- ind.wav.var.p[, mean(variance), by = scale]
setnames(ind_mean_wav_var.p, "V1", "variance")

ggplot(ind.wav.var.p, aes(x=as.numeric(scale), y=variance)) + 
  theme(legend.position = "none") + 
  geom_line(aes(group = ID, color = ID, alpha = 0.8)) +
  geom_line(data = ind_mean_wav_var.p, color = "red", size = 1) + 
  labs(x = "Scale log2 (N x Morgans)", y = "Wavelet variance")  







# rerun wavelet decomp on population mean
pop.mean.modwt.p <- pop.mean.permuted[,modwt(pop_mean, "haar", 
                                                 n.levels = n.levels)]
pop.mean.wav.var.p <- rbindlist(list(pop.mean.modwt.p))
pop.mean.wav.var.p <- pop.mean.wav.var.p[,-c("s15")][,lapply(.SD,wav_var)]
setnames(pop.mean.wav.var.p, paste0("d",1:15), as.character(1:15))

pop.mean.wav.var.p <- melt(pop.mean.wav.var.p, measure.vars = as.character(1:15),
                         variable.name = "scale", value.name = "variance")

ggplot(pop.mean.wav.var.p, aes(x = as.numeric(scale),y=variance)) + geom_line()




# DWT ----------------
w2 <- wd(a$mean_frq_interp[1:2^15], filter.number = 1, family = "DaubExPhase")
v <- NA
for(i in 0:14){
  v[i] <- sum(accessD(w2,level = i)^2)/2^15
}
points(rev(v), col = "red")
# sum((a$mean_frq_interp - mean(a$mean_frq_interp))^2)/length(a$mean_frq_interp)


par(mar = c(0,0,0,0), mfcol=c(5,3), pty="m", mar=c(5-2,4,4-2,2))

for(i in 1:15){
  plot.ts(w[[i]], ylab=names(w)[i], mar = c(0,0,0,0))
}, boundary = "reflection"









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
