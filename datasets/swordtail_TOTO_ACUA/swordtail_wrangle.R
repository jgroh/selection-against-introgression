library(data.table)

# 1. Read and Format Data Files ==========

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]
scaff <- args[2]

# read genotype files 
# probability homozygous for either ancestry
par1 <- fread(paste0("ancestry-probs-par1_allchrs_ACUA_historical_",year,".tsv"))
par2 <- fread(paste0("ancestry-probs-par2_allchrs_ACUA_historical_",year,".tsv"))

# read recomb map for focal chromosome 
r.chrom1 <- fread(paste0("LD_recMap/LD_map_xbirchmanni-COAC-10x-", scaff,".post.txt_mod.bed")) 
setnames(r.chrom1, c("chr","left_bp","right_bp","mean","V1","median","V3"))

# merge by ID 
gnom <- merge(par1, par2, by = "V1", suffixes = c(":gen11", ":gen22"))
setnames(gnom, "V1", "ID")

# subset to focal chromosome
chrom1.cols <- c("ID", names(gnom)[grep(paste0(scaff,":"), names(gnom))])
chrom1 <- gnom[, ..chrom1.cols]

# move genotype probs for all loci to a single column
chrom1 <- melt(chrom1, id.vars = "ID", value.name = "post_prob")

# separate by chromosome, position, genotype
chrom1[, c("chr", "phys_pos", "genotype") := tstrsplit(variable,":",fixed=TRUE)]
chrom1[,variable:=NULL]
chrom1[, phys_pos := as.integer(phys_pos)]

# move genotypes to separate columns
chrom1 <- dcast(chrom1, ... ~ genotype, value.var = "post_prob")

# add column for frequency of allele 1 in individual weighted by post. probs. 
chrom1[, frq := (gen11 + 0.5*(1 - sum(gen11, gen22))), by = .(ID,chr,phys_pos)]


# 2. Impute ancestry at SNP locations ==========

# impute ancestry at each SNP position in all individuals
chrom1[, c("ind_frq") := 
         approx(x = phys_pos, y = frq, xout = phys_pos, rule = 2)$y, 
       by = .(ID,chr)]

# average over individuals to get sample mean
chrom1[, "smpl_mean" := mean(ind_frq), by = .(chr, phys_pos)]

# examine ancestry pattern for a couple of individuals
#a[V1 %in% unique(a$V1)[1:2]] %>% ggplot(aes(x = phys_pos, y = frq, group = V1, color = V1)) + geom_point()


# 3. Create recombination map ==========

# generate vector of genetic distance using recombination LD map
# propogate median values for intervals to all intervening bps 
r.chrom1.allbp <- r.chrom1[, .(rho = rep(median, times = (right_bp - left_bp ))), by = right_bp][, -c("right_bp")]

# extend per bp values to ends of chromosome
if(min(r.chrom1$left_bp) > 1){
  front <- data.table(rho = rep(r.chrom1[left_bp == min(left_bp), median], times=min(r.chrom1$left_bp)))
} else{front <- data.table()}
if(max(r.chrom1$right_bp) < max(chrom1$phys_pos)){
  back <- data.table(
    rho = rep(r.chrom1[right_bp == max(right_bp), median],
               times=(max(chrom1$phys_pos)-max(r.chrom1$right_bp))))
} else{back <- data.table()}

# stick pieces together 
r.chrom1.allbp <- rbindlist(list(front, r.chrom1.allbp, back))

# get genetic distance, divide out 2Ne
Ne2 <- 97739 # see script scaff_lengths.R which compares LD map lengths to cM map lengths

r.chrom1.allbp[, r := rho/Ne2] # rho estimated by LDHelmet is 2Ner
#r.chrom1.allbp[, mb_per_cM := 1/(r*1e2*1e6)]
r.chrom1.allbp[, Morgan := cumsum(r)]

# assign genetic position of SNPS present in the data and per bp rate
chrom1[, c("Morgan", "r") := 
         list(r.chrom1.allbp[,Morgan][phys_pos], r.chrom1.allbp[,r][phys_pos]),
       by = .(chr)]

# sanity check
# dim(r.chrom1.vec)[1] == max(chrom1$phys_pos)
# plot gen distance against physical
# ind <- chrom1[ID == ID[1]]$phys_pos
# plot(y = cumsum(r.chrom1.vec$rate)[ind], x = chrom1[ID == ID[1]]$phys_pos)

# plot population mean of minor parent ancestry along chromosome
#chrom1[ID == ID[1]] %>% ggplot(aes(x = phys_pos, y = 1-sample_mean)) + geom_point()

# vs. genetic distance (non-interpolated genetic distances)
#chrom1[ID == ID[1]] %>% ggplot(aes(y = 1-pop_mean, x = gen_pos)) + 
#  geom_point() + labs(x = "Genetic distance (2N Morgans)", 
#                      y = "Avg. minor parent ancestry")



# Average rec. and ancestry in 50 kb windows -----------------

# calculate number of windows by chromosome
#chrom1[, num.win := ceiling(max(phys_pos)/50000), by = chr]

# add variable for window number
#chrom1[, window_50kb := 
#         cut(phys_pos, breaks = num.win, 
#             include.lowest = TRUE, labels = 1:num.win),
#       by = chr]
# get average by window    
#chrom1.50kb <- chrom1[, lapply(.SD, mean), 
#                      by = window_50kb, .SDcols = c("rec_rate", "pop_mean")]
#head(chrom1.50kb)

# plot window averages binned by recombination quantile
#chrom1.50kb[, .SD, by = cut(rec_rate, breaks=c(quantile(rec_rate, probs = seq(0, 1, by = 0.20))), 
#                            labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"), include.lowest = T)] %>%
#  ggplot(aes(x = cut, y = 1 - pop_mean)) + geom_jitter(width = 0.2, alpha = 0.5) +
#  labs(title = "Chromosome 1", x= "Recomb. rate quantile", y="Minor parent freq") +
#  stat_summary(fun.y=mean, geom="point", shape=18,
#               size=3, color="red") 


# 4. Interpolate ancestry on genetic scale ==========

# what unit to interpolate at?
# check # of SNPs / Morgan:
#length(unique(chrom1$phys_pos))/max(chrom1$Morgan) #roughly 50000 SNPs per M across chromosomes
# so we can interpolate to unit distance such that we get approx 1 SNP per unit
# this distance is roughly M*2^-15

xout = seq(0, max(chrom1$Morgan), by=2^-15) 

# interpolate individual ancestry at genetic coordinates
chrom1.interp.gen <- chrom1[, approx(x = Morgan, 
                                 y = ind_frq, 
                                 xout = xout, rule = 2), 
                        by = .(ID,chr)]
setnames(chrom1.interp.gen, c("x", "y"), c("Morgan","ind_frq_interp"))

# compute sample mean
chrom1.interp.gen[, smpl_mean_interp := mean(ind_frq_interp), 
              by = .(chr, Morgan)]


# 5. Interpolate ancestry on physical scale ==========

# approximately 1 SNP per kb so interpolate at this resolution

xoutbp = seq(0, max(chrom1$phys_pos), by=1e3) 

chrom1.interp.phys <- chrom1[, approx(x = phys_pos, 
                                 y = ind_frq, 
                                 xout = xoutbp, rule = 2), 
                        by = .(ID,chr)]
setnames(chrom1.interp.phys, c("x", "y"), c("position","ind_frq_interp"))

# compute sample mean
chrom1.interp.phys[, smpl_mean_interp := mean(ind_frq_interp), 
              by = .(chr, position)]

# interpolate recombination at same resolution
r.interp.phys <- chrom1[, approx(x = phys_pos, 
                y = r, 
                xout = xoutbp, rule = 2), by = .(ID,chr)]
setnames(r.interp.phys, c("x", "y"), c("position","r_interp"))

# combine these
chrom1.interp.phys <- merge(r.interp.phys,chrom1.interp.phys, by = c("ID","chr","position"))

# look at SNP density along genome
# u <- chrom1[,seq(min(phys_pos),max(phys_pos),length.out=512)]
# v <- density(unique(chrom1$phys_pos))$y
# plot(y=v, x = u)

# compare correlation at original data points to interpolated locations
# ggplot(chrom1[ID==ID[1]], aes(x = log(r) , y = 1-pop_mean)) + geom_point() + geom_smooth(method = "lm")
# cor.test(chrom1[ID==ID[1],log(r)], chrom1[ID==ID[1], 1-pop_mean])
# 
# ggplot(chrom1.interp.phys[ID==ID[1]], aes(x = log(r_interp),y = 1-pop_mean_interp)) + 
#    geom_point() + geom_smooth(method = "lm")
# cor.test(chrom1.interp.phys[ID==ID[1],log(r_interp)], chrom1.interp.phys[ID==ID[1],1-pop_mean_interp])

# output
save(chrom1, chrom1.interp.gen, chrom1.interp.phys, file = paste0("ACUA_",year,"/",scaff,".RData"))


# Examine effect of interpolation on genetic scale to correlation -----------

# intepolate rec map to get rec values at corresponding M distances
# r.chrom1.interp <- r.chrom1.allbp[, approx(x = Morgan, 
#                                  y = r, 
#                                  xout = xout, rule = 2)]

# setnames(r.chrom1.interp, c("x", "y"), c("Morgan","r"))

# add this to chrom1 interp
# chrom1.interp[,r := r.chrom1.interp$r, by = .(ID,chr)]

# ggplot(chrom1[ID==ID[1]], aes(x = log(r) , y = 1-pop_mean)) + geom_point() + geom_smooth(method = "lm")
# cor.test(chrom1[ID==ID[1],log(r)], chrom1[ID==ID[1], 1-pop_mean])

# ggplot(chrom1.interp[ID==ID[1]], aes(x = log(r),y = 1-pop_mean_interp)) + 
#   geom_point() + geom_smooth(method = "lm")
# cor.test(chrom1.interp[ID==ID[1],log(r)], chrom1.interp[ID==ID[1],1-pop_mean_interp])

# the correlation appears much weaker

# # Interpolation method #2 recombination rates 
# # for each SNP position, get cM per kb comparing to previous SNP
# cmPerKb <- function(x){
#   if ( x[,phys_pos] == min(chrom1[,phys_pos]) ){
#     lowerSNP <- 0
#     lowerM <- 0
#   } else{
#     lowerSNP <- max(chrom1[phys_pos < x[,phys_pos], phys_pos])
#     lowerM <- chrom1[phys_pos == lowerSNP, Morgan][1]
#   }
#   return(100*1000*(x[,phys_pos] - lowerSNP)/(x[,Morgan] - lowerM))
# }

# add identifier number to each snp
# chrom1[, num := seq_len(nrow(.SD)), by=ID]
# chrom1[ID==ID[1], cmPerKb(.SD),by=num]


# # Interpolation method #3 recombination rates 
# getRecRate <- function(x){
#   # x should be the interpolated genetic distance in chrom1.interp
#   # get the next smallest M distance and bp for real SNP locations 
#   if (min(chrom1[,Morgan]) < x[,Morgan]){
#     maxLess_M <- max(chrom1[Morgan < x[,Morgan], Morgan])
#     maxLess_Bp <- chrom1[Morgan == maxLess_M, phys_pos][1] 
#   } else{ 
#     maxLess_M <- 0
#     maxLess_Bp <- 0
#   }
#   
#   # next largest value of genetic distance
#   minGreater_M <- min(chrom1[Morgan >= x[,Morgan], Morgan])
#   minGreater_Bp <- chrom1[Morgan == minGreater_M, phys_pos][1] 
#   
#   # multiply by 1000 and 100 to get cM/kb from M/bp
#   return(1000*100*(minGreater_M-maxLess_M)/(minGreater_Bp-maxLess_Bp))
# }
#   


#chrom1.interp[ID == ID[1]] %>% ggplot(aes(x = gen_pos, y = 1-pop_mean_interp)) + geom_point()



