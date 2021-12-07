library(data.table)
library(ggplot2)
library(magrittr)

# 1. Read Data Files ==========

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]
chr <- args[2]
#year <- 2018
#chr <- "ScyDAA6-10-HRSCAF-60"; chr <- "ScyDAA6-11-HRSCAF-73"; chr <- "ScyDAA6-1107-HRSCAF-1306"; chr <- "ScyDAA6-1196-HRSCAF-1406"; chr <- "ScyDAA6-

  
# read genotype files 
# probability homozygous for either ancestry
par1 <- fread(paste0("ancestry-probs-par1_allchrs_ACUA_historical_",year,".tsv"))
par2 <- fread(paste0("ancestry-probs-par2_allchrs_ACUA_historical_",year,".tsv"))

# read chrlengths
chrLen <- fread("xbir10x_chrlengths.txt", col.names = c("chromosome", "len"))

# read recomb map for focal chromosome 
rhoMap <- fread(paste0("LD_recMap/LD_map_xbirchmanni-COAC-10x-", chr,".post.txt_mod.bed"), col.names = c("chr","start","end","mean","V1","median","V3")) 
rhoMap <- rhoMap[, c("chr", "start", "end", "median")]
setnames(rhoMap, "median", "rho")


# 2. Reformat ancestry files =====

# merge by ID 
gnom <- merge(par1, par2, by = "V1", suffixes = c(":gen11", ":gen22"))
setnames(gnom, "V1", "ID")

# subset to focal chromosome
chromCols <- c("ID", names(gnom)[grep(paste0(chr,":"), names(gnom))])
chrom <- gnom[, ..chromCols]

# move genotype probs for all loci to a single column
chrom <- melt(chrom, id.vars = "ID", value.name = "postProb")

# separate by chromosome, position, genotype
chrom[, c("chr", "position", "genotype") := tstrsplit(variable,":",fixed=TRUE)]
chrom[,variable:=NULL]
chrom[, position := as.integer(position)]

# move genotypes to separate columns
chrom <- dcast(chrom, ... ~ genotype, value.var = "postProb")

# add column for frequency of allele in individual weighted by post. probs. 
chrom[, indivFreq := (gen11 + 0.5*(1 - sum(gen11, gen22))), by = .(ID,chr,position)]


# 3. Interpolate ancestry at evenly space physical distances =====

# midpoints of 1kb intervals
xout <- chrLen[chromosome==chr, seq(500, len, by = 500)]
xout <- xout[xout %% 1000 != 0]

chromAnc1kb <- chrom[, approx(x = position, y = indivFreq, xout=xout, rule=2), by = .(chr,ID)]
setnames(chromAnc1kb, c("x","y"), c("position", "indivFreq"))

# average over individuals to get sample mean
chromAnc1kb[, "meanFreq" := mean(indivFreq), by = .(chr, position)]

# plot thinned sample mean ancestry
# chromAnc1kb[position %% 555 == 0,] %>% ggplot(aes(x = position, y = meanFreq)) + geom_point()


# 4. Create recombination map ==========

if(min(rhoMap$start) != 0){
  rhoMap <- rbind(list(chr=rhoMap[1,chr],
                       start=0,
                       end=min(rhoMap$start),
                       rho=rhoMap[start==min(start),rho]),
                  rhoMap)
  }
if(max(rhoMap$end) < chrLen[chromosome==chr,len]){
  rhoMap <- rbind(rhoMap,
                  list(chr=rhoMap[1,chr],
                       start=max(rhoMap$end),
                       end=chrLen[chromosome==chr,len+1], # add +1 bc bed is zero-indexed
                       rho=rhoMap[end==max(end),rho]))
  }


Ne2 <- 27447 # see script xbir_makeBed_LDRecMap.R for Ne calculation

# cap outlier values 
rhoMap[rho >= 0.005, rho:= 0.005]

MorganVec <- cumsum(rhoMap[, rep(rho/Ne2, end-start)])


# 5. Interpolate ancestry at evenly space genetic distances =====

# assign genetic position of SNPS present in the data
chrom[, "Morgan" := MorganVec[position]]

# check number of SNPs per Morgan
# chrom[,length(unique(position))/max(Morgan)]
# roughly across chromosomes there are 30,000 SNPs per Morgan. So if we interpolate to M*2^-14, this gives roughly 1-2 SNPs per interpolation window

xoutMorgan <- seq(0,max(MorganVec), by = 2^-14)

# interpolate individual ancestry at genetic coordinates
chromAncInterpMorgan <- chrom[, approx(x = Morgan, 
                                     y = indivFreq, 
                                     xout = xoutMorgan, rule = 2), 
                            by = .(ID,chr)]
setnames(chromAncInterpMorgan, c("x", "y"), c("Morgan","indivFreq"))

# compute sample mean
chromAncInterpMorgan[, meanFreq := mean(indivFreq), 
                  by = .(chr, Morgan)]

# 5. Write Output =====

save(chromAnc1kb, chromAncInterpMorgan, file = paste0("ACUA_",year,"/",chr,".RData"))
