library(data.table)
#library(ggplot2)

# 1. Read Data Files ==========

args <- commandArgs(trailingOnly = TRUE)
year <- args[1]
recfile <- args[2]
chr <- args[3]

# if running interactively:
#year <- 2018
#recfile <- "xbir_LDRecMap.bed"
#chr <- "ScyDAA6-10-HRSCAF-60"; chr <- "ScyDAA6-11-HRSCAF-73"; chr <- "ScyDAA6-1107-HRSCAF-1306"; chr <- "ScyDAA6-1196-HRSCAF-1406"; chr <- "ScyDAA6-

# read genotype files 
# probability homozygous for either ancestry
par1 <- fread(paste0("ancestry-probs-par1_allchrs_ACUA_historical_",year,".tsv"))
par2 <- fread(paste0("ancestry-probs-par2_allchrs_ACUA_historical_",year,".tsv"))

# read recomb map for focal chromosome 
recmap <- fread(recfile, col.names = c("chromosome", "start", "end", "region_id", "r"))
recmap <- recmap[chromosome == chr]


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
chrom[, variable := NULL]
chrom[, position := as.integer(position)]

# move genotypes to separate columns
chrom <- dcast(chrom, ... ~ genotype, value.var = "postProb")

# add column for frequency of allele in individual weighted by post. probs. 
chrom[, indivFreq := (gen11 + 0.5*(1 - sum(gen11, gen22))), by = .(ID,chr,position)]


# 3. Interpolate ancestry at evenly space physical distances =====

# midpoints of 1kb intervals. The first interval beginning at the first SNP position in the data
xout <-  seq(min(chrom$position) + 500, max(chrom$position), by = 1000)

# do interpolation
chromAnc1kb <- chrom[, approx(x = position, y = indivFreq, xout = xout, rule=1), by = .(chr, ID)]
setnames(chromAnc1kb, c("x","y"), c("position", "indivFreq"))

# average over individuals to get sample mean
chromAnc1kb[, "meanFreq" := mean(indivFreq), by = .(chr, position)]

# plot thinned sample mean ancestry
# chromAnc1kb[seq(1, nrow(chromAnc1kb), by = 1000)] %>% ggplot(aes(x = position, y = meanFreq)) + geom_point()


# 4. Interpolate ancestry at evenly spaced genetic distances =====

# assign genetic position of SNPS present in the data
MorganVec <- cumsum(recmap[, rep(r, end-start)])
chrom[, "Morgan" := MorganVec[position]]

# check number of SNPs per Morgan
# chrom[,length(unique(position))/max(Morgan)]
# roughly across chromosomes there are 30,000 SNPs per Morgan. So if we interpolate to M*2^-14, this gives roughly 1-2 SNPs per interpolation window

# midpoints of windows of length 2^-14 Morgans, left boundary of first window determined by leftmost SNP
xoutMorgan <- seq(min(chrom$Morgan) + 2^-15, max(chrom$Morgan), by = 2^-14)

# interpolate individual ancestry at genetic coordinates
chromAncInterpMorgan <- chrom[, approx(x = Morgan, 
                                     y = indivFreq, 
                                     xout = xoutMorgan, rule = 1), 
                            by = .(ID, chr)]
setnames(chromAncInterpMorgan, c("x", "y"), c("Morgan","indivFreq"))

# compute sample mean
chromAncInterpMorgan[, meanFreq := mean(indivFreq), 
                     by = .(chr, Morgan)]

# 5. ===== Calculate Recombination Rate of genetic windows 
# interpolate physical position at right endpoint of genetic windows
xoutMorganRight <- xoutMorgan + 2^-15


# If the endpoint of the final window is beyond the Morgan position of the final SNP, take the SNP's position as the endpoint
if( max(xoutMorganRight) > max(chrom$Morgan) ) {
  xoutMorganRight[length(xoutMorganRight)] <- max(chrom$Morgan)
}

chromBpInterp <- chrom[, approx(x = Morgan, 
                                       y = position, 
                                       xout = xoutMorganRight, rule = 1), 
                              by = .(ID, chr)]
setnames(chromBpInterp, c("x", "y"), c("Morgan","end"))
chromBpInterp[, end := ceiling(end)]

# get # number of bp between successive interpolation points
chromBpInterp[, bp_span := end - shift(end), by = ID]

# first value is NA from above operation, manually put in the value
chromBpInterp[, bp_span := c( min(end) - min(chrom$position), bp_span[-1]), by = ID]

# merge interpolated ancestry and bp size of windows
# adjust Morgan values back to center of windows for the sake of merging the tables together
chromBpInterp[, Morgan := Morgan - 2^-15]

# Since the final window is forced to end at the final SNP's genetic position, the above subtraction
# is not correct for the final window, so just set it directly. 
# technically, the ancestry estimate in the final window is not centered in that window,
# since that window ends at the last SNP position in the data, but who cares
chromBpInterp[Morgan == max(Morgan), Morgan := max(xoutMorgan)]

# the two tables should now have the same entries for the Morgan column, where the Morgan value
# is located in the center of the window (except for the very last window) and the bp_span column
# gives the measure of the recombination rate
chromAncInterpMorgan <- merge(chromBpInterp, chromAncInterpMorgan,  
      by = c("ID", "chr", "Morgan"))


# ===== make bed file of physical windows

if (year == 2018){
  # only need to do this for each chromosome once, so in one year
  # (these will need to be combined across chromosomes in a later step)
  
  physical_windows_bed <- chromAnc1kb[ID==ID[1], .(chr, position)]
  physical_windows_bed[, start := position - 500][, end := position + 500][, position := NULL][]
  
  if( max(physical_windows_bed$end) > max(chrom$position) ){
    physical_windows_bed[end == max(end), end :=  max(chrom$position) + 1]
  }
  
  write.table(physical_windows_bed, file = paste0(chr, "_physical_windows.bed"), quote = F, sep = "\t", col.names = F, row.names = F) 
  
  # genetic windows
  # recall chromBpInterp has the first row per ID with the endpoint being the starting position of the first SNP,
  # which we don't need
  genetic_windows_bed <- chromBpInterp[ID==ID[1], .(chr = chr, start = end - bp_span, end = end)]
  genetic_windows_bed[nrow(genetic_windows_bed), end := end + 1]
  
  write.table(genetic_windows_bed, file = paste0(chr, "_genetic_windows.bed"), quote = F, sep = "\t", col.names = F, row.names = F) 
  
}


# 5. Write out interpolated files =====

save(chromAnc1kb, chromAncInterpMorgan, file = paste0("ACUA_", year, "/", chr, ".RData"))
