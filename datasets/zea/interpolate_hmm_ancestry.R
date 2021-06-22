# ===== Input =====
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
indFile <- args[1]
ind <- args[2]
mapFile <- args[3]
chrLenFile <- args[4]
outPath <- args[5]

# If running locally:
#indFile <- "HILO2.posterior"
#ind <- "HILO2"
#mapFile <- "ogut_2015_rmap_v2_to_v4_EXTENDED.txt"
#chrLenFile <- "Zea_mays.AFPv4.dna.chr.autosome.lengths"

# hmm post. prob. for individual 
gnom <- fread(indFile)
setnames(gnom, c('2,0','1,1','0,2','position','chrom'),
         c('p2.0','p1.1','p0.2','pos_bp','chr'))
gnom[, freqMex := 0.5*p1.1 + p0.2]

# recombination map
map <- fread(mapFile)

# chromosome physical lengths
chrLen <- fread(chrLenFile, col.names = c("chr","len"))

# ===== Interpolate Ancestry  =====

# ----- Genetic scale -----

# What scale to interpolate at?
#-log2(map[, max(pos_cM)-min(pos_cM),by=chr][,V1/100]/gnom[,length(unique(pos_bp)),by=chr][,V1])
# interpolate to 2^-14 of a Morgan so that on average there is >1 SNP per interval

m <- merge(gnom,map, by = c("chr","pos_bp"),all=T)
m[, pos_cM := approx(y=pos_cM,x=pos_bp,xout=pos_bp)$y, by = chr]

# get genetic position of SNPS
MexAncGenScale <- m[, approx(x=pos_cM, 
           y=freqMex, 
           xout=seq(min(pos_cM),max(pos_cM)/100,by=2^-14),
           rule=2),
  by=chr]
setnames(MexAncGenScale, c('x','y'),c('Morgan','freqMex'))
MexAncGenScale[, Morgan := Morgan - min(Morgan), by = chr] # shift minimum to zero
MexAncGenScale[, ID := ind]

# ----- Physical Scale -----
# What scale to interpolate at?

#chrLen[, len, by=chr]/gnom[,length(unique(pos_bp)),by=chr]
# roughly 1 SNP per 5-7 kb
# interpolate to 1kb 

gnom <- merge(gnom,chrLen)
MexAncPhysScale <- gnom[, approx(x=pos_bp,
                                 y=freqMex,
                                 xout=seq(500,len[1],by=1000),
                                 rule=2),
                        by=chr]
setnames(MexAncPhysScale, c('x','y'),c('position','freqMex'))
MexAncPhysScale[, ID := ind]

# ===== Output =====

write.table(MexAncGenScale,file = paste0(outPath,"/genetic/",ind,".txt"),quote=F,sep="\t",row.names=F)
write.table(MexAncPhysScale,file = paste0(outPath,"/physical/",ind,".txt"),quote=F,sep="\t",row.names=F)
