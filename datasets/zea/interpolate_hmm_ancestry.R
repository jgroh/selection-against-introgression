# ===== Input =====
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
indFile <- args[1]
ind <- args[2]
mapFile <- args[3]
chrLenFile <- args[4]
outPath <- args[5]

# # If running locally:
# indFile <- "HILO2.posterior"
# ind <- "HILO2"
# mapFile <- "ogut_2015_rmap_v2_to_v4_EXTENDED.txt"
# chrLenFile <- "Zea_mays.AFPv4.dna.chr.autosome.lengths"

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
#floor(-log2(max(map[, .(cm=max(pos_cM)-min(pos_cM)),by=chr][,cm/100]/gnom[, .(nSNP=length(unique(pos_bp))),by=chr][,nSNP])))
# interpolate to 2^-14 of a Morgan so that on average there is >1 SNP per interval

# get genetic position of SNPS
m <- merge(gnom,map, by = c("chr","pos_bp"),all=T)
m[, pos_cM := approx(y=pos_cM,x=pos_bp,xout=pos_bp)$y, by = chr]

# ----- logit transform ancestry -----
m[, freqMexTr := log(freqMex/(1-freqMex))]

# replace values of zero or 1 by small deviation so logit works
if(any(m$freqMex==0, na.rm=T) | any(m$freqMex ==1, na.rm=T)){
  epsilon <- m[freqMex > 0, min(freqMex)]/2
  m[freqMex == 0, freqMexTr := log(epsilon/(1-epsilon))]
  m[freqMex == 1, freqMexTr := log((1-epsilon)/epsilon)]
}

# interpolate on logit scale
MexAncGenScale <- m[, approx(x=pos_cM, 
           y=freqMexTr, 
           xout=seq(min(pos_cM),max(pos_cM)/100,by=2^-14),
           rule=2),
  by=chr]
setnames(MexAncGenScale, 'x', 'position') # note that units of 'position' are Morgans
MexAncGenScale[, freqMex := exp(y)/(1+exp(y))][, y := NULL] # inverse logit transform
MexAncGenScale[, Morgan := Morgan - min(Morgan), by = chr] # shift minimum to zero
MexAncGenScale[, ID := ind]

# ----- Physical Scale -----
# What scale to interpolate at?
chrLen[, len, by=chr]/gnom[,length(unique(pos_bp)),by=chr]
# roughly 1 SNP per 5-7 kb
# interpolate to 1kb for simplicity now

gnom <- merge(gnom,chrLen)

gnom[, freqMexTr := log(freqMex/(1-freqMex))]

# replace values of zero or 1 by small deviation so logit works
if(any(gnom$freqMex==0, na.rm=T) | any(gnom$freqMex ==1, na.rm=T)){
  epsilon <- gnom[freqMex > 0, min(freqMex)]/2
  gnom[freqMex == 0, freqMexTr := log(epsilon/(1-epsilon))]
  gnom[freqMex == 1, freqMexTr := log((1-epsilon)/epsilon)]
}


MexAncPhysScale <- gnom[, approx(x=pos_bp,
                                 y=freqMexTr,
                                 xout=seq(500,len[1],by=1000),
                                 rule=2),
                        by=chr]
setnames(MexAncPhysScale, 'x','position') # units are bps
MexAncPhysScale[, freqMex := exp(y)/(1+exp(y))][, y := NULL] # inverse logit transform
MexAncPhysScale[, ID := ind]

# ===== Output =====

write.table(MexAncGenScale,file = paste0(outPath,"/genetic/",ind,".txt"),quote=F,sep="\t",row.names=F)
write.table(MexAncPhysScale,file = paste0(outPath,"/physical/",ind,".txt"),quote=F,sep="\t",row.names=F)
