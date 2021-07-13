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
gnom[, freq_sp2 := 0.5*p1.1 + p0.2]

# recombination map
map <- fread(mapFile)

# chromosome physical lengths
chrLen <- fread(chrLenFile, col.names = c("chr","len"))


# ===== Interpolate Ancestry  =====

# ----- Genetic scale -----

# What scale to interpolate at?
# floor(-log2(max(map[, .(cm=max(pos_cM)-min(pos_cM)),by=chr][,cm/100]/gnom[, .(nSNP=length(unique(pos_bp))),by=chr][,nSNP])))
# interpolate to 2^-14 of a Morgan so that on average there is >1 SNP per interval

# get genetic position of SNPS
m <- merge(gnom, map, by = c("chr","pos_bp"),all=T)
m[, pos_cM := approx(y=pos_cM,x=pos_bp,xout=pos_bp)$y, by = chr]

# ----- logit transform ancestry -----
m[, freq_sp2_tr := log(freq_sp2/(1-freq_sp2))]

# replace values of zero or 1 by small deviation so logit works
if(any(m$freq_sp2==0, na.rm=T) | any(m$freq_sp2 ==1, na.rm=T)){
  epsilon <- m[freq_sp2 > 0, min(freq_sp2)]/2
  m[freq_sp2 == 0, freq_sp2_tr := log(epsilon/(1-epsilon))]
  m[freq_sp2 == 1, freq_sp2_tr := log((1-epsilon)/epsilon)]
}

# interpolate on logit scale
sp2AncGenScale <- m[, approx(x=pos_cM, 
           y=freq_sp2_tr, 
           xout=seq(min(pos_cM)/100,max(pos_cM)/100,by=2^-14),
           rule=2),
  by=chr]
setnames(sp2AncGenScale, 'x', 'position') # note that units of 'position' are Morgans
sp2AncGenScale[, freq_sp2 := exp(y)/(1+exp(y))][, y := NULL] # inverse logit transform
sp2AncGenScale[, position := position - min(position), by = chr] # shift minimum to zero
sp2AncGenScale[, ID := ind]


# ----- Physical Scale -----
# What scale to interpolate at?
chrLen[, len, by=chr]/gnom[,length(unique(pos_bp)),by=chr]
# roughly 1 SNP per 5-7 kb
# interpolate to 5kb

m <- merge(m,chrLen, by = "chr", all=T) # add chrom length as we'll seq to the end of the chrom

sp2AncPhysScale <- m[, approx(x=pos_bp,
                                 y=freq_sp2_tr,
                                 xout=seq(2500,len[1],by=5000),
                                 rule=2),
                        by=chr]
setnames(sp2AncPhysScale, 'x','position') # units are bps
sp2AncPhysScale[, freq_sp2 := exp(y)/(1+exp(y))][, y := NULL] # inverse logit transform
sp2AncPhysScale[, ID := ind]


# ===== Output =====

write.table(sp2AncGenScale, file = paste0(outPath,"/genetic/",ind,".txt"), quote=F,sep="\t",row.names=F)
write.table(sp2AncPhysScale, file = paste0(outPath,"/physical/",ind,".txt"), quote=F,sep="\t",row.names=F)

