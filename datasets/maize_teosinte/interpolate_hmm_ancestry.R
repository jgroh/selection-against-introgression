library(data.table)

# ===== Read Data =====
args <- commandArgs(trailingOnly = TRUE)
indFile <- args[1]
ind <- args[2]
mapFile <- args[3]

# hmm post. prob. for individual 
gnom <- fread(indFile)
setnames(gnom, c('2,0','1,1','0,2','position','chrom'),
         c('p2.0','p1.1','p0.2','pos_bp','chr'))
gnom[, freqMex := 0.5*p1.1 + p0.2]

# recombination map
map <- fread(mapFile)


# ===== Interpolate Ancestry in Genetic Distance =====

# What scale to interpolate at?
#-log2(map[, max(pos_cM)-min(pos_cM),by=chr][,V1/100]/gnom[,length(unique(pos_bp)),by=chr][,V1])
# interpolate to 2^-14 of a Morgan so that on average there is >1 SNP per interval

m <- merge(gnom,map, by = c("chr","pos_bp"),all=T)
m[, pos_cM := approx(y=pos_cM,x=pos_bp,xout=pos_bp)$y, by = chr]


MexAncGenScale <- m[, approx(x=pos_cM, 
           y=freqMex, 
           xout=seq(0,max(pos_cM)/100,by=2^-14),
           rule=2),
  by=chr]

setnames(MexAncGenScale, c('x','y'),c('Morgan','freqMex'))
MexAncGenScale[, ID := ind]

# write output to console
write.table(MexAncGenScale,file = "",quote=F,sep="\t",row.names=F)





