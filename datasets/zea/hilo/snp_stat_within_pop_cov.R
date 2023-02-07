library(data.table)
setDTthreads(threads = 3)
#source("~/workspace/gnomwav/R/multi_modwts.R")
#source("~/workspace/gnomwav/R/correlation_decomp.R")

source('/home/jgroh/gnomwav/R/multi_modwts.R')
source('/home/jgroh/gnomwav/R/correlation_decomp.R')

args <- commandArgs(trailingOnly = T)
maize_frqs <- fread(args[1])
mex_frqs <- fread(args[2])
rmap <- fread(args[3])
IDS <- fread(args[4], header=F, col.names='ID')


#maize_frqs <- fread("hilo/allopatric_maize.mafs.gz")
#mex_frqs <- fread("hilo/allopatric_mexicana.mafs.gz")
#rmap <- fread("hilo/ogut_2015_rmap_v2_to_v4_EXTENDED.txt")
#IDS <- fread('hilo/allopatric_mexicana_ids.list', header=F, col.names = 'ID')

pairs <- data.table(t(combn(IDS[,ID], 2))) #;pairs <- pairs[1:2]
setnames(pairs, c("ID1", "ID2"))

# function to do wavelet covariance for individuals
f <- function(x, y, dir){
  
  # read data for two individuals
  id1_frqs <- fread(paste0('individual_alleleFrqs_thinnedSNPs/', x, '.mafs.gz'))
  id2_frqs <- fread(paste0('individual_alleleFrqs_thinnedSNPs/', y, '.mafs.gz'))
  
  # merge with reference frqs
  d <- merge(merge(merge(maize_frqs[,.(chromo, position, maize_frq = phat)],
                         mex_frqs[,.(chromo, position, mex_frq = phat)]),
                   id1_frqs[,.(chromo, position, id1_frq = phat)]),
             id2_frqs[,.(chromo, position, id2_frq = phat)])

  d[, id1_snp_stat := (id1_frq-maize_frq)/(mex_frq-maize_frq)]
  d[, id2_snp_stat := (id2_frq-maize_frq)/(mex_frq-maize_frq)]
  
  # interpolate
  # get genetic position of informative SNPs
  snp_cm <- rmap[, approx(x=pos_bp, y=pos_cM, xout=d[chromo == .BY, position]), by = chr]
  setnames(snp_cm, c("x","y"),c("position","cm"))
  
  d <- merge(d, snp_cm[, .(chromo=chr,position,cm)])
  d[, Morgan := cm/100]
  
  interp1 <- d[, approx(x=Morgan, y=id1_snp_stat, xout=seq(min(Morgan), max(Morgan), by = 2^-14) ), by = chromo]
  setnames(interp1, c("x", "y"), c("Morgan", "id1_snp_stat"))
  
  interp2 <- d[, approx(x=Morgan, y=id2_snp_stat, xout=seq(min(Morgan), max(Morgan), by = 2^-14) ), by = chromo]
  setnames(interp2, c("x", "y"), c("Morgan", "id2_snp_stat"))
  interp <- merge(interp1, interp2, by = c("chromo", "Morgan"))
  
  wavcov <- interp[, cov_tbl(data=.SD, chromosome = 'chromo', signals = c('id1_snp_stat', 'id2_snp_stat'))]
  wavcov[, level := factor(level, levels = c(paste0('d', 1:15), 's14', 's15', 'chromo'))]
  wavcov[, ID1 := x][, ID2 := y][]
}

covs <- pairs[, f(x=ID1, y=ID2), by = seq_len(nrow(pairs))]
covs[, seq_len := NULL][]

fwrite(x=covs, file = "", quote = F, sep = "\t", row.names = F, col.names = T)
