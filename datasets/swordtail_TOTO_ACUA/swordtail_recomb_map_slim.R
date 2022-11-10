library(data.table)
options(scipen=999)

# argument 1 is chromosome length file
# argument 2 is one of the files containing snp positions

args <- commandArgs(trailingOnly = TRUE)
chrLen <- fread(args[1], col.names = c("chr", "len")) # physical lengths
#chrLen <- fread("xbir10x_chrlengths.txt", col.names = c("chr", "len"))  

par1 <- fread(args[2], nrows=1)
#par1 <- fread("ancestry-probs-par1_allchrs_ACUA_historical_2006.tsv", nrows=1)

# merge by ID 
setnames(par1, "V1", "ID")

# concatenate SNP recombination distances for each chromosome

recVec <- NULL
chrVec <- NULL

for(chrz in chrLen$chr){
  chromCols <- names(par1)[grep(paste0(chrz,":"), names(par1))]
  
  l <- strsplit(chromCols, ":", fixed = TRUE)
  SNP_positions <- as.numeric(sapply(l, '[',2))
  
  rChrom <- fread(
    paste0("LD_recMap/", 
           "LD_map_xbirchmanni-COAC-10x-", 
           chrz, ".post.txt_mod.bed"), col.names = c("chr","start","end","mean_2Ner","V1","median_2Ner","V3"))
  
  rChrom <- rChrom[, c("chr", "start", "end", "median_2Ner")]
  
  if(min(rChrom$start) != 0){
    rChrom <- rbind(list(chr=rChrom[1,chr],
                         start=0,
                         end=min(rChrom$start),
                         median_2Ner=rChrom[start==min(start),median_2Ner]),
                    rChrom)
  }
  if(max(rChrom$end) < chrLen[chr==rChrom[1,chr],len]){
    rChrom <- rbind(rChrom,
                    list(chr=rChrom[1,chr],
                         start=max(rChrom$end),
                         end=chrLen[chr==rChrom[1,chr],len+1], # add +1 bc bed is zero-indexed
                         median_2Ner=rChrom[end==max(end),median_2Ner]))
  }
  
  # outliers: set max median_2Ner to 0.005 (see script xbir_makeBed_LDRecMap.R)
  rChrom[median_2Ner >= 0.005, median_2Ner:= 0.005]

  #MorganVec <- cumsum(rChrom[, rep(median_2Ner/75981, end-start)])
  MorganVec <- cumsum(rChrom[, rep(median_2Ner/27447, end-start)])
  
  
  MorganPositions <- MorganVec[SNP_positions]
  
  #print(c(chrz, length(SNP_positions)/max(MorganVec)))
  
  recVec <- c(recVec, 0.5, diff(MorganPositions))
  chrVec <- c(chrVec, rep(chrz, length(diff(MorganPositions)) + 1))
  #hist(log10(diff(MorganPositions)), xlab = "log10(Morgans) between adjacent SNPs")
}


fwrite(list(recVec[-1]), 
            file = "swordtail_SNP_recMap_slim_rhoCap0.005.txt")
fwrite(data.table(chrVec, recVec)[-1], 
       file = "swordtail_SNP_recMap_slim_rhoCap0.005_verbose.txt")

