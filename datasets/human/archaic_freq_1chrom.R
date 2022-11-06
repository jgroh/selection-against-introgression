library(data.table)


if(Sys.getenv("RSTUDIO") == "1"){
  archaic_all <- fread("Neanderthal_files/41586_2020_2225_MOESM3_ESM.txt")
  rmap_all <- fread("aau1043_datas3")
  #rmap_all <- fread("recomb-hg38/genetic_map_GRCh38_merged.tab")
  chromosome <- 20
} else {
  args <- commandArgs(trailingOnly = TRUE)
  chromosome <- args[1]
  archaic_all <- fread(args[2])
  rmap_all <- fread(args[3])
}


# subset to focal chrom
archaic <- archaic_all[chrom == chromosome]
rmap <- rmap_all[Chr == paste0("chr", chromosome)]

archaic[, freq := freq/55132] # this is the number of haploid genomes used 
rmap[, Morgan := cM/100]


# ===== Allele frequency in at physical points =====
# get archaic allele frequency at 1kb intervals
pos1kb <- data.table(pos = seq(rmap[,min(End)], rmap[, max(End)], by = 1e3))

# count overlapping fragments at each position
frq1kb <- pos1kb[, freq := archaic[start < pos & end >= pos, sum(freq)], by = seq_len(nrow(pos1kb))][]
frq1kb[, chr := chromosome]
frq1kb[, rec := approx(xout = frq1kb$pos, x = rmap$End, y = rmap$cMperMb)$y]
#ggplot(frq1kb, aes(x = pos, y = freq)) + geom_point()



# ===== Allele frequency in genetic windows =====

# get physical coordinates of genetic windows 
# window size? windows should be roughly the same resolution as the recombination map, or more course
# this is between 2^-16 and 2^-17
# rmap[, mean(Morgan - shift(Morgan), na.rm=T)]
#x <- merge(rmap_all[, seq(min(Morgan), max(Morgan), by = 2^-16), by = chrom][, .(nw = .N), by = chrom], rmap_all[, .(nw2 = .N), by = chrom])
#ggplot(x, aes(x = nw, y = nw2)) + geom_point()

xout <- rmap[, seq(min(Morgan), max(Morgan), by = 2^-16)]

posM <- rmap[, approx(x = Morgan, y = End, xout = xout)]
setnames(posM, c("Morgan", "pos"))
posM[, pos := round(pos)]
ggplot(posM, aes(x = pos, y = 1/(pos-shift(pos)))) + geom_point()

# get overlapping fragments
frqM <- posM[, freq := archaic[start < pos & end >= pos, sum(freq)], by = seq_len(nrow(posM))][]
frqM[, chr := chromosome]

frqM[, rec1 := 1/(pos-shift(pos))]
frqM[, rec1 := c(rec[2], rec[2:nrow(.SD)])]
# gives same answer
#frqM[, rec2 := approx(xout = frqM$pos, x = rmap$End, y = rmap$cMperMb)$y] 
#plot(frqM$rec1, frqM$rec2) 
#abline(0,1)




# === output ===== 
fwrite(frq1kb, file = paste0("archaic_freqs/chr", chromosome, "_frq_physical_windows.txt"), quote = F, sep = "\t")
fwrite(frqM, file = paste0("archaic_freqs/chr", chromosome, "_frq_genetic_windows.txt"), quote = F, sep = "\t")


#ggplot(frq, aes(x = wstart, y = freq)) + geom_point()
