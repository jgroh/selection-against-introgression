library(data.table)


if(Sys.getenv("RSTUDIO") == "1"){
  chromosome <- 22
  skov_fragments <- fread("Skov_etal_2020_data/41586_2020_2225_MOESM3_ESM.txt")
  sank_calls <- fread("Sankararaman_etal_2014_data/chr-22.thresh-90.length-0.00.gz")
  sank_snps <- fread("Sankararaman_etal_2014_data/sankHg19ToHg38Snps.bed", col.names = c("chr", "start", "end", "hg19_id"))
  stein_calls <- fread("Steinrucken_etal_2018/March2018/CEU_lax_chr22/chr22_frqs.txt", col.names = c("chr", "pos", "freq"))
  stein_snps <- fread("Steinrucken_etal_2018/March2018/steinHg19ToHg38Snps.bed", col.names = c("chr", "start", "end", "hg19_id"))
  rmap_all <- fread("Halldorsson_etal_2019_data/aau1043_datas3")
} else {
  args <- commandArgs(trailingOnly = TRUE)
  chromosome <- args[1]
  skov_fragments <- fread(args[2])
  sank_calls <- fread(args[3])
  sank_snps <- fread(args[4], col.names = c("chr", "start", "end", "hg19_id"))
  stein_calls <- fread(args[5], col.names = c("chr", "pos", "freq"))
  stein_snps <- fread(args[6], col.names = c("chr", "start", "end", "hg19_id"))
  rmap_all <- fread(args[7])
}

# subset skov fragments and recomb map to focal chrom
archaic <- skov_fragments[chrom == chromosome, .(chrom, start, end, freq)]
rmap <- rmap_all[Chr == paste0("chr", chromosome)]

# convert count to frequency, this is the number of haploid genomes used 
archaic[, freq := freq/55132] 

# Morgans
rmap[, Morgan := cM/100]



# ===== Skov frequency physical windows  =====

# ----- method 1, weighted average of fragment frequencies -----
windows <- data.table(chrom=chromosome, pos = seq(rmap[,min(End)], rmap[, max(End)], by = 1e3))
windows[, start := pos - 500][, end := pos + 500]


# find all fragments that overlap grid windows
setkey(archaic, chrom, start, end)
overlaps <- foverlaps(windows, archaic, by.x=c("chrom","start","end"), by.y=c("chrom","start","end"))

# weight frequency by length of fragment overlap
overlaps[, weighted_freq := freq*(pmin(end, i.end) - pmax(start, i.start))/1e3 , by=.(chrom,start,end)]

# windows with no overlap get frequency zero
overlaps[is.na(freq), weighted_freq := 0]

# sum weighted frequencies of fragments in each window
frq1kb <- overlaps[, .(skov_freq = sum(weighted_freq)), by = .(chrom, i.start, pos, i.end)]
setnames(frq1kb, c("i.start", "i.end"), c("start", "end"))

frq1kb[, Morgan_end := approx(xout = frq1kb$end, x = rmap$End, y = rmap$Morgan)$y]
frq1kb[, Morgan_start := approx(xout = frq1kb$start, x = rmap$End, y = rmap$Morgan)$y]

frq1kb[, rec := (Morgan_end-Morgan_start)/1e3]
frq1kb[, rec := c(.SD[2, rec], .SD[2:nrow(.SD), rec])]

frq1kb[, Morgan_end := NULL][, Morgan_start := NULL][, start := NULL][, end := NULL]

# ---- method 2: count overlapping fragments at bp in center of window
# frq1kb[, skov_freq2 := archaic[start < pos & end >= pos, sum(freq)], by = seq_len(nrow(frq1kb))][]


# ===== Skov allele frequency in genetic windows =====

# get physical coordinates of genetic windows 
# window size? windows should be roughly the same resolution as the recombination map, or more course
# this is between 2^-16 and 2^-17
# rmap[, mean(Morgan - shift(Morgan), na.rm=T)]

windowsM <- data.table(chrom = chromosome, Morgan = rmap[, seq(min(Morgan), max(Morgan), by = 2^-16)])
windowsM[, 
         start := round(approx(x = rmap$Morgan, y = rmap$End, xout = windowsM$Morgan - 2^-17)$y)][
           , pos := round(approx(x = rmap$Morgan, y = rmap$End, xout = windowsM$Morgan)$y)][
             , end := round(approx(x = rmap$Morgan, y = rmap$End, xout = windowsM$Morgan + 2^-17)$y)
           ]
windowsM <- windowsM[!is.na(start) & !is.na(end)]


# find all fragments that overlap grid windows
overlapsM <- foverlaps(windowsM, archaic, by.x=c("chrom","start","end"), by.y=c("chrom","start","end"))

# weight frequency by length of fragment overlap
overlapsM[, weighted_freq := freq*(pmin(end, i.end) - pmax(start, i.start))/(end-start), by=.(chrom,start,end)]

# windows with no overlap get frequency zero
overlapsM[is.na(freq), weighted_freq := 0]

#ggplot(overlapsM, aes(x = freq, y = weighted_freq)) + geom_point() + geom_abline()


# sum weighted frequencies of fragments in each window
frqM <- overlapsM[, .(skov_freq = sum(weighted_freq)), by = .(chrom, Morgan, i.start, pos, i.end)]
setnames(frqM, c("i.start", "i.end"), c("start", "end"))

frqM[, rec := 2^-16/(end-start)]


# ggplot(frq1kb, aes(x = pos)) + geom_point(aes(y = skov_freq)) + 
#  geom_point(data=frqM, aes(x = end, y = skov_freq, color = 'red'))


# ===== Frequency estimates from Sankararaman et al 2014 =====
sank <- merge(sank_snps, 
      sank_calls[, hg19_id := paste0("chr", V2, "_", V4)][, .(sank_freq = V11, hg19_id)], 
      by = "hg19_id")
sank[, start := NULL]
setnames(sank, "end", "pos")

frq1kb[, sank_freq := approx(x = sank$pos, y = sank$sank_freq, xout = frq1kb$pos)$y]

#ggplot(frq1kb, aes(x = pos, y=skov_freq)) + geom_point() + geom_point(aes(x = pos, y= sank_freq, color = 'red'))

frqM[, sank_freq := approx(x = sank$pos, y = sank$sank_freq, xout = frqM$pos)$y]

#ggplot(frqM, aes(x = pos, y=skov_freq)) + geom_point() + geom_point(aes(x = pos, y= sank_freq, color = 'red', alpha = 0.5))


# ===== Frequency estimates from Steinrucken et al 2014 =====
stein_calls[, hg19_id := paste0("chr", chr, "_", pos)]
stein_snps[, start := NULL]
setnames(stein_snps, "end", "pos")

stein <- merge(stein_snps, stein_calls[, .(freq, hg19_id)], by = 'hg19_id')

frq1kb[, stein_freq := approx(x = stein$pos, y = stein$freq, xout = frq1kb$pos)$y]
frqM[, stein_freq := approx(x = stein$pos, y = stein$freq, xout = frqM$pos)$y]

#ggplot(frq1kb, aes(x = pos, y = sank_freq)) + geom_point() + geom_point(aes(y = stein_freq), color = 'blue')
#ggplot(frq1kb, aes(x = stein_freq, y = sank_freq)) + geom_point()


# === output ===== 
fwrite(frq1kb, file = paste0("archaic_freqs/chr", chromosome, "_frq_physical_windows.txt"), quote = F, sep = "\t")
fwrite(frqM, file = paste0("archaic_freqs/chr", chromosome, "_frq_genetic_windows.txt"), quote = F, sep = "\t")


#ggplot(frq, aes(x = wstart, y = freq)) + geom_point()
