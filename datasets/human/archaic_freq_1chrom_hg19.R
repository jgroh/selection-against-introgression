library(data.table)


if(Sys.getenv("RSTUDIO") == "1"){
  chromosome <- 4
  skov_fragments <- fread("Skov_etal_2020_data/41586_2020_2225_MOESM3_ESM.txt")
  skov_loc <- fread("Skov_etal_2020_data/Skov_fragments_hg19.bed", col.names = c("chr", "start", "end", "hg38_loc"))
  sank_calls <- suppressWarnings(fread("Sankararaman_etal_2014_data/chr-4.thresh-90.length-0.00.gz"))
  stein_calls <- fread("Steinrucken_etal_2018/March2018/CEU_lax_chr4/chr4_frqs.txt", col.names = c("chr", "pos", "freq"))
} else {
  args <- commandArgs(trailingOnly = TRUE)
  chromosome <- args[1]
  sank_calls <- fread(args[2])
  stein_calls <- fread(args[3], col.names = c("chr", "pos", "freq"))
  skov_fragments <- fread(args[4])
  skov_loc <- fread(args[5], col.names = c("chr", "start", "end", "hg38_loc"))
}

# ===== formatting 

# column names of sank_calls
# col2: chromosome, col3: genetic position, col4: physical position, col11: average posterior probability of N
setnames(sank_calls, paste0("V", 1:17))
setnames(sank_calls, c("V2", "V3","V4", "V11"), c("chr", "cM", "pos", "freq"))

# icelander neanderthal fragment locations in hg19 coordinates
skov_fragments[, hg38_loc := paste0(chrom, "_", start, "_", end, "_", .I)]
skov_fragments <- merge(skov_fragments[chrom == chromosome, .(hg38_loc, freq)], skov_loc[, .(start, end, hg38_loc)], by = "hg38_loc")
skov_fragments[, freq := freq/55132] 
skov_fragments[, chrom := as.character(chromosome)]

# ===== Neanderthal frequency in physical windows
windows <- data.table(chrom=as.character(chromosome), pos = seq(min(sank_calls$pos), max(sank_calls$pos), by = 1e3))
windows[, start := pos - 500][, end := pos + 500]

# ---- skov estimates
# find all fragments that overlap grid windows
setkey(skov_fragments, chrom, start, end)
overlaps <- foverlaps(windows, skov_fragments, by.x=c("chrom","start","end"), by.y=c("chrom","start","end"))

# weight frequency by length of fragment overlap
overlaps[, weighted_freq := freq*(pmin(end, i.end) - pmax(start, i.start))/1e3 , by=.(chrom,start,end)]

# windows with no overlap get frequency zero
overlaps[is.na(freq), weighted_freq := 0]

# sum weighted frequencies of fragments in each window
frq1kb <- overlaps[, .(skov_freq = sum(weighted_freq)), by = .(chrom, i.start, pos, i.end)]
setnames(frq1kb, c("i.start", "i.end"), c("start", "end"))


# ---- add Sankararaman and Steinrucken estimates 
frq1kb[, sank_freq := approx(x = sank_calls$pos, y = sank_calls$freq, xout = pos)$y, by = chrom]
frq1kb[, stein_freq := approx(x = stein_calls$pos, y = stein_calls$freq, xout = pos)$y, by = chrom]


# ------ add recomb info 
frq1kb[, Morgan := approx(x = sank_calls$pos, y = sank_calls$cM, xout = pos)$y/100, by = chrom]

frq1kb[, Morgan_start := approx(xout = pos-500, y = sank_calls$cM/100, x = sank_calls$pos)$y]
frq1kb[, Morgan_end := approx(xout = pos+500, y = sank_calls$cM/100, x = sank_calls$pos)$y]

frq1kb <- frq1kb[!is.na(Morgan_start) & !is.na(Morgan_end)]
frq1kb[, rec := (Morgan_end-Morgan_start)/1e3]

frq1kb[, start := NULL][, end := NULL]


# ===== Neanderthal frequency in genetic windows

windowsM <- data.table(chrom = as.character(chromosome), Morgan = sank_calls[, seq(min(cM)/100, max(cM)/100, by = 2^-16)])
windowsM[, 
         start := round(approx(x = sank_calls$cM/100, y = sank_calls$pos, xout = windowsM$Morgan - 2^-17)$y)][
           , pos := round(approx(x = sank_calls$cM/100, y = sank_calls$pos, xout = windowsM$Morgan)$y)][
             , end := round(approx(x = sank_calls$cM/100, y = sank_calls$pos, xout = windowsM$Morgan + 2^-17)$y)
           ]
windowsM <- windowsM[!is.na(start) & !is.na(end)]

# ---- Skov et al data

# find all fragments that overlap grid windows
overlapsM <- foverlaps(windowsM, skov_fragments, by.x=c("chrom","start","end"), by.y=c("chrom","start","end"))

# weight frequency by length of fragment overlap
overlapsM[, weighted_freq := freq*(pmin(end, i.end) - pmax(start, i.start))/(i.end-i.start), by=.(chrom,i.start,i.end)]

# windows with no overlap get frequency zero
overlapsM[is.na(freq), weighted_freq := 0]

# sum weighted frequencies of fragments in each window
frqM <- overlapsM[, .(skov_freq = sum(weighted_freq)), by = .(chrom, Morgan, i.start, pos, i.end)]
setnames(frqM, c("i.start", "i.end"), c("start", "end"))

frqM[, rec := 2^-16/(end-start)]

# ----- add Sankararaman and Steinrucken estimates

frqM[, sank_freq := approx(x=sank_calls$pos, y = sank_calls$freq, xout=pos)$y]
frqM[, stein_freq := approx(x=stein_calls$pos, y = stein_calls$freq, xout=pos)$y]


# === output ===== 
fwrite(frq1kb, file = paste0("archaic_freqs_hg19/chr", chromosome, "_frq_physical_windows.txt"), quote = F, sep = "\t")
fwrite(frqM, file = paste0("archaic_freqs_hg19/chr", chromosome, "_frq_genetic_windows.txt"), quote = F, sep = "\t")


