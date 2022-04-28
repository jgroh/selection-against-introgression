library(data.table)

args <- commandArgs(trailingOnly = TRUE)
chromosome <- args[1]
archaic_all <- fread(args[2])
rmap_all <- fread(args[3])

#archaic_all <- fread("Neanderthal_files/41586_2020_2225_MOESM3_ESM.txt")

#rmap_all <- fread("recomb-hg38/genetic_map_GRCh38_merged.tab")
#rmap_all[, max(pos_cm)/100, by = chrom]

#chromosome <- 20

# subset to focal chrom
archaic <- archaic_all[chrom == chromosome]
rmap <- rmap_all[chrom == paste0("chr", chromosome)]

#rmap[, my_rate := (pos_cm-shift(pos_cm))/(pos-shift(pos))]
#ggplot(rmap, aes(x = recomb_rate, y = my_rate)) + geom_point() + geom_pointdensity()
archaic[, freq := freq/55132]

#rmap[, Morgan_width := (stdrate*0.0116)/100]

rmap[, Morgan := pos_cm/100]


# ===== Functions for getting allele frequency =====

# not giving expected result atm 
# frac_overlap <- function(x.start, x.end, y.start, y.end){
#   # returns fraction of overlap of window covered by an archaic block
#   x.len <- abs(x.end - x.start)
#   max.start = max(c(x.start, y.start))
#   min.end = min(c(x.end, y.end))
#   overlap = min.end - max.start
#   overlap = ifelse(overlap <= 0, 0, overlap)
#   return(overlap / x.len)
# }
# 
# window_freq <- function(wstart, wend){
#   # gets overlapping archaic blocks, returns average over blocks of frequency weighted by overlap
#   fragments <- archaic[(start > wstart & start < wend) | (end > wstart & end < wend)]
#   if(nrow(fragments) == 0){
#     return(0)
#   } else{
#     return(fragments[, freq*frac_overlap(wstart, wend, start, end), by = seq_len(nrow(fragments)) ][, sum(V1)] )
#   }
# }


# ===== Allele frequency in at physical points =====
# what's the resolution of the recomb map? < 1 kb
#rmap_all[, pos-shift(pos), by = chrom][, mean(V1, na.rm=T)]
#rmap_all[, pos-shift(pos), by = chrom][, median(V1, na.rm=T)]


# get archaic allele frequency at these locations
pos1kb <- data.table(pos = seq(rmap[,min(pos)], rmap[, max(pos)], by = 1e3))

#frq1kb <- windows1kb[, .(freq = window_freq(wstart, wend)), by = wstart]
#ggplot(frq1kb, aes(x = wstart, y = freq)) + geom_point()

frq1kb <- pos1kb[, freq := archaic[start < pos & end >= pos, sum(freq)], by = seq_len(nrow(pos1kb))][]
frq1kb[, chr := chromosome]
frq1kb[, rec := approx(xout = frq1kb$pos, x = rmap$pos, y = rmap$recomb_rate)$y]
#ggplot(frq1kb, aes(x = pos, y = freq)) + geom_point()



# ===== Allele frequency in genetic windows =====

# get physical coordinates of genetic windows 
# window size? windows should be roughly the same resolution as the recombination map, or more course
# this is between 2^-16 and 2^-17
# rmap[, mean(Morgan - shift(Morgan), na.rm=T)]
#x <- merge(rmap_all[, seq(min(Morgan), max(Morgan), by = 2^-16), by = chrom][, .(nw = .N), by = chrom], rmap_all[, .(nw2 = .N), by = chrom])
#ggplot(x, aes(x = nw, y = nw2)) + geom_point()

xout <- rmap[, seq(min(Morgan), max(Morgan), by = 2^-16)]

posM <- rmap[, approx(x = Morgan, y = pos, xout = xout)]
setnames(posM, c("Morgan", "pos"))
posM[, pos := round(pos)]

frqM <- posM[, freq := archaic[start < pos & end >= pos, sum(freq)], by = seq_len(nrow(posM))][]
frqM[, chr := chromosome]
frqM[, rec := approx(xout = frqM$pos, x = rmap[, pos], y = rmap[, recomb_rate])$y]
#ggplot(frqM, aes(x = pos, y = freq)) + geom_point()




# === output ===== 
fwrite(frq1kb, file = paste0("archaic_freqs/chr", chromosome, "_frq_1kb_windows.txt"), quote = F, sep = "\t")
fwrite(frqM, file = paste0("archaic_freqs/chr", chromosome, "_frq_genetic_windows.txt"), quote = F, sep = "\t")


#ggplot(frq, aes(x = wstart, y = freq)) + geom_point()
