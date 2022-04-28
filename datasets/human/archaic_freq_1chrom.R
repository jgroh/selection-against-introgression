library(data.table)

args <- commandArgs(trailingOnly = TRUE)
chromosome <- args[1]
archaic_all <- fread(args[2])
rmap_all <- fread(args[3])

#archaic_all <- fread("Neanderthal_files/41586_2020_2225_MOESM3_ESM.txt")
#rmap_all <- fread("Kong_etal_recmaps/sex-averaged.rmap")
#chromosome <- 20

# subset to focal chrom
archaic <- archaic_all[chrom == chromosome]
rmap <- rmap_all[chr == paste0("chr", chromosome)]

archaic[, freq := freq/55132]

# convert to Morgans. using the average genetic length of 10kb bins cited in Kong et al averaged between sexes
rmap[, Morgan_width := (stdrate*0.0116)/100]
rmap[, Morgan_end := cumsum(Morgan_width)]


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


# ===== Allele frequency in physical windows =====
# do 10kb windows for now since this is the resolution of the recombination map 


# get archaic allele frequency at these locations
windows10kb <- data.table(wstart = rmap[,pos] - 5e3, wmid = rmap[, pos], wend = rmap[, pos + 5e3])
windows10kb[, rec := rmap[, stdrate]]

#frq1kb <- windows1kb[, .(freq = window_freq(wstart, wend)), by = wstart]
#ggplot(frq1kb, aes(x = wstart, y = freq)) + geom_point()

frq10kb <- windows10kb[, freq := archaic[start < wmid & end >= wmid, sum(freq)], by = seq_len(nrow(windows10kb))][]
frq10kb[, chr := chromosome]
#ggplot(frq10kb, aes(x = wmid, y = freq)) + geom_point()



# ===== Allele frequency in genetic windows =====

# get physical coordinates of genetic windows 
# window size? windows should be roughly the same resolution as the recombination map, or more course
# rmap[, mean(Morgan_width)] # this is between 2^-13 and 2^-12


# we'll treat the *endpoints of the physical windows as the *midpoints of the genetic windows
# (because we can associate the endpoints of the physical windows with their Morgan distance, given that we have window lengths in Morgans)

# desired Morgan positions - start and endpoints of windows
xout_end <- rmap[, seq(min(Morgan_end) + 0.5*2^-12, max(Morgan_end) + 0.5*2^-12, by = 2^-12)]
xout_mid <- xout_end - 0.5*2^-12
xout_start <- xout_end - 2^-12

# get corresponding bp positions
midpoint <- rmap[, approx(xout = xout_mid, x = Morgan_end, y = pos+5e3)] 
setnames(midpoint, c("x", "y"), c("Morgan_midpoint", "wmid"))

# since pos were given as midpoints of intervals, add 5kb
start <- rmap[, approx(xout = xout_start, x = Morgan_end, y = pos+5e3)] # throws a warning due to windows of 0 Morgans, but result is as desired
start[, x := x + 0.5*2^-12] # add 1/2 interval length here so that x corresponds to Morgan midpoint of interval
setnames(start, c("x", "y"), c("Morgan_midpoint", "wstart"))

end <- rmap[, approx(xout = xout_end, x = Morgan_end, y = pos+5e3)] # throws a warning due to windows of 0 Morgans, but result is as desired
end[, x := x - 0.5*2^-12]
setnames(end, c("x", "y"), c("Morgan_midpoint", "wend"))


windowsM <- merge(midpoint, merge(start, end, by = "Morgan_midpoint"), by = "Morgan_midpoint")
windowsM <- windowsM[!is.na(wstart) & !is.na(wend)]

# get frequency at midpoint
frqM <- windowsM[, freq := archaic[start < wmid & end >= wmid, sum(freq)], by = seq_len(nrow(windowsM))][]
#ggplot(frqM, aes(x = wmid, y = freq)) + geom_point()

frqM[, rec := 1/(wend - wstart)]
frqM[, chr := chromosome]
# ggplot(frqM, aes(x = wmid, y = rec)) + geom_point()

# === output ===== 
fwrite(frq10kb, file = paste0("chr", chromosome, "_frq_10kb_windows.txt"), quote = F, sep = "\t")
fwrite(frqM, file = paste0("chr", chromosome, "_frq_genetic_windows.txt"), quote = F, sep = "\t")


#ggplot(frq, aes(x = wstart, y = freq)) + geom_point()
