library(data.table)

args <- commandArgs(trailingOnly = TRUE)
chromosome <- args[1]
archaic_all <- fread(args[2])
rmap_all <- fread(args[3])

#archaic_all <- fread("Neanderthal_files/41586_2020_2225_MOESM3_ESM.txt")
#head(moes_archaic)
#rmap_all <- fread("Kong_etal_recmaps/sex-averaged.rmap")

# subset to focal chrom
archaic <- archaic_all[chrom == chromosome]
rmap <- rmap_all[chr == paste0("chr", chromosome)]

archaic[, freq := freq/55132]

# all fragments
#sum(moes_archaic$freq)

# convert to Morgans. using the average genetic length of 10kb bins cited in Kong et al averaged between sexes
rmap[, Morgan_width := (stdrate*0.0116)/100]
rmap[, Morgan_end := cumsum(Morgan_width)]


# get physical coordinates of genetic windows (need to pick sensible window size)

# we'll treat the *endpoints of the physical windows as the *midpoints of the genetic windows
# (because we can associate the endpoints of the physical windows with their Morgan distance, given that we have window lengths in Morgans)

# desired Morgan positions - start and endpoints of windows
xout_start <- rmap[, seq(min(Morgan_end) - 0.5*2^-15, max(Morgan_end) - 0.5*2^-15, by = 2^-15)]
xout_end <- rmap[, seq(min(Morgan_end) + 0.5*2^-15, max(Morgan_end) + 0.5*2^-15, by = 2^-15)]


# get corresponding bp positions
# since pos were given as midpoints of intervals, add 5kb
start <- rmap[, approx(xout = xout_start, x = Morgan_end, y = pos+5e3)] # throws a warning due to windows of 0 Morgans, but result is as desired
start[, x := x + 0.5*2^-15] # add 1/2 interval length here so that x corresponds to Morgan midpoint of interval
setnames(start, c("x", "y"), c("Morgan_midpoint", "wstart"))

end <- rmap[, approx(xout = xout_end, x = Morgan_end, y = pos+5e3)] # throws a warning due to windows of 0 Morgans, but result is as desired
end[, x := x - 0.5*2^-15]
setnames(end, c("x", "y"), c("Morgan_midpoint", "wend"))

windows <- merge(start, end, by = "Morgan_midpoint")
windows <- windows[!is.na(wstart) & !is.na(wend)]


# ===== functions to get windowed frequency of archaic segments

frac_overlap <- function(x.start, x.end, y.start, y.end){
  # returns fraction of overlap of window covered by an archaic block
  x.len <- abs(x.end - x.start)
  max.start = max(c(x.start, y.start))
  min.end = min(c(x.end, y.end))
  overlap = min.end - max.start
  overlap = ifelse(overlap <= 0, 0, overlap)
  return(overlap / x.len)
}

window_freq <- function(wstart, wend){
  # gets overlapping archaic blocks, returns average over blocks of frequency weighted by overlap
  fragments <- archaic[(start > wstart & start < wend) | (end > wstart & end < wend)]
  if(nrow(fragments) == 0){
    return(0)
  }
  return(fragments[, freq*frac_overlap(wstart, wend, start, end), by = seq_len(nrow(fragments))][, sum(V1)])
}

frq <- windows[, window_freq(wstart, wend), by = wstart]
setnames(frq, "V1", "freq")
frq[, chr := chromosome]

frq <- merge(frq, windows, by = "wstart")

fwrite(frq, file = "Neanderthal_files/test.txt", quote = F, sep = "\t")
save(frq, file = paste0("chr", chromosome, "_freq.RData"))

#ggplot(frq, aes(x = wstart, y = freq)) + geom_point()
