library(data.table)
#library(ggplot2)

chromosome <- commandArgs(trailingOnly = T)[1]

# recombination rate
rhofile <- paste0("VilgalysFogel_Amboseli_admixture-main/03_baboon-genomic-resources/Resources/Recombination_Rates/anubisSNPRC.", chromosome, ".txt.gz")
rho <- fread(rhofile, skip=3, col.names = c("start", "end", "mean", "p0.025", "p0.5", "p0.975"))

# apply threshold to outliers? see script rho_threshold.R
rho[p0.5 > 0.01, p0.5 := 0.01]
rho[mean > 0.01, mean := 0.01]
setkey(rho, start, end)

# convert to Morgan
rho[, r := p0.5/64443.48 ] # see rho_threshold.R for 2Ne estimate to convert to Morgans


# tracts
tracts <- fread("amboseli_LCLAE_tracts.txt")
tracts <- tracts[chrom == paste0("chr", chromosome)]
tracts[, anubis := state/2][, state := NULL][]
setnames(tracts, "table_s1_id", "id")

# metadata 
s1 <- fread("table_S1.txt")
s1 <- s1[, .(id, population , hybrid_status, genome_wide_anubis_ancestry, percent_reads_mapped)]

# merge tracts and metadata
d <- merge(s1[, .(id, hybrid_status, genome_wide_anubis_ancestry, population)],
      tracts[, .(id, chrom, start, end, anubis = anubis)])


# ===== ancestry in 50kb windows 

# create target grid windows along chromosome
bounds <- d[, .(max = max(end), min=min(start)), by = .(id,chrom)][, .(min=max(min), max=min(max))]
windows <- bounds[, .(start = seq(min, max, by = 50e3), end = seq(min, max, by = 50e3) + 50e3)][-.N]

#  overlaps of tracts with grid windows
setkey(d, start, end)
overlaps <- d[, foverlaps(windows, .SD, by.x=c("start","end"), by.y=c("start","end")), by = id]

# weight frequency of anubis by length of fragment overlap with target window
overlaps[, weighted_anubis := anubis*(pmin(end, i.end) - pmax(start, i.start))/50e3 , by=.(chrom, start,end)]

# what to do with windows with no overlap? doesn't appear these are in the data
# overlaps[is.na(anubis), ]

# sum weighted frequencies of fragments in each window
anubis_50kb <- overlaps[, .(anubis = sum(weighted_anubis)), by = .(id, chrom, i.start, i.end, hybrid_status)]
setnames(anubis_50kb, c("i.start", "i.end"), c("start", "end"))

# ===== recombination rate in target windows
overlaps_rho <- foverlaps(windows, rho, by.x=c("start","end"), by.y=c("start","end"))

# weight rho by length of fragment overlap with target window
#overlaps_rho[, weighted_rho_mn := mean*(pmin(end, i.end) - pmax(start, i.start))/50e3 , by=.(start,end)]
overlaps_rho[, weighted_r := r*(pmin(end, i.end) - pmax(start, i.start))/50e3 , by=.(start,end)]

rho_50kb <- overlaps_rho[, .(r=sum(weighted_r)), by = .(i.start, i.end)]

setnames(rho_50kb, c("i.start", "i.end"), c("start", "end"))

anubis_50kb <- merge(rho_50kb, anubis_50kb)
setkey(anubis_50kb, id, chrom, start, end)

# for chrom 1, validate against fig 3b from paper
#anubis_50kb[id=="AMB_202", plot(anubis ~ start)]
#anubis_50kb[id=="AMB_044", plot(anubis ~ start)]


# ===== anubis ancestry in genetic windows
# create windows of length 2^-12 Morgans
Morgan_vec <- rho[, cumsum(rep(r, end-start))]
bp_vec <- rho[, seq(min(start), max(end)-1)]

windowsM <- data.table(start_Morgan = seq(min(Morgan_vec), max(Morgan_vec), by = 2^-12))
windowsM[, end_Morgan := start_Morgan + 2^-12]
windowsM[, start_bp := round(approx(x = Morgan_vec, y = bp_vec, xout = start_Morgan)$y)]
windowsM[, end_bp := shift(start_bp, n = -1)]
windowsM <- windowsM[!is.na(end_bp)]
windowsM <- windowsM[start_bp > bounds[, min] & end_bp < bounds[, max]] # don't extend past beyond calls region
setkey(windowsM, start_bp, end_bp)
#windowsM[, mean(end_bp-start_bp)]

# overlapping anubis ancestry tracts
overlapsM <- d[, foverlaps(windowsM, .SD, by.x=c("start_bp","end_bp"), by.y=c("start","end")), by = id]

# weight frequency of anubis by length of fragment overlap with target window
overlapsM[, weighted_anubis := anubis*(pmin(end, end_bp) - pmax(start, start_bp))/(end_bp-start_bp), by=.(chrom, start, end)]

# sum weighted frequencies of fragments in each window
anubis_frqM <- overlapsM[, .(anubis = sum(weighted_anubis)), by = .(id, chrom, start_Morgan, end_Morgan, start_bp, end_bp, hybrid_status)]
setnames(anubis_frqM, c("start_bp", "end_bp"), c("start", "end"))

anubis_frqM[, r := (end_Morgan - start_Morgan)/(end-start)]
#anubis_frqM[id=="AMB_202", plot(anubis ~ start)]
#anubis_frqM[id=="AMB_044", plot(anubis ~ start)]


# ===== output

fwrite(anubis_50kb, file = paste0("anubis_freqs/chr", chromosome, "_50kb_windows.txt.gz"), 
       col.names = T, row.names = F, quote = F, sep = "\t")
fwrite(anubis_frqM, file = paste0("anubis_freqs/chr", chromosome, "_genetic_windows.txt.gz"), 
       col.names = T, row.names = F, quote = F, sep = "\t")



