library(data.table)
library(wCorr)

if(interactive()){
  setwd("~/workspace/selection-against-introgression/datasets/baboon/")
  source("~/workspace/gnomwav/R/multi_modwts.R")
  source("~/workspace/gnomwav/R/variance_decomp.R")
  source("~/workspace/gnomwav/R/correlation_decomp.R")
} else{
  source("/Users/brogroh/gnomwav/R/multi_modwts.R")
  source("/Users/brogroh/gnomwav/R/variance_decomp.R")
  source("/Users/brogroh/gnomwav/R/correlation_decomp.R")
}


# ===== read data ====
dem <- fread("VilgalysFogel_Amboseli_admixture-main/07_landscape-of-introgression/07.2.anubis_ancestry_over_time/amboseli.demographic.info.txt")
setnames(dem, "table_s1_id", "id")

# read 50kb window anubis freqs
gnomP <- rbindlist(lapply(list.files("anubis_freqs/", pattern = "50kb", full.names=T), fread))
gnomP <- merge(gnomP, dem, by = "id")

# read genetic window anubis freqs
gnomG <- rbindlist(lapply(list.files("anubis_freqs/", pattern = "genetic", full.names=T), fread))
gnomG <- merge(gnomG, dem, by = "id")

if(interactive()){
  gnomP <- gnomP[id %in% c("AMB_001", "AMB_202", "AMB_044", "AMB_002", "AMB_003") & chrom %in% paste0("chr", 15:20)]
  gnomG <- gnomG[id %in% c("AMB_001", "AMB_202", "AMB_044") & chrom %in% paste0("chr", 15:20)]
}

# ===== Wavelet Variances =====

# ----- 50 kb windows
wv_ind_P <- gnomP[, gnom_var_decomp(.SD, chromosome = "chrom", signals = "anubis"), by = .(id, hybrid_status, genome_wide_anubis_ancestry)]

# for recombination, only need to do once
wv_ind_P <- merge(wv_ind_P, gnomP[id == id[1], gnom_var_decomp(.SD, chromosome = "chrom", signals = "r")])
setkey(wv_ind_P, id)
wv_ind_P[, units := '50kb']



# ----- genetic windows
wv_ind_G <- gnomG[, gnom_var_decomp(.SD, chromosome = "chrom", signals = "anubis"), by = .(id, hybrid_status, genome_wide_anubis_ancestry)]

# for recombination, only need to do once
wv_ind_G <- merge(wv_ind_G, gnomP[id == id[1], gnom_var_decomp(.SD, chromosome = "chrom", signals = "r")])
setkey(wv_ind_G, id)
wv_ind_G[, units := 'genetic']

wv_ind_all <- rbind(wv_ind_P, wv_ind_G)


# ===== calculate mean anubis frequencies ===== 


# -- average frq in all individuals
anubis_frq_all <- gnomP[, .(anubis_frq = mean(anubis),
                        r = mean(r)), 
                        by = .(chrom, start, end)]
anubis_frq_all[, group := 'all']

# -- by hybrid status
anubis_frq_hyb_status <- gnomP[, .(anubis_frq = mean(anubis),
                              r = mean(r)), 
                          by = .(chrom, start, end, hybrid_status)]
setnames(anubis_frq_hyb_status, "hybrid_status", "group")

# -- by anubis quintile
gnomP[, anubis_qntl := cut(genome_wide_anubis_ancestry,
                           breaks = quantile(genome_wide_anubis_ancestry, probs = 0:5/5),
                           labels = FALSE, include.lowest = TRUE) ]
anubis_frq_anubis_qntl <- gnomP[, .(anubis_frq = mean(anubis),
                                      r = mean(r)), 
                                  by = .(chrom, start, end, anubis_qntl)]
anubis_frq_anubis_qntl[, anubis_qntl := paste0("qntl", anubis_qntl)]
setnames(anubis_frq_anubis_qntl, "anubis_qntl", "group")

# combine
anubis_frqs_by_grp <- rbind(anubis_frq_all, anubis_frq_hyb_status, anubis_frq_anubis_qntl)


# ===== wavelet correlations =====
cor_frq_rec <- anubis_frqs_by_grp[, gnom_cor_decomp(.SD, chromosome = "chrom",
                                                    signals = c("anubis_frq", "r")), by = group]
cor_frq_rec[, units := "50kb"]

# ----- wavelet correlation anubis freq, B values
load("VilgalysFogel_Amboseli_admixture-main/VilgalysFogel_main_data_file.250kb_windows.RData")
baboon250 <- setDT(to_analyze)

cor_frq_B <- baboon250[, gnom_cor_decomp(.SD, chromosome = "chr", signals = c("mean_ancestry", "B"))]


save(wv_ind_all, cor_frq_rec, cor_frq_B, file = "baboon_wavelet_results.RData")

