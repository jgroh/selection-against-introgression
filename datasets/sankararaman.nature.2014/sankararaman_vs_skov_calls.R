# this script is to compare frequency estimates of sankararaman and skov

library(data.table)

if(Sys.getenv("RSTUDIO") == "1"){
  chromosome <- 22
  skov_file <- "Neanderthal_files/41586_2020_2225_MOESM3_ESM.txt"
  hg19file <- paste0("summaries.release/EUR.hapmap/summaries/chr-", chromosome, ".thresh-90.length-0.00.gz")
  snp_file <- "hg19ToHg38_snps.bed"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  chromosome <- args[1]jkkkk
  snp_file <- args[2] 
  hg19file <- args[3]
  skov_file <- args[4]
}

# ------ sankararaman SNP positions mapped onto hg38 coordinates:
snps <- fread(snp_file, col.names = c("chr_hg38", "start_hg38", "end_hg38", "id_hg19"))
snps[, chr_hg38 := gsub("chr", "", chr_hg38)]
snps[, start_hg38 := NULL]
setnames(snps, "end_hg38", "pos_hg38") # position of SNP now 1-based

# subset to focal chrom
snps_sub <- snps[chr_hg38 == chromosome]

# ------- skov archaic fragments
archaic_all <- fread(skov_file)

# subset to focal chrom
archaic_sub <- archaic_all[chrom == chromosome]

# ------ sankararaman calls in hg19
hg19calls <- fread(hg19file)
setnames(hg19calls, paste0("V", 1:18))
setnames(hg19calls, c("V1", "V11"), c("id_hg19", "frq_hg19"))

# merge calls with new coordinates
hg19calls[, id_hg19 := gsub(":","_",id_hg19)]
hg19calls[, id_hg19 := paste0("chr", id_hg19)]

snps_sub <- merge(snps_sub, hg19calls[, .(id_hg19, frq_hg19)], by = 'id_hg19')
setkey(snps_sub, pos_hg38)

# get skov frequencies at sankararaman snps
hg38frqs <- snps_sub[, frq_hg38 := archaic_sub[start < pos_hg38 & end >= pos_hg38, sum(freq)], by = seq_len(nrow(snps_sub))][]
hg38frqs[, frq_hg38 := frq_hg38/55132]

fwrite(hg38frqs, file = paste0('chr', chromosome, "_sankararaman_vs_skov_frqs.txt"))
