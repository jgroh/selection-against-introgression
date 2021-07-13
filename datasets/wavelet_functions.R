
# Unbiased estimator of wavelet variance
wav_var <- function(x){sum(x^2,na.rm=TRUE)/(length(x[!is.na(x)]))} 

# ===== Haar dwt variance decomp for non-dyadic signals =====

haar_dwt_nondyadic_var <- function(data,variable,max.level){
  u <- data[, get(variable)]
  u <- u-mean(u) # center to mean zero
  
  # pad with zeros to next highest power of 2
  M <- length(u)
  N <- 2^(ceiling(log2(M)))
  x <- c(u, rep(mean(u), N - M)) 
  
  # haar dwt
  y <- dwt(x, "haar", n.levels = ceiling(log2(M)))
  
  # variance decomp 
  # var(u) = var(x)*N/M
  s <- lapply(y, function(v){sum(v^2)/length(u)})
  
  # the last component of s will be zero since we subtracted off the mean, so we ignore it
  s <- s[-length(s)] 
  # the remaining terms now give the wavelet spectrum.
  
  # convert to long format
  result <- melt(as.data.table(s), measure=1:length(s), variable.name="level", value.name="variance")
  result[, level := as.numeric(gsub("d","",level))]
  
  # add number of wavelets as weighting for averaging over chromosomes
  result[, n.wavelets := M/2^level]
  
  # add zero variance for levels not present up to maxLevel
  if (!max.level %in% result[,level]){
    result <- rbind(result,
                    data.table(level=setdiff(1:max.level, result[,level]),
                               variance=0,
                               n.wavelets=0))
  }
  
  return(result)
  # the sum of the variances over scales should now exactly decompose the total
  # result[, sum(variance)]
  # sum((x-mean(x))^2)/length(x)
}

# ===== Haar DWT wavelet coeffs for lm analysis =====

haar_dwt_nondyadic_coeffs <- function(x){
  # used for lm 
  x <- x-mean(x)
  M <- length(x)
  N <- 2^(ceiling(log(M, 2)))
  xx <- c(x, rep(mean(x), N - M)) # pad with zeros to next highest power of 2
  y <- dwt(xx, "haar", n.levels = ceiling(log(M, 2)))
  
  # we don't need the scaling coefficient
  y <- y[-length(y)]
  
  J <- length(y) # max level of decomp
  for(j in 1:J){
    # truncate to coefficients not affected by padding 
    # (retain coefficients partially affected)
    y[[j]] <- y[[j]][1:ceiling(M/2^j)]
    
    # set rest of coefficients to NA to equal max num of coeffs at any scale
    y[[j]] <- c(y[[j]], rep(NA, max(sapply(y, length))-length(y[[j]]) ))
  }
  
  # reformat output
  y <- as.data.table(y)
  y[, k := seq_len(.N)] # position index of wavelet
  y <- melt(y, id.vars = "k", value.name="w", variable.name="level")
  y[, level := gsub("d","",level)]
  
  return(y)
}
