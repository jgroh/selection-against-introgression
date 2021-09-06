# Expected wavelet variance: approximate integrand ------------------------------------------
# assumes infinite population
wav_var_approx <- function(x, u, n.sample, alpha, t.gens) {
  (1/n.sample)*(
    alpha*exp(-t.gens*u*abs(x[2]-x[1])) + 
      (alpha^2)*(1-exp(-t.gens*u*abs(x[2]-x[1])))
  ) + ((n.sample-1)/n.sample)*alpha^2
}

# Expected wavelet variance: exact  integrand ------------------------------------------

wav_var_exact <- function(x, expected.crossovers.per.unit.dist,
                          n.pop, n.sample, alpha, t.gens) {
  u <- expected.crossovers.per.unit.dist
  
  v <- n.pop*u*abs(x[2]-x[1])
  w <- exp(-(t.gens/n.pop)*(1+v))
  (
    (1/n.sample)*(
      alpha*(1+v*w)/(1+v) + alpha^2*(v*(1-w))/(1+v)
    )
    + 
      ((n.sample-1)/n.sample)*(
        alpha*(1-w)/(1+v) + alpha^2*(v+w)/(1+v)) 
  )
}

# Continuous Haar Wavelet  -----------------------------------------------
haarCts <- function(y_H,j_H){
  # define wavlet support 
  upper <- 2^j_H
  lower <- 0
  mid <- upper/2
  return((y_H < lower)*0 + 
           (y_H >= lower & y_H < mid)*2^(-j_H/2) + 
           (y_H >= mid & y_H <= upper)*(-2^(-j_H/2)))
}






# Unbiased estimator of wavelet variance
wav_var <- function(x){sum(x^2,na.rm=TRUE)/(length(x[!is.na(x)]))} 

# ===== Compute expected wavelet variance for population at equilibium under neutrality =====
wavelet_variance_equilbrium <- function(n.pop, n.sample, scale, gen, alpha){
  genlist <- list()
  
  # loop over generations 
  for(i in 1:length(gen)){
    t.gens <- as.numeric(gen[i])
    
    # make grid of parameters over which we evaluate the function
    grd <- expand.grid(n.sample=n.sample, n.pop=n.pop, scale=scale, stringsAsFactors = F)
    grd <- grd[grd$n.pop >= grd$n.sample,] # we only want evaluation where the sample is less than or equal to the population size
    
    grd$gen <- rep(t.gens, nrow(grd)) # rep since we are inside the loop for a specific generation
    
    grd$variance <- vector(length = nrow(grd)) # this is the vector we fill in the calculation
    
    for(q in 1:nrow(grd)){
      j <- grd[q,]$scale
      ns <- grd[q,]$n.sample
      np <- grd[q,]$n.pop
      
      if(n.pop == Inf){ # use infinite population approximation
        part1 <- adaptIntegrate(wav_var_approx, n.sample = ns, expected.crossovers.per.unit.dist=1/1024, alpha=alpha, t.gens = t.gens, lowerLimit = c(0,0), 
                                upperLimit = c(2^(j-1),2^(j-1)))
        part2 <- adaptIntegrate(wav_var_approx, n.sample = ns, expected.crossovers.per.unit.dist=1/1024, alpha=alpha, t.gens = t.gens, lowerLimit = c(0,2^(j-1)),
                                upperLimit = c(2^(j-1),(2^j)))
      } else { # use exact formula
        part1 <- adaptIntegrate(wav_var_exact, n.sample = ns, n.pop = np, expected.crossovers.per.unit.dist=1/1024, alpha=alpha, t.gens = t.gens, lowerLimit = c(0,0), 
                                upperLimit = c(2^(j-1),2^(j-1)))
        part2 <- adaptIntegrate(wav_var_exact, n.sample = ns, n.pop = np, expected.crossovers.per.unit.dist=1/1024, alpha=alpha, t.gens = t.gens, lowerLimit = c(0,2^(j-1)),
                                upperLimit = c(2^(j-1),(2^j)))
      }
      grd$variance[q] <- ((part1$integral - part2$integral)/(2^(2*j-1)))
    }
    
    genlist[[i]] <- grd
  }
  
  return(data.table(do.call(rbind.data.frame, genlist)))
}










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


# ===== Haar MODWT variance estimation w/ brick wall boundary condition =====

haar_modwt_var <- function(data,variable,max.levels){
  u <- data[, get(variable)]
  u <- u-mean(u) # center to mean zero
  
  # wavelet variance 
  m <- brick.wall(
    modwt(u, "haar", n.levels = floor(log2(length(u)))),
    "haar")
  
  n.coeffs <- sapply(m, function(x)length(x[!is.na(x)]))
  
  wv <- data.table(wave.variance(m)[1], 
    keep.rownames = T)
  setnames(wv, "rn", "level")
  
  wv[,n.coeffs := n.coeffs]
  
  #n coeffs
  
  
  if(floor(log2(length(u))) == max(max.levels)){
    
    result <- rbind(wv, 
                    data.table(level = paste0("s", max.levels[max.levels != max(max.levels)]),
                               variance = NA, 
                               n.coeffs = NA))
  }
}
  
  
haar_dwt_coeffs <- function(data,variable){
  
  u <- data[, get(variable)]
  # haar dwt
  y <- dwt(u, "haar", n.levels = log2(length(u)))
  
  # we don't need the scaling coefficient
  y <- y[-length(y)]
  
  J <- length(y) # max level of decomp
  for(j in 1:J){
    # set rest of coefficients to NA to equal max num of coeffs at any scale
    y[[j]] <- c(y[[j]], rep(NA, max(sapply(y, length))-length(y[[j]]) ))
  }
  
  y <- as.data.table(y)
  y[, k := seq_len(.N)] # position index of wavelet
  y <- melt(y, id.vars = "k", value.name="w", variable.name="level")
  y[, level := gsub("d","",level)]
  
  return(y)
}


# ===== Haar DWT wavelet coeffs for lm analysis =====

haar_dwt_nondyadic_coeffs <- function(data, variable, max.level){
  
  u <- data[, get(variable)]
  u <- u-mean(u) # center to mean zero
  
  # pad with zeros to next highest power of 2
  M <- length(u)
  N <- 2^(ceiling(log2(M)))
  x <- c(u, rep(mean(u), N - M)) 
  
  # haar dwt
  y <- dwt(x, "haar", n.levels = ceiling(log2(M)))
  
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
