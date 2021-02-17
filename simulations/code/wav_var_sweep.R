library(cubature)
#library(tidyverse)

# Continuous Haar wavelet function -----------------------------------------------
haarCts <- function(x, j){
  (x <= 0)*0 + (x > 0 & x <= 2^(j-1))*2^(-j/2) + (x > 2^(j-1) & x <= 2^j)*(-2^(-j/2))
}

# Expected wavelet variance: single sweep ----------------------------------------------
# assumes infinite population. x[1] and x[2] are positions of l1,l2, ls assumed=0
wav_var_sweep <- function(x, j, r, n.sample, alpha, s, t, N, l.s) {
  
  #frequency of resident allele through time
  if(t > (4/s)*log(2*N)){
    p <- 1
  } else{ 
    p <- (1-alpha)*exp(s*t/2) / ( alpha + (1-alpha)*exp(s*t/2) )
  }
  q <- 1-p
  
  # recomb. probabilities (assume large enough N to ignore coal.)
  f_prime <- function(a, b){
    exp(-r*abs(a-b) * 2*log(alpha + (1-alpha)*exp(s*t/2)) / s )
  }
  g_prime <- function(a,b){ 1 - f_prime(a,b) }
  f <- function(a,b){ exp(-t*r*abs(a-b)) }
  g <- function(a,b){1-f(a,b)}
  
  # integrand, condition on locus configuration
  if(x[1] > x[2]){
    # order neutral loci, need to multiply by 2 later down
    cov_ii <- 0
  } else if(x[1] == x[2]){ 
    cov_ii <- alpha
    } else {
      if(x[1] >= l.s){
        # ls,l1,l2
        cov_ii <-  q*(f_prime(x[2],l.s) + alpha*(g_prime(x[1],l.s)*f(x[1],x[2]) + 
                                                 f_prime(x[1],l.s)*g(x[1],x[2]) +
                                                 alpha*g_prime(x[1],l.s)*g(x[1],x[2]))) +
        p*(alpha*g_prime(x[1],l.s)*(f(x[1],x[2]) + alpha*g(x[1],x[2])))
      } else if(x[1] < l.s && x[2] >= l.s){
        #l1,ls,l2
        cov_ii <-  q*(f_prime(x[1],l.s) + alpha*g_prime(x[1],l.s))*(f_prime(x[2],l.s) + alpha*g_prime(x[2],l.s))
        + p*alpha^2*g_prime(x[1],l.s)*g_prime(x[2],l.s)
      } else if(x[1] < l.s && x[2] < l.s){
        #l1,l2,ls
        cov_ii <-  q*(f_prime(x[1],l.s) + alpha*(g_prime(x[2],l.s)*f(x[1],x[2]) + 
                                                 f_prime(x[2],l.s)*g(x[1],x[2]) +
                                                 alpha*g_prime(x[2],l.s)*g(x[1],x[2]))) +
        p*(alpha*g_prime(x[2],l.s)*(f(x[1],x[2]) + alpha*g(x[1],x[2])))
      } else{cov_ii <- 0}
  } 
  haarCts(x[1], j=j)*haarCts(x[2],j=j)*(1/n.sample)*cov_ii
}

# Calculate expected wavelet variance with single sweep -------------------------

# integrate
genlist <- list()
gen <- c("0005","0010","0050","0100","0500","1000")
names(gen) <- gen
gen <- as.list(gen)

compute_all_wv <- function(x){
  t <- as.numeric(x)
  scale_var <- NA
  for(j in 1:10){
    upper <- 2^j
    vals_j <- NULL
    for(l.s in seq(from=0,to=1024,length.out=1000)){
      # avg. over location of selected site
      h <- hcubature(wav_var_sweep,c(0,0),c(upper,upper),j=j,
                   alpha=0.5,s=0.1,t=t,r=0.01,
                   n.sample=1,N=10000,l.s=l.s,
                   maxEval = 1e5)
      vals_j <- c(vals_j, h$integral*(1/2^(j-1))) # this is where mult by 2 comes in
    }
    scale_var[j] <- mean(vals_j)
  }
  return(scale_var)
}

dfs <- as.data.frame(sapply(gen,compute_all_wv,USE.NAMES = T))
dfs$scale <- 1:10
save(dfs, file = "/Users/brogroh/selection-against-introgression/simulations/results/WV_single_sweep.RData")

#dfs <- gather(dfs, key = "gen", value = "variance", -scale)
#save(dfs, file = "~/workspace/selection-against-introgression/simulations/results/wav_var_single_sweep.RData")
#dfs %>% ggplot(aes(x = scale, y=variance,group=gen,color=gen)) + 
#  geom_point() + geom_line() + labs(x = "Scale (log2 cM)", y = "Wavelet variance")


