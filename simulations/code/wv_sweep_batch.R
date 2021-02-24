library(cubature)
library(data.table)
# This script computes an expected squared wavelet detail coefficient
# for a particular scale and location in the sequence, under a model
# of a single selective sweep against an introgressing allele. 

# read in generation, scale, location from command line
args <- commandArgs(trailingOnly = TRUE)
gen <- as.numeric(args[1])
j <- as.numeric(args[2])
k <- as.numeric(args[3])

# define wavlet support 
upper <- k*(2^j)
lower <- k*(2^j) - 2^j
mid <- mean(c(lower,upper))

# Continuous Haar wavelet function -----------------------------------------------
haarCts <- function(x){
  (x < lower)*0 + 
  (x >= lower & x < mid)*2^(-j/2) + 
  (x >= mid & x <= upper)*(-2^(-j/2))
}

# Expected wavelet variance with sweep ----------------------------------------------
# ignores coal. x[1] and x[2] are positions of l1, l2 (neutral loci).
# ls must be provided (we integrate over ls below)
wav_var_sweep <- function(x, j=j, r=0.01, n.sample=1, alpha=0.5, s=0.1, N=10000, t=gen, l.s) {
  
  # expected time to fixation 
  ts <- (2/s)*log( ((1-1/(2*N))*alpha ) / ( (1/(2*N)) * alpha) )
  
  #frequency of resident allele through time
  if(t > ts){
    t <- ts
    p <- 1
  } else{ 
    p <- (1-alpha)*exp(s*t/2) / ( alpha + (1-alpha)*exp(s*t/2) )
  }
  q <- 1-p
  
  # recomb. probabilities (assume large enough N to ignore coal.)
  f <- function(a,b){ exp(-t*r*abs(a-b)) }
  g <- function(a,b){1-f(a,b)}
  
  f_prime <- function(a,b){
    exp(2*r*abs(a-b)*log(alpha + (1-alpha)*exp(s*t/2))/s - r*abs(a-b)*t)
  }
  g_prime <- function(a,b){ 1 - f_prime(a,b) }
  
  u_prime <- function(a, b){
    exp(-r*abs(a-b) * 2*log(alpha + (1-alpha)*exp(s*t/2)) / s )
  }
  v_prime <- function(a,b){ 1 - u_prime(a,b)}
  
  # integrand, condition on locus configuration
  if(x[1] > x[2]){
    # order neutral loci, need to multiply by 2 later down
    cov_ii <- 0
  } else if(x[1] == x[2]){ 
    cov_ii <- alpha
  } else {
    if(x[1] >= l.s){
      # ls,l1,l2
      cov_ii <-  q*(u_prime(x[2],l.s) + alpha*(v_prime(x[1],l.s)*f(x[1],x[2]) + 
                                                 u_prime(x[1],l.s)*g(x[1],x[2]) +
                                                 alpha*v_prime(x[1],l.s)*g(x[1],x[2]))) +
        p*(alpha*g_prime(x[1],l.s)*(f(x[1],x[2]) + alpha*g(x[1],x[2])))
    } else if(x[1] < l.s && x[2] >= l.s){
      #l1,ls,l2
      cov_ii <-  q*(u_prime(x[1],l.s) + alpha*v_prime(x[1],l.s))*(u_prime(x[2],l.s) + alpha*v_prime(x[2],l.s))
      + p*alpha^2*g_prime(x[1],l.s)*g_prime(x[2],l.s)
    } else if(x[1] < l.s && x[2] < l.s){
      #l1,l2,ls
      cov_ii <-  q*(u_prime(x[1],l.s) + alpha*(v_prime(x[2],l.s)*f(x[1],x[2]) + 
                                                 u_prime(x[2],l.s)*g(x[1],x[2]) +
                                                 alpha*v_prime(x[2],l.s)*g(x[1],x[2]))) +
        p*(alpha*g_prime(x[2],l.s)*(f(x[1],x[2]) + alpha*g(x[1],x[2])))
    } else{cov_ii <- 0}
  } 
  haarCts(x[1])*haarCts(x[2])*(1/n.sample)*cov_ii
}

# Average calculation over 100 evenly spaced locations of the selected site
# (This is an approximation to the continuous integral)
vals <- NULL
for(l.s in seq(from=0,to=1024,length.out=100)){
  # avg. over location of selected site
  h <- hcubature(wav_var_sweep,c(lower,lower),c(upper,upper),
                 l.s=l.s,
                 maxEval = 1e4)
  vals <- c(vals, h$integral*2) # this is where mult by 2 comes in bc we ordered the neutral loci
}

# output gen, scale, k, expected squared wavelet coefficient
cat(gen, j, k, mean(vals))
cat("\n")


# process output 
#wvss <- fread("~/workspace/selection-against-introgression/simulations/results/single-sweep/wvss_coeffs.txt")
#setnames(wvss, c("gen", "scale", "shift", "coeff"))
#head(wvss)

#wvss <- wvss[, sum(coeff)/1024, by = .(gen,scale)]
# order scale as factor for plotting
#wvss[,genlab := paste("Gen", gen)]
#wvss$genlab <- factor(wvss$gen, levels = as.character(c(5,10,50,100,500,1000)))

# str(wvss)
# setnames(wvss, "V1", "variance")
# wvss %>% ggplot(aes(x = scale, y = variance)) + 
#   geom_point() + geom_line() + facet_wrap(~genlab) +
#   labs(x = "Scale: log2 (cM)", y = "Variance") +
#   scale_x_continuous(breaks = 1:10, labels = 1:10, minor_breaks = 1:10) +
#   theme_bw()
