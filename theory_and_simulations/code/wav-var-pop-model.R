library(data.table)
library(ggplot2)

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

# Transition probability matrix for haplotype state --------------------------------
m <- function(l_M,N_M,r_M=1/1024,tau_M){
  # l is a two vector of locus positions
  # N_M is population size
  # r_M is recombination distance per unit distance
  # tau_M is duration for which we want transition probabilities
  
  # only for l1 != l2
  # for l1==l2, haplotype state probs calculated separately in integrand funciton
  l <- as.numeric(l_M)
  L <- abs(l[2]-l[1]) 
  
  # matrix of right eigenvectors (A)
  A <- rbind(c(1,1),c(-1/(N_M*r_M*L),1))
  # matrix exponential of diagonal matrix of eigenvalues 
  Dexp <- rbind( c( exp(-tau_M*((1 + N_M*r_M*L)/N_M)), 0), c(0,1))
  return(A %*% Dexp %*% solve(A))
}


wvBottleneckIntegrand <- function(x_I,J_I,Nvec_I,genvec_I,r_I=1/1024,n.sample_I,alpha_I=0.5){
  # x is a 2-vector of neutral locus positions
  
  # Nvec is number of haploid chromosomes, a vector through time 
  # genvec are the time intervals, ordered toward the present. *These are durations of intervals.
  # Note these can be scalar values for a constant pop size.
  
  # r is recombination distance in M per unit distance of the signal
  # n.sample is the number of sampled chromosomes
  # alpha is the initial admixture proportion

  if(x_I[1] == x_I[2]){ 
    # l1==l2
    ii <- alpha_I
    
    # probability not coalesced:
    p_nc <- exp(-sum(genvec_I/Nvec_I))
    p_c <- 1 - p_nc
    ij <- alpha_I*p_c + (alpha_I^2)*p_nc
  
  } else {
    # l1! = l2: 2-state markov model for haplotype state probabilities
  
    # matrix product of transition probability matrices
    # leftmost matrix is starting from the present
    P <- m(l_M=x_I, N_M=Nvec_I[length(Nvec_I)], tau_M=genvec_I[length(genvec_I)])
  
    if(length(genvec_I) > 1){
     for(i in (length(genvec_I)-1):1){
        P <- P %*% m(l_M=x_I, N_M=Nvec_I[i], tau_M=genvec_I[i])
       }
    }

    ii_ij <- P %*% matrix(c(alpha_I,alpha_I^2), nrow=2)
    ii <- ii_ij[1,]
    ij <- ii_ij[2,]
  }
  
  return(haarCts(y_H=x_I[1],j_H=J_I)*haarCts(y_H=x_I[2],j_H=J_I)*(
      ((1/n.sample_I)*ii + ((n.sample_I-1)/n.sample_I)*ij) ) )
}


dblRiemann <- function(scale_R, nMesh=100, Nvec_R, genvec_R, n.sample_R){
  # right Riemann sum approx of integral over a square region
  pnts <- seq(0,2^scale_R,length.out=nMesh+1) 
  xy <-  expand.grid(pnts[-1], pnts[-1]) # take off 1st value for right riemann sum
  
  totalVol <- 0
  for(i in 1:nrow(xy)){
    x <- as.numeric(xy[i,])
    height <- wvBottleneckIntegrand(x, Nvec_I=Nvec_R, genvec_I=genvec_R, n.sample_I=n.sample_R, J_I=scale_R)
    totalVol <- totalVol + height*(2^scale_R/nMesh)^2
  }
  return(totalVol/2^scale_R)
}

# wvBottleneck evaluation function --------
# operates on a vector of parameters in a table
wvBottleneck <- function(d, popSizeModel,epochs){
  epochCum <- c(0,cumsum(epochs))
  n.sample <- d[,n.sample]
  scl <- d[,scale]
  g <- d[,gen] 
  
  # construct Nvec and genvec based on gen
  t0 <- g - max(epochCum[epochCum < g]) # gives remaining duration of interval spent in first epoch toward past
  if(max(epochCum[epochCum < g]) == 0){
    Tau <- t0
    N <- popSizeModel[1]
  } else{
    tn <- epochs[1:(which(epochCum == max(epochCum[epochCum < g])) - 1)] # gets duration of subsequent intervals
    Tau <- c(tn, t0) # this vec
    N <- popSizeModel[1:(which(epochCum == max(epochCum[epochCum < g])))]
  }
  
  a <- dblRiemann(scale_R=scl, nMesh=100, Nvec_R=N, genvec_R=Tau, n.sample_R=n.sample)
  return(a)
}

# Perform Computations ------------------------------------------

# make grid of parameters over which we evaluate the function
ns <- 1
sc <- 1:10
gn <- c(5,10,50,100,500,1000) 
grd1 <- data.table(expand.grid(gen=gn, n.sample=ns, scale=sc, stringsAsFactors = F))

# calculate wavelet variance under constant pop size
grd1[, wvEquil := wvBottleneck(.SD, popSizeModel=10000, epochs=1000), by = seq_len(nrow(grd1))]

# population size model of bottleneck. vector of pop size and epochs from time of admixture to present
popSizeModel <- c(10000,10,10000)
epochs <- c(50,50,900) # durations of epochs after admixture. these should go all the way to the present

# calculate wavelet variance under above model 
grd1[, wvBottleneck := wvBottleneck(.SD, popSizeModel=popSizeModel, epochs=epochs), by = seq_len(nrow(grd1))]

# Visualize results ---------------------------------
grdP <- melt(grd1, measure.vars = c("wvEquil", "wvBottleneck"),
     variable.name = "popModel", 
     value.name = "var")

ggplot(grdP, aes(x = scale, y = var, group = popModel)) + 
  geom_point(aes(shape = popModel)) + 
  geom_line(aes(linetype = popModel)) +
  facet_wrap(~gen) + 
  theme_classic() + 
  labs(x = expression(Scale: log[2](Morgans)), 
       y = "Variance", shape = "", linetype = "") +
  scale_color_manual(values = c("black","red")) +
  scale_x_continuous(breaks = 1:10, labels = -10:-1) +
  scale_linetype_discrete(labels = c("Equilibrium","Bottleneck\ngens 50-100")) +
  scale_shape_discrete(labels = c("Equilibrium","Bottleneck\ngens 50-100")) +
  theme(aspect.ratio=1,
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5)) 
  


