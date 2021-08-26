library(data.table)
library(ggplot2)
library(magrittr)
library(waveslim)
library(cubature)

# Expected wavelet variance: approximate integrand ------------------------------------------
# assumes infinite population
wav_var_approx <- function(x, u, n.sample, alpha) {
  (1/n.sample)*(
    alpha*exp(-t*u*abs(x[2]-x[1])) + 
      (alpha^2)*(1-exp(-t*u*abs(x[2]-x[1])))
  ) + ((n.sample-1)/n.sample)*alpha^2
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

# Compute expected wavelet variance on genetic scale ---------------------------------------------

genlist <- list()
gen <- c("0010","0100","1000")

# loop over generations 
for(i in 1:3){
  t <- as.numeric(gen[i])
  
  # loop over different population size, sample size, scale
  
  n.sample <- 20000
  n.pop <- Inf
  scale <- 1:10
  
  # make grid of parameters over which we evaluate the function
  grd <- expand.grid(n.sample=n.sample, n.pop=n.pop, scale=scale, stringsAsFactors = F)
  grd <- grd[grd$n.pop >= grd$n.sample,] # we only want evaluation where the sample is less than or equal to the population size
  
  grd$gen <- rep(t, nrow(grd)) # rep since we are inside the loop for a specific generation
  
  grd$variance <- vector(length = nrow(grd)) # this is the vector we fill in the calculation
  
  for(q in 1:nrow(grd)){
    j <- grd[q,]$scale
    n.sample <- grd[q,]$n.sample
    n.pop <- grd[q,]$n.pop
    
      part1 <- adaptIntegrate(wav_var_approx, n.sample = n.sample, u=1/1024, alpha=0.5, lowerLimit = c(0,0), 
                              upperLimit = c(2^(j-1),2^(j-1)))
      part2 <- adaptIntegrate(wav_var_approx, n.sample = n.sample, u=1/1024, alpha=0.5, lowerLimit = c(0,2^(j-1)),
                              upperLimit = c(2^(j-1),(2^j)))
    
    grd$variance[q] <- ((part1$integral - part2$integral)/(2^(2*j-1)))
  }
  genlist[[i]] <- grd
}

dfg <- do.call(rbind.data.frame, genlist)


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
ns <- c(20000) # number of sampled haplotypes
sc <- 1:10 # scales
gn <- c(5,10,50,100,500,1000) # generations 
#gn <- c(10,100,1000)
grd1 <- data.table(expand.grid(gen=gn, n.sample=ns, scale=sc, stringsAsFactors = F))

# calculate wavelet variance under constant pop size
grd1[, wvEquil := wvBottleneck(.SD, popSizeModel=20000, epochs=1000), by = seq_len(nrow(grd1))]

# population size model of bottleneck. vector of pop size and epochs from time of admixture to present
popSizeModel <- c(200,20000)
epochs <- c(10,990) # durations of epochs after admixture. these should go all the way to the present

# calculate wavelet variance under above model 
grd1[, wvBottleneck := wvBottleneck(.SD, popSizeModel=popSizeModel, epochs=epochs), by = seq_len(nrow(grd1))]
#grd1[, wvBottleneck := NA]

# Visualize results ---------------------------------
grdP <- melt(grd1, measure.vars = c("wvEquil", "wvBottleneck"),
     variable.name = "popModel", 
     value.name = "var")

grdP[, propVar := var/sum(var), by = .(gen,popModel)]

grdP[gen %in% c(5,10,50,100,500,1000)] %>% 
  ggplot(aes(x = scale, y = var, group = popModel, color = popModel)) + #interaction(popModel, n.sample))) + 
  geom_point(size=2)+ # aes(shape = popModel)) + 
  geom_line(linetype = 2) + #geom_line(aes(linetype = as.factor(n.sample))) +
  facet_wrap(~gen) + 
  theme_classic() + 
  labs(x = expression(Scale: log[2](Morgans)), 
       y = "Variance", shape = "",
       color = "") + #, linetype = "") +
  scale_color_manual(values = c("black","red"), labels = c("equilibrium\n2N=20000\n","bottleneck\ngen1-10: 2N=200\ngen 11-1000: 2N=20000")) +
  scale_x_continuous(breaks = 1:10, labels = -10:-1) +
  #scale_linetype_discrete(labels = c("n=20")) +
  scale_shape_manual(values = c(1,2), labels = c("Equilibrium","Bottleneck\ngens 50-100")) +
  theme(aspect.ratio=1,
        text=element_text(size=15),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5))


# ===== Simulated Data =====
file1 <- "theory_and_simulations/neutral-const-recomb2/ancestry_master.txt"

a1 <- read.table(file1, row.names = 1)
rownames(a1) <- paste0("neutral-const-recomb_", rownames(a1))

d <- as.data.table(t(a1))
d$position <- 1:1024

# add genetic distance corresponding to variable recombination rate (expected number of crossovers, binomial n*p)
d$pos_gen <- cumsum(rep(1/1024,1024))

# tidy data
b  <-  d %>%
  gather(key = sim_rep_gen,
         value = freq, -c(position,pos_gen)) %>% 
  separate(sim_rep_gen, c("sim_rep.id", "gen"), sep = "_gen") %>% 
  separate(sim_rep.id, c("sim", "rep.id"), sep = "_replicate")
setDT(b)

# compute wavelet variance from modwt from simulated data
wv <- function(x){
  w <- wave.variance(
    brick.wall(modwt(x[,freq],
                     "haar", n.levels = 10),
               wf = "haar"
               )
  )
  w <- as.data.table(w[1:10,])
  w[, scale := 1:10]

}

wvSim <- b[,wv(.SD), by = .(rep.id, gen)]  

wvSim[,propVar := wavevar/sum(wavevar), by = .(gen, rep.id)]
wvSim %>% 
  ggplot(aes(x = scale, y = propVar)) + 
  geom_point() + 
  facet_wrap(~gen)

# combine expectation and sim data

head(grdP)
head(wvSim)

grdP[, propVar := var/sum(var), by = .(gen, n.sample, popModel)]
grdP[, wavevar := var]
wvSim[,gen := as.numeric(gen)]
grdP[, type := "expectation"]
wvSim[, type := "simulation"]

avgSim <- wvSim[, lapply(.SD, mean), .SDcols = c("wavevar", "propVar"), by = .(gen, scale, type)]


cols <- c("gen", "scale", "type", "wavevar", "propVar")

setDT(dfg)

dfg[,type := "expectation"]
dfg[, wavevar := variance]
dfg[, propVar := wavevar/sum(wavevar), by = .(gen)]
#allDat <- rbind(grdP[popModel == "wvEquil", ..cols], avgSim[, ..cols])
allDat <- rbind(dfg[,..cols], avgSim[,..cols])

allDat[type == "expectation"] %>%
  ggplot(aes(x = scale, y = wavevar, group = type)) + 
  geom_point() + 
  geom_line(aes(color = type)) + 
  facet_grid(~gen)



