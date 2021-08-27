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

# Expected wavelet variance: exact  integrand ------------------------------------------

wav_var_exact <- function(x, expected.crossovers.per.unit.dist, n.pop, n.sample, alpha) {
  u <- expected.crossovers.per.unit.dist
  
  v <- n.pop*u*abs(x[2]-x[1])
  w <- exp(-(t/n.pop)*(1+v))
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

# ----- Calculate expected wavelet variance for equilibrium population -----

genlist <- list()
gen <- c("0001","0002","0003","0004","0005","0010","0025","0050","0100","0500","1000")

# loop over generations 
for(i in 1:length(gen)){
  t <- as.numeric(gen[i])
  
  # loop over different population size, sample size, scale
  
  #n.sample <- 20000
  #n.pop <- 20000
  
  #n.sample <- 2*c(0.5, 10, 100, 1000, 10000, 100000)
  #n.pop <- 2*c(100, 1000, 10000, 100000, Inf)
  
  n.pop <- 2000
  n.sample <- c(1,2000)
  
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
    
    if(n.pop == Inf){ # use infinite population approximation
      part1 <- adaptIntegrate(wav_var_approx, n.sample = n.sample, expected.crossovers.per.unit.dist=1/1024, alpha=0.5, lowerLimit = c(0,0), 
                              upperLimit = c(2^(j-1),2^(j-1)))
      part2 <- adaptIntegrate(wav_var_approx, n.sample = n.sample, expected.crossovers.per.unit.dist=1/1024, alpha=0.5, lowerLimit = c(0,2^(j-1)),
                              upperLimit = c(2^(j-1),(2^j)))
    } else { # use exact formula
      part1 <- adaptIntegrate(wav_var_exact, n.sample = n.sample, n.pop = n.pop, expected.crossovers.per.unit.dist=1/1024, alpha=0.5, lowerLimit = c(0,0), 
                              upperLimit = c(2^(j-1),2^(j-1)))
      part2 <- adaptIntegrate(wav_var_exact, n.sample = n.sample, n.pop = n.pop, expected.crossovers.per.unit.dist=1/1024, alpha=0.5, lowerLimit = c(0,2^(j-1)),
                              upperLimit = c(2^(j-1),(2^j)))
    }
    grd$variance[q] <- ((part1$integral - part2$integral)/(2^(2*j-1)))
  }
  genlist[[i]] <- grd
}

dfg <- setDT(do.call(rbind.data.frame, genlist))


dfg.long <- dcast(dfg, scale + gen ~ n.sample, value.var= "variance")
setnames(dfg.long, c("1", "2000"), c("n1","n2000"))
dfg.long[, d := n1-n2000]

dfg.long[, .(varInd = sum(d)), by = gen] %>% ggplot(aes(x = gen, y = varInd)) +
  geom_line(size=2) + labs(y = "Variance of individuals' mean ancestry") + theme_classic()

dfg.long[, scale := as.factor(scale)]
dfg.long %>% ggplot(aes(x = log10(as.numeric(gen)), y = d, group = scale, color= scale)) + 
  geom_line(size=2) + 
  scale_color_viridis_d(direction = -1) + 
  theme_classic() +
  labs(color = "Level", 
       y = "Contribution to variance of individual mean ancestry", x = "Log 10 (Generation)")

# examine variance across individuals
dcast(dfg())


setDT(dfg)
ggplot(dfg, aes(x = scale, y  = variance)) + geom_point() + geom_line() + facet_wrap(~gen) + theme(aspect.ratio=1)


# ----- Read and format simulation data -----
file1 <- "theory_and_simulations/results/equilibrium/ancestry_master.txt"
file2 <- "theory_and_simulations/results/bottleneck/ancestry_master.txt"
file3 <- "theory_and_simulations/results/add-sel-const-recomb/ancestry_master.txt"
file4 <- "theory_and_simulations/results/add-sel-periodic-recomb/ancestry_master.txt"


a1 <- read.table(file1, row.names = 1)
rownames(a1) <- paste0("wvEquil_", rownames(a1))

a2 <- read.table(file2, row.names = 1)
rownames(a2) <- paste0("wvBottleneck_", rownames(a2))

a3 <- read.table(file3, row.names = 1)
rownames(a3) <- paste0("wvSelConstantRec_", rownames(a3))

a4  <- read.table(file4, row.names = 1)
rownames(a4) <- paste0("wvSelPeriodicRec_", rownames(a4))

d <- as.data.table(cbind(t(a1), t(a2), t(a3), t(a4)))
d$position <- 1:1024

# add genetic distance corresponding to variable recombination rate (expected number of crossovers, binomial n*p)
d$pos_gen <- cumsum(rep(1/1024,1024))

# tidy data
library(tidyr)
b  <-  d %>%
  gather(key = sim_rep_gen,
         value = freq, -c(position,pos_gen)) %>% 
  separate(sim_rep_gen, c("sim_rep.id", "gen"), sep = "_gen") %>% 
  separate(sim_rep.id, c("popModel", "rep.id"), sep = "_replicate")
setDT(b)

# compute wavelet variance from modwt from simulated data
wv <- function(x){
  w <- wave.variance(brick.wall(modwt(x[,freq],
                                      "haar", n.levels = 10), "haar")
  )
  w <- as.data.table(w[1:10,])
  w[, scale := 1:10]
}

wvSim <- b[,wv(.SD), by = .(popModel,rep.id, gen)]  

wvSim[, variance := wavevar]
wvSim[, gen := as.double(gen)]
wvSim[, type := "sim"]
dfg[, type := "theory1"]
dfg[,popModel := "wvEquil"]

plotData <- rbind(wvSim[popModel=="wvEquil",.(gen,scale,variance,type)], dfg[,.(gen,scale,variance,type)])

ggplot(plotData, aes(x = scale, y = variance,color = type)) + facet_wrap(~gen) + geom_point()

wvSim[, propVar := wavevar/sum(wavevar), by = .(gen,popModel, rep.id)]

ggplot(wvSim[gen %in% c(10,100,1000) & popModel != "wvBottleneck"], aes(x = scale, y = variance, color = popModel)) + 
  geom_point() + geom_line(aes(group=interaction(popModel, rep.id))) + facet_wrap(~gen) + 
  theme_classic() + 
  theme(aspect.ratio=1,
        text=element_text(size=15),
        legend.position = "none",
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5)) +
  labs(x = expression(Scale: log[2](Morgans)), y = "Variance") + 
  scale_x_continuous(breaks = 1:10, labels = -10:-1) +
  scale_color_manual(values = c("darkcyan", "gray","lightsalmon")) 
 
# Perform Computations ------------------------------------------

# make grid of parameters over which we evaluate the function
ns <- c(20000) # number of sampled haplotypes
sc <- 1:10 # scales
gn <- c(10,100,1000) # generations 
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
     value.name = "variance")
grdP[, type := "theory2"]
grdP[, propVar := variance/sum(variance), by = .(gen,popModel)]

grdP[gen %in% c(10,100,1000)] %>% 
  ggplot(aes(x = scale, y = variance, group = popModel, color = popModel)) + #interaction(popModel, n.sample))) + 
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



plotData <- rbind(plotData,grdP[popModel=="wvEquil", .(gen,scale,variance,type)])
ggplot(plotData, aes(x = scale, y = variance,color = type)) + facet_wrap(~gen) + geom_point()





# ----- combine simulation and theory ----
head(wvSim)
head(grdP)

wvSim[, type := "sim"]
grdP[, type := "theory"]
wvSim[,gen:= as.double(gen)]
allPlotDat <- merge(wvSim, grdP, by = c("gen", "popModel", "scale"))

setDT(dfg)

wvSim[popModel %in% c("wvEquil", "wvBottleneck") & gen %in% c(10,100,1000)] %>%
  ggplot(aes(x = scale, y = variance)) + #interaction(popModel, n.sample))) + 
  geom_point(size=1, alpha = 0.3, aes(color = popModel)) + # aes(shape = popModel)) + 
  geom_line(linetype = 1, aes(group = interaction(rep.id, popModel), color = popModel), alpha = 0.3) + #geom_line(aes(linetype = as.factor(n.sample))) +
  facet_wrap(~gen) + 
  theme_classic() + 
  labs(x = expression(Scale: log[2](Morgans)), 
       y = "Variance", shape = "",
       color = "") + #, linetype = "") +
  scale_color_manual(values = c("lightsalmon","darkcyan")) +
  scale_x_continuous(breaks = 1:10, labels = -10:-1) +
  #scale_linetype_discrete(labels = c("n=20")) +
  theme(aspect.ratio=1,
        text=element_text(size=15),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.5)) +
  #geom_point(data = dfg, aes(x = scale, y=var), size = 2) +
  geom_point(data = grdP, aes(x = scale, y = variance, color = popModel, group = popModel), size = 2) +
  geom_line(data = grdP, aes(x = scale,y=variance, group = popModel, color = popModel), size=1) 



