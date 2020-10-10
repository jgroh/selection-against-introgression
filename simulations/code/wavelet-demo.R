library(wavethresh)

# create a signal with two levels of oscillation

a <- rep(c(rep(10,4),rep(5,4)),16)
b <- rep(c(5,5,0,0),32)

a <- rep(c(10,10,15,15,1,1,5,5), 16)
lines(rep(c(10,10,15,15,1,1,5,5), 16))

# two levels of variation within chromosomes

# this vector times the baserate gives for two chromosomes
baseRate <- 1e-8

twoChroms <- c(c(baseRate*c(0.1,0.1,1,1)*c(10,1,10,5),0.5), c(baseRate*c(0.01,0.01,0.1,0.1)*c(10,1,10,5),0.5))
allChroms <- rep(twoChroms,4)
allChroms <- allChroms[-length(allChroms)]

plot(y = log(allChroms), x = breaks)
lines(y = log(allChroms), x = breaks)

breaks <-  c((1:4*32)-1,128)
for(i in 1:7){
 breaks <-  c(breaks, c((1:4*32)-1,128)+128*i)
}
breaks <- breaks[-length(breaks)]








a <- rep(c(rep(10,4),rep(5,4)),16)
b <- rep(c(5,5,0,0),32)


a-b
signal1 <- a-b + rnorm(128)
plot(y = signal1, x = 1:128, type = "n")
lines(y = signal1, x = 1:128, col = "red", lwd = 2)

# create another signal with added noise
signal2 <- signal1 +  rnorm(128)
lines(signal2, col = "blue", add = T, lwd = 2)

# compute wavelet transform using Haar wavelet
w <- wd(signal1, family = "DaubExPhase", filter.number = 1)
w2 <- wd(signal2, family= "DaubExPhase", filter.number = 1)

plot(w)
plot(w2)

# examine correlation at two scale levels with signal
# higher 'level' means finer scale

s1l4 <- accessD(w, level = 4) 
s2l4 <- accessD(w2, level = 4)
s1l5 <- accessD(w, level = 5)
s2l5 <- accessD(w2, level = 5)

plot(s1l4 ~ s2l4)
plot(s1l5 ~ s2l5)
# the correlation is stronger at the broader scale 
#because the added noise has comparatively less effect


