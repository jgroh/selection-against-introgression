L <- 1e8
x=seq(0, L-1, by = 1e5)
signal = rep(0,length(x))
baseRate=1e-8

# 10 scales of variation in recombination.
# finest scale has period of 1.65 Mb. We'll sample every 100kb 
for(i in 27:18){
  signal = signal + cos(x*2*pi/2^i)                                                                                         
}

plot(signal, cex = .5)
lines(signal)

signal = signal - min(signal)
b = baseRate/mean(signal)
r = b*signal

plot(r, cex = 0.5); lines(r)
length(r)

length(x)
x + 1e5-1

out <- c(paste(x+1e5-1, collapse=","), paste(r, collapse=","))
writeLines(out, con = "~/workspace/selection-against-introgression/theory_and_simulations/variable_rec_map.txt")
