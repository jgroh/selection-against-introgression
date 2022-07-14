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
ends <- x + 1e5-1

out_slim <- c(paste(ends, collapse=","), paste(r, collapse=","))
out_msprime <- c(paste( c(0, ends[1:(length(ends) - 1)], ends[length(ends)] + 1), collapse = "," ), paste(c(r,0.0), collapse=","))
writeLines(out_slim, con = "variable_rec_map_slim.txt")
writeLines(out_msprime, con = "variable_rec_map_msprime.txt")
