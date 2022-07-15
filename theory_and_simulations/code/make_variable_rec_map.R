options(scipen=999)
L <- 1e8
x=seq(0, L-1, by = 1e5)
signal = rep(0,length(x))
baseRate=1e-8

# 10 scales of variation in recombination.
# finest scale has period of 1.65 Mb. We'll sample every 100kb 
for(i in 27:18){
  signal = signal + cos(x*2*pi/2^i)                                                                                         
}

signal = signal - min(signal)
b = baseRate/mean(signal)
r = b*signal

ends <- x + 1e5 -1

slim_ends <- ends
slim_rates <- r

out_slim <- c(paste(slim_ends,collapse=','), paste(slim_rates, collapse=','))

msprime_ends <- c(0, ends[1:(length(ends)-1)], ends[length(ends)]+1)
msprime_rates <- c(slim_rates, 0) # pyslim tutorial is confusing on this
out_msprime <- c(paste(msprime_ends,collapse=','), paste(msprime_rates, collapse=','))

writeLines(out_slim, con = "variable_rec_map_slim.txt")
writeLines(out_msprime, con = "variable_rec_map_msprime.txt")
