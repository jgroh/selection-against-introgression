import msprime, pyslim, tskit
import numpy as np
import sys
import os

# Load the .trees file
ts = tskit.load(sys.argv[1])
rep = float(os.path.basename(sys.argv[1]).lstrip('replicate').rstrip('.trees'))
gens = [int(g) for g in sys.argv[2:]]

frqs = []
L = int(ts.sequence_length)
for gen in gens:
    tm = 1002 - gen #1002 is the last generation of output in the human chrom1 sims
    sset=ts.samples(0, time=tm) # in the hg1 sims the population id is p0 so 0 corresponds to this
    ts2 = ts.simplify(samples=sset, keep_input_roots=True) #keep_input_roots essential here
    ancestry_all_seq = np.zeros(L)
    for tree in ts2.trees():
        subpop_sum, subpop_weights = 0, 0
        for root in tree.roots:
           # leaves_count = tree.num_samples(root) - 1 # the root is a sample 
            leaves_count = tree.num_samples(root) # the root is NOT a sample after the simplification
            subpop_sum += tree.population(root) * leaves_count 
            # for population 0, this will be 0. for population 1, it will sum the number of individuals descended from that population. 
            subpop_weights += leaves_count
        ancestry_all_seq[ list(range (int(tree.interval[0]), int(tree.interval[1]) ))  ] = subpop_sum / subpop_weights

    frqsout = np.c_[np.tile(rep, L), np.tile(gen, L), ancestry_all_seq]
    frqsout.tolist()
    for a in frqsout:
        print(*a)
