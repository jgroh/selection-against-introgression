# arguments to this script are (in order) file name, followed by generations of output from slim simulation
import msprime, pyslim
import tskit
import numpy as np
import pandas as pd
import sys
import os


def allele_frequencies(ts, sample_sets=None):
    if sample_sets is None:
       sample_sets = [ts.samples()] 
    n = np.array([len(x) for x in sample_sets])
    def f(x):
       return x / n
    return ts.sample_count_stat(sample_sets, f, len(sample_sets), windows='sites', polarised=True, mode='site', strict=False)


# Load the .trees file
#ts = tskit.load(sys.argv[1])
ts = tskit.load("/Users/jeff/workspace/selection-against-introgression/theory_and_simulations/results/neutral_sims/equilibrium/replicate0.trees")
#rep = float(os.path.basename(sys.argv[1]).lstrip('replicate').rstrip('.trees'))
rep = float(0)


haps = []
frqs = []

for gen in sys.argv[2:]:
    tm = 1000-int(gen)
    hapset = [np.random.choice(ts.samples(population=1, time = tm), 10, replace=False)]

    for s in np.arange(len(hapset[0])):
        sset = [np.array([hapset[0][s]])] # many parentheses bc it likes it I guess
        af1 = allele_frequencies(ts, sset) # allele frequency for single individual gives their ancestry state
        hapout = np.c_[np.tile(s, 1024), np.tile(gen, 1024), np.tile(rep, 1024), af1]
        haps.append(hapout) 
        
    popset = [ts.samples(1, time = tm)]
    afn = allele_frequencies(ts, popset)
    frqsout = np.c_[np.tile(gen, 1024), np.tile(rep, 1024), 1:1025, afn]
    frqs.append(frqsout)
    



fn = os.path.abspath(sys.argv[1]).rstrip('.trees')
np.savetxt(fname = fn + "_haps.txt", X=np.vstack(haps), fmt = "%s")
np.savetxt(fname = fn + "_frqs.txt", X=np.vstack(frqs), fmt = "%s")
