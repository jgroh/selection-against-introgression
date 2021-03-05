import msprime, pyslim
import numpy as np
import pandas as pd
import sys
import os

# this function computes allele frequency for the derived allele at all segregating sites
def allele_frequencies(ts, sample_sets=None):
    if sample_sets is None:
        sample_sets = [ts.samples()]
    n = np.array([len(x) for x in sample_sets])
    def f(x):
        return x / n
    return ts.sample_count_stat(sample_sets, f, len(sample_sets), windows='sites', polarised=True, mode='site', strict=False)


# Load the .trees file
ts = pyslim.load(sys.argv[1]).simplify()
allset = [ts.samples(p) for p in range(ts.num_populations)]
popset = [allset[0]]

# define set of samples for calculation
allset = [ts.samples(p) for p in range(ts.num_populations)]
popset = [allset[0]]

# calculate allele frequency for segregating sites and appropriate positions
val = allele_frequencies(ts, sample_sets=popset).tolist()
val = [item for sublist in val for item in sublist]
pos = [int(s.position) for s in ts.sites()] # sites of mutations, should be in same order as reported by the allele frequency calc
freq = [0]*int(ts.sequence_length) # list that will hold allele frequencies

# fill non-zero values of freq at appropriate positions
for i in range(0, len(pos)):
    freq[pos[i]] = val[i]

# output
print(os.path.basename(sys.argv[1]).rstrip('.trees'), *freq)
