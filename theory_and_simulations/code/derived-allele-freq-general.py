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
ts = pyslim.load(sys.argv[1])
#ts = tskit.load("/Users/jeff/workspace/selection-against-introgression/theory_and_simulations/results/neutral_sims/equilibrium/replicate0.trees")
rep = float(os.path.basename(sys.argv[1]).lstrip('replicate').rstrip('.trees'))

frqs = []

for gen in sys.argv[2:]: 
    tm = 1000-int(gen)
    allsets = [ts.samples(1, time = tm)]
    af = allele_frequencies(ts, allsets)
    allsites = np.array([s.position for s in ts.sites()])
    frqsout = np.c_[np.tile(rep, len(allsites)), np.tile(gen, len(allsites)), allsites, af]
    frqs.append(frqsout)

fn = os.path.abspath(sys.argv[1]).rstrip('.trees')
np.savetxt(fname=fn + "_ancestry.txt", X=np.vstack(frqs))

