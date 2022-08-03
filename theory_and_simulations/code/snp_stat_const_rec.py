import pyslim, tskit, msprime
import numpy as np
import os
import sys

orig_ts = pyslim.load(sys.argv[1]) 
rep = float(os.path.basename(sys.argv[1]).lstrip('replicate').rstrip('.trees'))
N = int(sys.argv[2])
k = int(sys.argv[3])

def allele_frequencies(ts, sample_sets=None):
    if sample_sets is None:
        sample_sets = [ts.samples()]
    n = np.array([len(x) for x in sample_sets])
    def f(x):
        return x / n
    return ts.sample_count_stat(sample_sets, f , len(sample_sets), windows='sites', polarised=True, mode='site', strict=False)

# recapitation
rts = pyslim.recapitate(orig_ts, recombination_rate=1e-8, ancestral_Ne=N, random_seed=123)

# overlay mutations
ts = pyslim.SlimTreeSequence(msprime.sim_mutations(rts, rate=1e-8, model=msprime.SLiMMutationModel(type=0) ))

haps = []
frqs = []
# loop over desired generations
for gen in [0,10,100,500,1000]:
    tm = 1000-gen

    #haplotypes from each pop in this generation
    hapset = [np.random.choice(ts.samples(population=p, time = tm), k, replace=False) for p in [0,1,2]]
    ts2 = ts.simplify(samples=np.concatenate(hapset))
    G = ts2.genotype_matrix()
    sites = [site.position for site in ts2.sites()]
    hapsout = np.c_[np.tile(rep, len(sites)), np.tile(gen, len(sites)), sites, G]
    haps.append(hapsout)

    #output allele frequencies in this generation
    allsets = [ts.samples(p, time = tm) for p in [0,1,2]]
    af = allele_frequencies(ts, allsets)
    allsites = np.array([s.position for s in ts.sites()])
    frqsout = np.c_[np.tile(rep, len(allsites)), np.tile(gen, len(allsites)), allsites, af]
    frqs.append(frqsout)

fn = os.path.abspath(sys.argv[1]).rstrip('.trees')
np.savetxt(fname=fn + "_haps.txt", X=np.vstack(haps))
np.savetxt(fname=fn + "_frqs.txt", X=np.vstack(frqs))
