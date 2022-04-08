import pyslim, tskit, msprime
import numpy as np
import os
import sys

orig_ts = tskit.load(sys.argv[1])

# recapitation
rts = pyslim.recapitate(orig_ts,
            recombination_rate=1e-8,
            ancestral_Ne=10000, random_seed=5)

# overlay mutations
ts = pyslim.SlimTreeSequence(
       msprime.sim_mutations(
           rts,
           rate=1e-8,
           model=msprime.SLiMMutationModel(type=0),
           keep=True,
           )
    )

def allele_frequencies(ts, sample_sets=None):
    if sample_sets is None:
       sample_sets = [ts.samples()] 
    n = np.array([len(x) for x in sample_sets])
    def f(x):
       return x / n
    return ts.sample_count_stat(sample_sets, f, len(sample_sets), windows='sites', polarised=True 
                                , mode='site', strict=False)


ssets_gen1 = [ts.samples(p, time = 1000) for p in range(ts.num_populations)]
af_gen1 = allele_frequencies(ts, ssets_gen1[0:3])

ssets_gen10 = [ts.samples(p, time = 990) for p in range(ts.num_populations)]
af_gen10 = allele_frequencies(ts, ssets_gen10[0:3])

ssets_gen100 = [ts.samples(p, time = 900) for p in range(ts.num_populations)]
af_gen100 = allele_frequencies(ts, ssets_gen100[0:3])

ssets_gen1000 = [ts.samples(p, time = 0) for p in range(ts.num_populations)]
af_gen1000 = allele_frequencies(ts, ssets_gen1000[0:3])


sites = np.array([s.position for s in ts.sites()])

#https://stackoverflow.com/questions/8486294/how-to-add-an-extra-column-to-a-numpy-array
gen1 = np.c_[np.tile(1, len(sites)), sites, af_gen1]
gen10 = np.c_[np.tile(10, len(sites)), sites, af_gen10]
gen100 = np.c_[np.tile(100, len(sites)), sites, af_gen100]
gen1000 = np.c_[np.tile(1000, len(sites)), sites, af_gen1000]

output = np.vstack((gen1, gen10, gen100, gen1000))

[print(*line) for line in output]
