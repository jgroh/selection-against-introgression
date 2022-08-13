import msprime, pyslim, tskit
import numpy as np
import sys
import os

# Load the .trees file

orig_ts = tskit.load(sys.argv[1]) 
rep = float(os.path.basename(sys.argv[1]).lstrip('replicate').rstrip('.trees'))
N = int(sys.argv[2]) 
k = int(sys.argv[3])
L = int(orig_ts.sequence_length)


# derived allele frequency function
def allele_frequencies(ts, sample_sets=None):
    if sample_sets is None:
        sample_sets = [ts.samples()]
    n = np.array([len(x) for x in sample_sets])
    def f(x):
        return x / n
    return ts.sample_count_stat(sample_sets, f , len(sample_sets), windows='sites', polarised=True, mode='site', strict=False)


# recapitation
rts = pyslim.recapitate(orig_ts, recombination_rate=1e-8, ancestral_Ne=N)

# overlay mutations
ts = msprime.sim_mutations(rts, rate=1e-8, model=msprime.SLiMMutationModel(type=0) )



# lists that will contain: haplotype variants, population frequencies, haplotype ancestries
haps = []
frqs = []
tru_hap_anc = []

# definitions used in loop
gens = [0,2,11,101,1001] # gen 0 accesses parental individuals in first gen of admixed pop. gen 2 is 1 generation removed from F1, gen 11 is 10 gen removed from F1, etc
L = int(ts.sequence_length)
ancestors = ts.samples(time = 1001 + N)
TableCollection = orig_ts.dump_tables()
nodes_table = TableCollection.nodes

# loop over gens to get hap variants, pop faqs, hap ancestry

for gen in gens:
    tm = 1001-gen

    # variants on haplotypes from each pop in this generation
    hapset = [np.random.choice(ts.samples(population=p, time = tm), k, replace=False) for p in [0,1,2]]
    ts2 = ts.simplify(samples=np.concatenate(hapset))
    G = ts2.genotype_matrix()
    sites = [site.position for site in ts2.sites()]
    hapsout = np.c_[np.tile(gen, len(sites)), np.tile(rep, len(sites)), sites, G]
    haps.append(hapsout)
    
    # allele frequencies in this generation
    allsets = [ts.samples(p, time = tm) for p in [0,1,2]]
    af = allele_frequencies(ts, allsets)
    allsites = np.array([s.position for s in ts.sites()])
    frqsout = np.c_[np.tile(gen, len(allsites)), np.tile(rep, len(allsites)),allsites, af]
    frqs.append(frqsout)

    # true ancestry of haplotypes
    p2samples = hapset[2]
    ancestor_table = TableCollection.link_ancestors(samples=p2samples, ancestors=ancestors)
    anc_tbl = ancestor_table[np.lexsort((ancestor_table.left, ancestor_table.child))]
    source_pop = []
    for i in np.arange(len(anc_tbl)):
        parent_node_id = anc_tbl[i].parent
        source_pop.append(nodes_table[parent_node_id].population)
    
    anc_out = np.c_[np.tile(gen, len(anc_tbl)), np.tile(rep, len(anc_tbl)), anc_tbl.child, anc_tbl.left, anc_tbl.right, source_pop]
    tru_hap_anc.append(anc_out)



fn = os.path.abspath(sys.argv[1]).rstrip('.trees')

# outputs gen, rep, sites variable in haplotype set, genotype matrix for individuals selected from pops 0, 1, 2 in order
np.savetxt(fname=fn + "_haps.txt", X=np.vstack(haps))

# outputs gen, rep, variable sites, allele frqs for pops 0, 1 (parental) and 2 (hybrid)
np.savetxt(fname=fn + "_frqs.txt", X=np.vstack(frqs))

# outputs columns gen, rep, haplotype id, left coordinate (inclusive), right coordinate (exclusive), source pop
np.savetxt(fname=fn + "_ancestry.txt", X=np.vstack(tru_hap_anc))

