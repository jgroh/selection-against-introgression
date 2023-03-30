import msprime, pyslim, tskit
import numpy as np
import sys
import os

def ts_get_times(ts):
    has_indiv = ts.tables.nodes.individual >=0
    which_indiv = ts.tables.nodes.individual[has_indiv]
    individual_times = np.zeros(ts.num_individuals)
    individual_times[which_indiv] = ts.tables.nodes.time[has_indiv]
    return np.unique(individual_times)

# Load the .trees file

orig_ts = tskit.load(sys.argv[1]) 
rep = float(os.path.basename(sys.argv[1]).lstrip('replicate').rstrip('.trees'))
N = int(sys.argv[2]) 
k = int(sys.argv[3])
L = int(orig_ts.sequence_length)

# list that will contain haplotype ancestries
tru_hap_anc = []

times = ts_get_times(orig_ts)

# definitions used in loops
ancestors = orig_ts.samples(time = max(times))
TableCollection = orig_ts.dump_tables()
nodes_table = TableCollection.nodes


# loop over gens to get hap ancestry

for tm in times:
    gen = max(times) - tm

    hapset = [np.random.choice(orig_ts.samples(population=p, time = tm), k, replace=False) for p in [0]]
    ts2 = orig_ts.simplify(samples=np.concatenate(hapset))
    G = ts2.genotype_matrix()
    
    # true ancestry of haplotypes
    p0samples = hapset[0]
    ancestor_table = TableCollection.link_ancestors(samples=p0samples, ancestors=ancestors)
    anc_tbl = ancestor_table[np.lexsort((ancestor_table.left, ancestor_table.child))]
    source_pop = []
    for i in np.arange(len(anc_tbl)):
        parent_node_id = anc_tbl[i].parent
        source_pop.append(nodes_table[parent_node_id].population)
    
    anc_out = np.c_[np.tile(gen, len(anc_tbl)), np.tile(rep, len(anc_tbl)), anc_tbl.child, anc_tbl.left, anc_tbl.right, source_pop]
    tru_hap_anc.append(anc_out)


fn = os.path.abspath(sys.argv[1]).rstrip('.trees')

# outputs columns gen, rep, haplotype id, left coordinate (inclusive), right coordinate (exclusive), source pop
np.savetxt(fname=fn + "_hap_ancestry.txt", X=np.vstack(tru_hap_anc))
