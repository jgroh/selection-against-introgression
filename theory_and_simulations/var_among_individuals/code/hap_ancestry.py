import msprime, pyslim, tskit
import numpy as np
import sys
import os


# Load the .trees file
orig_ts = tskit.load(sys.argv[1]) 
bn = os.path.basename(sys.argv[1])
alpha = float(bn.split('_')[0].lstrip('alpha'))
rep = float(bn.split('_')[1].lstrip('replicate').rstrip('.trees'))
N = int(sys.argv[2]) 
k = int(sys.argv[3])
L = int(orig_ts.sequence_length)


# lists that will contain: haplotype ancestries
tru_hap_anc = []


# definitions used in loops
gens = [1,2,11,101,1001] # gen 1 accesses F1s in admixed pop.  gen 2 is 1 generation removed from F1, gen 11 is 10 gen removed from F1, etc

ancestors = orig_ts.samples(time = 1001 + 2) # the +2 corresponds to the fact that migration begins in gen 2
TableCollection = orig_ts.dump_tables()
nodes_table = TableCollection.nodes


# loop over gens to get hap ancestry

for gen in gens:
    tm = 1001-gen

    # variants on haplotypes from each pop in this generation
    hapset = np.random.choice(orig_ts.samples(population=2, time = tm), k, replace=False)

    # true ancestry of haplotypes
    p2samples = hapset[2]
    ancestor_table = TableCollection.link_ancestors(samples=p2samples, ancestors=ancestors)
    anc_tbl = ancestor_table[np.lexsort((ancestor_table.left, ancestor_table.child))]
    source_pop = []
    for i in np.arange(len(anc_tbl)):
        parent_node_id = anc_tbl[i].parent
        source_pop.append(nodes_table[parent_node_id].population)
    
    anc_out = np.c_[np.tile(gen, len(anc_tbl)), np.tile(alpha, len(anc_tbl)),  np.tile(rep, len(anc_tbl)), anc_tbl.child, anc_tbl.left, anc_tbl.right, source_pop]
    tru_hap_anc.append(anc_out)

fn = os.path.abspath(sys.argv[1]).rstrip('.trees')

# outputs columns gen, rep, alpha, haplotype id, left coordinate (inclusive), right coordinate (exclusive), source pop
np.savetxt(fname=fn + "_hap_ancestry.txt", X=np.vstack(tru_hap_anc))



