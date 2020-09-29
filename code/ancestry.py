import msprime, pyslim
import numpy as np
import pandas as pd
import sys
import os

# Load the .trees file
ts = pyslim.load(sys.argv[1]).simplify()

#breaks = np.zeros(ts.num_trees + 1)
#ancestry = np.zeros(ts.num_trees + 1)
ancestry_all_seq = np.zeros(int(ts.sequence_length))
for tree in ts.trees():
    subpop_sum, subpop_weights = 0, 0
    for root in tree.roots:
        leaves_count = tree.num_samples(root) - 1 # the root is a sample 
        subpop_sum += tree.population(root) * leaves_count 
        # for population 0, this will be 0. for population 1, it will sum the number of individuals descended from that population. 
        subpop_weights += leaves_count
#    breaks[tree.index] = tree.interval[0]
#    ancestry[tree.index] = subpop_sum / subpop_weights
    
    ancestry_all_seq[ list(range (int(tree.interval[0]), int(tree.interval[1]) ))  ] = subpop_sum / subpop_weights
    
#breaks[-1] = ts.sequence_length
#ancestry[-1] = ancestry[-2]
print(os.path.basename(sys.argv[1]).rstrip('.trees'), *ancestry_all_seq)


