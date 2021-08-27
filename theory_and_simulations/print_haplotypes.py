import msprime, pyslim
import numpy as np
import pandas as pd
import sys

ts = pyslim.load(sys.argv[1]).simplify()

G = ts.genotype_matrix(isolated_as_missing=False) # isolated as missing = False means that the genotypes with the ancestral state (having not had any mutations above them in the tree) will be assigned zero rather than -1. 


for h in np.arange(0, np.shape(G)[1]):
    hap = G[:,h]
    print(*hap, sep=", ")


