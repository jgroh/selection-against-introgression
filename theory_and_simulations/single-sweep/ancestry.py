#!/usr/bin/env python3
# coding: utf-8

import subprocess, msprime, pyslim
import numpy as np
import pandas as pd
import os, sys

# Load the .trees file
ts = pyslim.load(sys.argv[1])

# genotypes of 1000 random haplotypes
G = ts.genotype_matrix(isolated_as_missing=False)
index = np.random.choice(range(20000),100,replace=False)

#output to stdout
np.savetxt(sys.stdout.buffer, G[:,index], delimiter=" ")
