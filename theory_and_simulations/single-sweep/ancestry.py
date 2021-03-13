#!/usr/bin/env python
# coding: utf-8

import subprocess, msprime, pyslim
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os, sys

# Load the .trees file
ts = pyslim.load(sys.argv[1])

# genotypes of 30 random haplotypes
G = ts.genotype_matrix(isolated_as_missing=False)
index = np.random.choice(range(20000),30,replace=False)

#output to stdout
np.savetxt(sys.stdout.buffer, G[:,index], delimiter=" ")
