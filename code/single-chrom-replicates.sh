#!/bin/bash
# run script from top-level project directory

# modify number in tail command to run different schemes
read scheme rates ends <<< $(cat code/recomb-schemes.txt | tail -n +1 | head -n1) 
export scheme rates ends

# where to write output
export outPath=/Users/jeff/workspace/selection-against-introgression/results/single-chrom/${scheme}
mkdir -p ${outPath}

# modify number after seq to run different number of replicates
# modify number after parallel to allow for x number of jobs to run at once
seq 2 | parallel -q -j 2 \
        slim \
        -l \
            -d seed={} \
            -d "scheme='${scheme}'" \
            -d baseRate=1e-8 \
            -d rates=${rates} \
            -d ends=${ends} \
            -d "outPath='${outPath}'" \
            code/single-chrom.slim

