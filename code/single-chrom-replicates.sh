#!/bin/bash
# run script from top-level project directory

# modify number in tail command to run different schemes
read scheme rates ends <<< $(cat code/recomb-schemes.txt | tail -n +3 | head -n1) 
export scheme rates ends

# where to write output
export outPath=/Users/jeff/workspace/selection-against-introgression/results/single-chrom/${scheme}
mkdir -p ${outPath}

# modify number after seq to run different number of replicates
# modify number after parallel to allow for x number of jobs to run at once if multiple corers are available
seq 1 10 | parallel -q -j 1 \
        slim \
        -l \
            -d seed={} \
            -d "scheme='${scheme}'" \
            -d baseRate=1e-8 \
            -d rates=${rates} \
            -d ends=${ends} \
            -d "outPath='${outPath}'" \
            code/single-chrom.slim &&
            find ${outPath} -name "*.trees" | xargs -I{} python code/ancestry.py {} >> ${outPath}/ancestry-results.txt &&
            Rscript code/plot-ancestry-course.R ${scheme}


