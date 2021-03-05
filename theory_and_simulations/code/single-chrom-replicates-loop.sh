#!/bin/bash
# run script from top-level project directory

# modify number in tail command to run different schemes
read scheme rates ends <<< $(cat code/recomb-schemes.txt | tail -n +5 | head -n1) 
export scheme rates ends

# where to write output
export outPath=/Users/jeff/workspace/selection-against-introgression/results/single-chrom/${scheme}
mkdir -p ${outPath}


for SEED in {1..10}
do 
        slim \
            -d seed=${SEED} \
            -d "scheme='${scheme}'" \
            -d baseRate=1e-8 \
            -d rates=${rates} \
            -d ends=${ends} \
            -d "outPath='${outPath}'" \
            code/single-chrom.slim | sed -e '1,/startread/d' | paste -s - >> ${outPath}/avgs-over-time.txt
done #&&
     #   find ${outPath} -name "*.trees" | xargs -I{} python code/ancestry.py {} >> ${outPath}/ancestry-results.txt &&
     #   Rscript code/plot-ancestry-course.R ${scheme}


