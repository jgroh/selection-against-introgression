#!/bin/bash
# run script from top-level project directory

for SEED in {1..10}
do 
        slim -d seed=${SEED} code/multi-scale.slim   
done 
