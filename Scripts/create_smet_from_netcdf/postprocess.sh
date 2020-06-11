#!/bin/bash

# Retrieve script variables
model=$1

# Concatenate yearly files into a single timeseries
g=$(ls -d output/${model}_[0-9][0-9][0-9][0-9] | tail -1)
for f in ${g}/*smet
do
        g=$(basename ${f} .smet)
        cat output/${model}_*/${g}.smet | awk '{if(/^SMET/) {header=1}; if((data==1 && header==0) || (data==0 && header==1)) {print}; if(/\[DATA\]/) {data=1; header=0}}' > output/${g}.smet
done

# Zip output .smet files
zip output/smet_forcing output/VIR*.smet
